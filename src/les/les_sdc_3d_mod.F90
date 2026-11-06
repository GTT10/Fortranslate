module les_sdc_3d_mod
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_mol_3d_mod, only: advance_reactive_mol_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use reactive_3d_mod, only: advance_reactive_sdc_3d,compute_reactive_external_source_3d
  use reactive_3d_mod, only: reactive_sdc_halo_filler,reactive_sdc_source_halo_filler,reactive_sdc_stage_agreement
  use les_filter_3d_mod, only: les_filter_weights
  use les_smagorinsky_3d_mod, only: compute_les_rhs_3d,compute_les_ghosted_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_boundary_rhs_3d,compute_dynamic_les_rhs_ghosted_3d
  implicit none
  private
  public :: advance_les_sdc_3d,les_sdc_auxiliary_source
  abstract interface
    subroutine les_sdc_auxiliary_source(state,temperature,time,interval,stage,rate,ok)
      import dp
      real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),time,interval
      integer, intent(in) :: stage
      real(dp), intent(out) :: rate(:,:,:,:)
      logical, intent(out) :: ok
    end subroutine les_sdc_auxiliary_source
  end interface
contains
  ! LES is an old/new source in the SDC predictor and correction. Its surface
  ! inventory belongs to the boundary ledger, not the volumetric forcing ledger.
  subroutine advance_les_sdc_3d(species,state,temperature,spacing,dt,time,iterations,cs,prt,sct, &
      boundary_change,source_change,maximum_diffusivity,ok,model,ci,external_force,transport,transport_flags, &
      reactions,reaction_source,chemistry_relative_tolerance,chemistry_absolute_tolerance,difmag,boundary_faces, &
      inflow_state,inflow_temperature,inflow_mach,inflow_pressure,outlet_pressure,wall_temperature,wall_velocity, &
      test_filter,test_ratio,grid_ratio,harmonic_mean,flux_x,flux_y,flux_z,fill_halos,fill_source_halos,new_source_estimate, &
      evaluate_reaction,constant_coefficients,auxiliary_source,auxiliary_change,heat_release,agree_stage, &
      mol_mode,mol_lax_friedrichs)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    ! Distributed callers agree metadata before entry; local failures are then
    ! agreed before entering another halo or auxiliary-source callback.
    procedure(reactive_sdc_stage_agreement), optional :: agree_stage
    procedure(les_sdc_auxiliary_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs
    type(eb_integrated_flux_3d) :: mol_flux
    logical :: mol,chemistry
    real(dp) :: rtol,atol
    integer :: auxiliary_families
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(reactive_sdc_halo_filler), optional :: fill_halos
    procedure(reactive_sdc_source_halo_filler), optional :: fill_source_halos
    ! Families are force, LES, optional auxiliary exchange, and molecular transport. For external
    ! iterations, reset gas/T to the interval start and retain this estimate.
    real(dp), optional, intent(inout) :: new_source_estimate(:,:,:,:,:)
    real(dp), allocatable :: estimate(:,:,:,:,:)
    ! Time-averaged full-face flux density, including hydro, transport and LES.
    real(dp), optional, intent(inout) :: flux_x(:,:,:,:),flux_y(:,:,:,:),flux_z(:,:,:,:)
    real(dp), allocatable :: x(:,:,:,:),y(:,:,:,:),z(:,:,:,:),lx(:,:,:,:),ly(:,:,:,:),lz(:,:,:,:)
    real(dp), allocatable :: old_x(:,:,:,:),old_y(:,:,:,:),old_z(:,:,:,:)
    integer :: n(3),nv
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp) :: old_auxiliary(size(state,1)),new_auxiliary(size(state,1)),auxiliary_inventory(size(state,1))
    real(dp), intent(in) :: spacing(3),dt,time,cs,prt,sct
    integer, intent(in) :: iterations
    real(dp), intent(inout) :: boundary_change(:),source_change(:),maximum_diffusivity
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: model,test_filter,boundary_faces(6)
    real(dp), optional, intent(in) :: ci,external_force(3),difmag
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: transport_flags(4),harmonic_mean
    type(elementary_reaction), optional, intent(in) :: reactions(:)
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:)
    real(dp), optional, intent(inout) :: heat_release(:,:,:)
    real(dp), optional, intent(in) :: chemistry_relative_tolerance,chemistry_absolute_tolerance
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),wall_temperature(:),wall_velocity(:,:)
    integer, optional, intent(in) :: test_ratio,grid_ratio
    real(dp), allocatable :: candidate(:,:,:,:),ct(:,:,:),ir(:,:,:,:),heat(:,:,:)
    real(dp) :: boundary(size(state,1)),source(size(state,1)),old_rate(size(state,1)),new_rate(size(state,1))
    real(dp) :: force(3),isotropic,diffusion,les_boundary(size(state,1))
    integer :: calls,ratio
    character(len=32) :: closure,filter
    logical :: valid
    ok=.false.;calls=0;diffusion=0;old_rate=0;new_rate=0
    mol=.false.;if(present(mol_mode)) mol=mol_mode
    chemistry=.false.;if(present(reactions)) chemistry=size(reactions)>0
    rtol=2e-7_dp;atol=1e-12_dp
    if(present(chemistry_relative_tolerance)) rtol=chemistry_relative_tolerance
    if(present(chemistry_absolute_tolerance)) atol=chemistry_absolute_tolerance
    old_auxiliary=0;new_auxiliary=0;auxiliary_inventory=0
    auxiliary_families=merge(1,0,present(auxiliary_source))
    call prepare(valid)
    if(.not.stage_ready(valid)) return
    if(mol) then
      call advance_reactive_mol_3d(species,candidate,ct,spacing,dt,time,iterations,chemistry,rtol,atol,valid, &
        reactions=reactions,reaction_source=ir,heat_release=heat,transport=transport,transport_flags=transport_flags, &
        source_count=2+auxiliary_families,evaluate_sources=evaluate,source_change=source,boundary_change=boundary, &
        integrated_flux=mol_flux,boundary_faces=boundary_faces,inflow_state=inflow_state, &
        inflow_temperature=inflow_temperature,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
        outlet_pressure=outlet_pressure,wall_temperature=wall_temperature,wall_velocity=wall_velocity, &
        constant_coefficients=constant_coefficients,fill_halos=fill_halos,agree_stage=agree_stage, &
        evaluate_reaction=evaluate_reaction,lax_friedrichs=mol_lax_friedrichs)
      if(valid.and.present(flux_x)) then
        x=mol_flux%x/dt;y=mol_flux%y/dt;z=mol_flux%z/dt
      endif
    else
    call advance_reactive_sdc_3d(species,candidate,ct,spacing,dt,time,iterations, &
      2+auxiliary_families+merge(1,0,present(transport)),evaluate,valid,difmag=difmag, &
      source_change=source,transport=transport,transport_flags=transport_flags,reactions=reactions,reaction_source=ir, &
      chemistry_relative_tolerance=chemistry_relative_tolerance,chemistry_absolute_tolerance=chemistry_absolute_tolerance, &
      boundary_faces=boundary_faces,boundary_change=boundary,outlet_pressure=outlet_pressure, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,wall_velocity=wall_velocity,flux_x=x,flux_y=y,flux_z=z, &
      fill_halos=fill_halos,fill_source_halos=fill_source_halos,new_source_estimate=estimate,evaluate_reaction=evaluate_reaction, &
      constant_coefficients=constant_coefficients,heat_release=heat,agree_stage=agree_stage)
    endif
    if(.not.valid) return
    if(mol.and.chemistry.and.iterations>1) then
      new_rate=0;new_auxiliary=0
      if(present(flux_x)) then
        lx=0;ly=0;lz=0
      endif
    endif
    les_boundary=.5_dp*dt*(old_rate+new_rate)
    boundary=boundary+les_boundary;source=source-les_boundary
    if(present(auxiliary_source)) then
      auxiliary_inventory=.5_dp*dt*(old_auxiliary+new_auxiliary)
      source=source-auxiliary_inventory
    endif
    valid=all(ieee_is_finite([boundary,source,diffusion,auxiliary_inventory]))
    if(present(flux_x)) then
      x=x+.5_dp*(old_x+lx);y=y+.5_dp*(old_y+ly);z=z+.5_dp*(old_z+lz)
      valid=valid.and.all(ieee_is_finite(x)).and.all(ieee_is_finite(y)).and.all(ieee_is_finite(z))
    end if
    if(.not.stage_ready(valid)) return
    if(present(flux_x)) then
      flux_x=x;flux_y=y;flux_z=z
    endif
    if(present(new_source_estimate)) new_source_estimate=estimate
    if(present(auxiliary_change)) auxiliary_change=auxiliary_inventory
    state=candidate;temperature=ct
    if(present(reaction_source)) reaction_source=ir
    if(present(heat_release)) heat_release=heat
    boundary_change=boundary;source_change=source;maximum_diffusivity=diffusion;ok=.true.
  contains
    logical function stage_ready(local_ok) result(ready)
      logical, intent(in) :: local_ok
      ready=local_ok
      if(present(agree_stage)) call agree_stage(local_ok,ready)
    end function stage_ready
    subroutine prepare(good)
      logical, intent(out) :: good
      good=.false.
      if(mol.and.present(new_source_estimate)) return
      if(present(auxiliary_source).neqv.present(auxiliary_change)) return
      if(present(auxiliary_change)) then
        if(size(auxiliary_change)/=size(state,1)) return
      endif
      if(present(fill_halos).neqv.present(fill_source_halos)) return
      if(size(boundary_change)/=size(state,1).or.size(source_change)/=size(state,1)) return
      closure='smagorinsky';filter='box';ratio=2;isotropic=0;force=0
      if(present(model)) closure=model
      if(present(test_filter)) filter=test_filter
      if(present(test_ratio)) ratio=test_ratio
      if(present(ci)) isotropic=ci
      if(present(external_force)) force=external_force
      if(.not.all(ieee_is_finite([cs,prt,sct,isotropic,force]))) return
      if(cs<0.or.cs>2.or.min(prt,sct)<1e-6_dp.or.max(prt,sct)>1e6_dp) return
      if(isotropic<0.or.isotropic>1) return
      select case(closure)
      case('smagorinsky','wale','vreman')
        if(present(harmonic_mean)) then
          if(harmonic_mean) return
        end if
      case('dynamic')
        if(isotropic/=0) return
      case default
        return
      end select
      n=shape(temperature);nv=size(state,1)
      if((present(flux_x).neqv.present(flux_y)).or.(present(flux_x).neqv.present(flux_z))) return
      if(present(flux_x)) then
        if(any(shape(flux_x)/=[nv,n(1)+1,n(2),n(3)])) return
        if(any(shape(flux_y)/=[nv,n(1),n(2)+1,n(3)])) return
        if(any(shape(flux_z)/=[nv,n(1),n(2),n(3)+1])) return
        allocate(x,mold=flux_x);allocate(y,mold=flux_y);allocate(z,mold=flux_z)
        allocate(lx,mold=x);allocate(ly,mold=y);allocate(lz,mold=z)
      end if
      if(present(fill_halos).and..not.allocated(lx)) then
        allocate(lx(nv,n(1)+1,n(2),n(3)),ly(nv,n(1),n(2)+1,n(3)),lz(nv,n(1),n(2),n(3)+1))
      end if
      if(present(new_source_estimate)) allocate(estimate,source=new_source_estimate)
      candidate=state;ct=temperature
      if(present(heat_release)) then
        if(any(shape(heat_release)/=n)) return
        allocate(heat,mold=heat_release)
      endif
      if(present(reaction_source)) then
        if(.not.present(reactions).and..not.mol) return
        allocate(ir,source=reaction_source)
      end if
      good=.true.
    end subroutine prepare
    subroutine evaluate(u,t,stage_time,rates,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time
      real(dp), intent(out) :: rates(:,:,:,:,:)
      logical, intent(out) :: good
      real(dp) :: rate(size(state,1)),viscosity
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),weights(:)
      integer :: width,c
      rates=0;good=.false.;rate=0;viscosity=0
      if(.not.stage_ready(ieee_is_finite(stage_time).and.size(rates,5)==2+auxiliary_families)) return
      call compute_reactive_external_source_3d(u,force,rates(:,:,:,:,1),good)
      if(.not.stage_ready(good)) then
        good=.false.;return
      endif
      if(present(fill_halos)) then
        width=2
        if(closure=='dynamic') then
          call les_filter_weights(filter,ratio,weights,good)
          if(.not.stage_ready(good)) then
            good=.false.;return
          endif
          width=(size(weights)-1)/2+5
        end if
        if(.not.stage_ready(width<=128)) then
          good=.false.;return
        end if
        allocate(ghost(nv,n(1)+2*width,n(2)+2*width,n(3)+2*width),gt(n(1)+2*width,n(2)+2*width,n(3)+2*width))
        call fill_halos(u,t,stage_time,width,.true.,ghost,gt,good)
        if(.not.stage_ready(good)) then
          good=.false.;return
        endif
        if(closure=='dynamic') then
          call compute_dynamic_les_rhs_ghosted_3d(species,ghost,gt,spacing,filter,ratio,width, &
            rates(:,:,:,:,2),lx,ly,lz,good,grid_ratio,harmonic_mean,viscosity,boundary_faces,wall_temperature)
        else
          call compute_les_ghosted_rhs_3d(species,ghost,gt,n,spacing,cs,prt,sct,rates(:,:,:,:,2),viscosity,good, &
            boundary_faces=boundary_faces,wall_temperature=wall_temperature,model=closure,ci=isotropic,fx=lx,fy=ly,fz=lz)
        end if
        if(.not.stage_ready(good)) then
          good=.false.;return
        endif
        do c=1,nv
          rate(c)=(sum(lx(c,1,:,:))-sum(lx(c,n(1)+1,:,:)))*spacing(2)*spacing(3)+ &
            (sum(ly(c,:,1,:))-sum(ly(c,:,n(2)+1,:)))*spacing(1)*spacing(3)+ &
            (sum(lz(c,:,:,1))-sum(lz(c,:,:,n(3)+1)))*spacing(1)*spacing(2)
        end do
      else if(closure=='dynamic') then
        call compute_dynamic_les_boundary_rhs_3d(species,u,t,spacing,filter,ratio,rates(:,:,:,:,2),viscosity,good, &
          grid_ratio=grid_ratio,harmonic_mean=harmonic_mean,boundary_faces=boundary_faces,boundary_rate=rate, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
          wall_velocity=wall_velocity,outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
          flux_x=lx,flux_y=ly,flux_z=lz)
      else
        call compute_les_rhs_3d(species,u,t,spacing,cs,prt,sct,rates(:,:,:,:,2),viscosity,good, &
          boundary_faces=boundary_faces,boundary_rate=rate,model=closure,ci=isotropic, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_temperature=wall_temperature, &
          wall_velocity=wall_velocity,outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
          fx=lx,fy=ly,fz=lz)
      end if
      if(.not.stage_ready(good)) then
        good=.false.;return
      endif
      if(present(auxiliary_source)) then
        call auxiliary_source(u,t,stage_time,dt,calls,rates(:,:,:,:,3),good)
        if(.not.stage_ready(good)) then
          good=.false.;return
        endif
        good=all(ieee_is_finite(rates(:,:,:,:,3)))
        if(.not.stage_ready(good)) then
          good=.false.;return
        endif
        do c=1,nv
          new_auxiliary(c)=sum(rates(c,:,:,:,3))*product(spacing)
        enddo
        good=all(ieee_is_finite(new_auxiliary))
        if(.not.good) return
        if(calls==0) old_auxiliary=new_auxiliary
      endif
      calls=calls+1
      if(calls==1) then
        old_rate=rate
        if(present(flux_x)) then
          allocate(old_x,source=lx);allocate(old_y,source=ly);allocate(old_z,source=lz)
        end if
      end if
      new_rate=rate;diffusion=max(diffusion,viscosity)
    end subroutine evaluate
  end subroutine advance_les_sdc_3d
end module les_sdc_3d_mod
