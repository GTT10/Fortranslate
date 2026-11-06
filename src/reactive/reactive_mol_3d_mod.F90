module reactive_mol_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use gas_transport_mod, only: gas_transport_species
  use reactive_transport_3d_mod, only: physical_transport_fluxes_3d,physical_transport_ghosted_fluxes_3d, &
    valid_constant_transport_coefficients
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: reactive_sdc_halo_filler,reactive_sdc_stage_agreement,reactive_sdc_source_evaluator, &
    react_sdc_state_3d,recover_reactive_temperatures_3d,compute_reactive_external_source_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d,rotating_momentum_source_3d, &
    convert_rotating_energy_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d, &
    fill_mach_inlet_boundary_ghosts_3d,fill_pressure_outlet_ghosts_3d,valid_reactive_boundaries_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  implicit none
  private
  public :: compute_reactive_mol_halo_3d,compute_reactive_mol_3d,advance_reactive_mol_3d
  public :: convert_regular_energy_3d
  public :: reactive_mol_blend_agreement
  abstract interface
    subroutine reactive_mol_blend_agreement(local_blend,global_blend,ok)
      import dp
      real(dp), intent(in) :: local_blend
      real(dp), intent(out) :: global_blend
      logical, intent(out) :: ok
    end subroutine reactive_mol_blend_agreement
  end interface
contains
  ! Regular Cartesian adapter to the shared MOL reconstruction/Riemann kernel.
  ! Input has exactly three same-resolution ghost cells on every side. Fluxes
  ! are dt-integrated per face area; rates are per volume per time. Domain lower
  ! denotes the first valid face, including for a distributed tile or rotation.
  subroutine compute_reactive_mol_halo_3d(species,ghost,ghost_t,n,spacing,dt,rate,fx,fy,fz,boundary,ok, &
      order,lax_friedrichs,domain_lower,frame,internal_energy_rate,internal_energy_integrated_flux, &
      internal_energy_ghost)
    real(dp), optional, intent(in) :: internal_energy_ghost(:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: ghost(:,:,:,:),ghost_t(:,:,:),spacing(3),dt
    integer, intent(in) :: n(3)
    real(dp), intent(inout) :: rate(:,:,:,:),fx(:,0:,:,:),fy(:,:,0:,:),fz(:,:,:,0:),boundary(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    logical, optional, intent(in) :: lax_friedrichs
    real(dp), optional, intent(in) :: domain_lower(3)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_integrated_flux
    type(eb_geometry_3d) :: geometry
    type(eb_integrated_flux_3d) :: flux,private_auxiliary_flux
    real(dp), allocatable :: private_auxiliary_rate(:,:,:)
    real(dp), allocatable :: private_rate(:,:,:,:)
    real(dp) :: lower(3),lo(3),hi(3),outside,private_boundary(size(boundary)),wall(size(boundary))
    integer :: nv,reconstruction
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));reconstruction=2
    if(present(order)) reconstruction=order
    if(any(n<1).or.any(n>huge(1)-6)) return
    if(.not.all(ieee_is_finite([spacing,dt])).or.any(spacing<=0).or.dt<=0) return
    if(any(shape(ghost)/=[nv,n+6]).or.any(shape(ghost_t)/=n+6)) return
    if(any(shape(rate)/=[nv,n]).or.size(boundary)/=nv) return
    if(any(shape(fx)/=[nv,n(1)+1,n(2),n(3)]).or.any(shape(fy)/=[nv,n(1),n(2)+1,n(3)]).or. &
      any(shape(fz)/=[nv,n(1),n(2),n(3)+1])) return
    if(present(internal_energy_ghost)) then
      if(any(shape(internal_energy_ghost)/=n+6)) return
      if(any(.not.ieee_is_finite(internal_energy_ghost))) return
    endif
    if(reconstruction/=1.and.reconstruction/=2) return
    lower=0
    if(present(domain_lower)) lower=domain_lower
    if(.not.all(ieee_is_finite(lower))) return
    lo=lower-3*spacing;hi=lower+(real(n,dp)+3)*spacing;outside=lo(1)-spacing(1)
    if(.not.all(ieee_is_finite([lo,hi,outside])).or.any(hi<=lo).or.outside>=lo(1)) return
    call build_axis_plane_eb_geometry_3d(n(1)+6,n(2)+6,n(3)+6,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3), &
      'x',outside,geometry,valid)
    if(.not.valid) return
    if(any(geometry%volume_fraction/=1)) return
    ! The global regular mesh owns these metrics. Reconstructing spacing from
    ! rank-local bounds introduces decomposition-dependent roundoff.
    geometry%dx=spacing(1);geometry%dy=spacing(2);geometry%dz=spacing(3)
    if(.not.geometry%is_valid()) return
    if(present(internal_energy_rate).neqv.present(internal_energy_integrated_flux)) return
    if(present(internal_energy_rate)) then
      if(any(shape(internal_energy_rate)/=n)) return
      allocate(private_auxiliary_rate(n(1),n(2),n(3)))
    endif
    allocate(private_rate(nv,n(1),n(2),n(3)))
    private_rate=0;private_boundary=0;wall=0
    if(present(internal_energy_rate)) then
      call compute_mol_eb_hydro_3d(species,geometry,ghost,ghost_t,dt,reconstruction,private_rate,flux, &
        private_boundary,wall,valid,lax_friedrichs,frame, &
        internal_energy_rate=private_auxiliary_rate,internal_energy_integrated_flux=private_auxiliary_flux, &
        internal_energy_density=internal_energy_ghost, &
        primitive_from_total_energy=.not.present(frame),preserve_species_ratios=.not.present(frame))
    else
      call compute_mol_eb_hydro_3d(species,geometry,ghost,ghost_t,dt,reconstruction,private_rate,flux, &
        private_boundary,wall,valid,lax_friedrichs,frame,internal_energy_density=internal_energy_ghost, &
        primitive_from_total_energy=.not.present(frame),preserve_species_ratios=.not.present(frame))
    endif
    if(.not.valid) return
    if(any(wall/=0).or..not.all(ieee_is_finite(private_rate)).or. &
      .not.all(ieee_is_finite(private_boundary))) return
    if(present(internal_energy_rate)) then
      internal_energy_rate=private_auxiliary_rate
      internal_energy_integrated_flux=private_auxiliary_flux
    endif
    rate=private_rate;fx=flux%x;fy=flux%y;fz=flux%z;boundary=private_boundary;ok=.true.
  end subroutine compute_reactive_mol_halo_3d

  ! A standalone regular domain fills physical ghosts itself. Distributed
  ! callers supply fill_halos and agree optional controls before entry; the
  ! agreement hook prevents a local failure from entering another halo call.
  subroutine compute_reactive_mol_3d(species,state,temperature,spacing,dt,time,rate,fx,fy,fz,boundary,ok, &
      order,lax_friedrichs,boundary_faces,inflow_state,inflow_temperature,inflow_mach,inflow_pressure, &
      outlet_pressure,wall_velocity,fill_halos,agree_stage,domain_lower,frame,internal_energy_rate, &
      internal_energy_integrated_flux, &
      internal_energy_density)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), allocatable :: held_ghost(:,:,:)
    real(dp), optional, intent(inout) :: internal_energy_rate(:,:,:)
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_integrated_flux
    type(eb_integrated_flux_3d) :: auxiliary_flux
    real(dp), allocatable :: auxiliary_rate(:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),dt,time
    real(dp), intent(inout) :: rate(:,:,:,:),fx(:,0:,:,:),fy(:,:,0:,:),fz(:,:,:,0:),boundary(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: order
    logical, optional, intent(in) :: lax_friedrichs
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),wall_velocity(:,:)
    procedure(reactive_sdc_halo_filler), optional :: fill_halos
    procedure(reactive_sdc_stage_agreement), optional :: agree_stage
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),private_rate(:,:,:,:),px(:,:,:,:),py(:,:,:,:),pz(:,:,:,:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    real(dp), allocatable :: ordinary(:,:,:,:)
    real(dp) :: lower(3),private_boundary(size(boundary))
    integer :: n(3),nv,reconstruction
    character(len=32) :: faces(6)
    logical :: valid
    ok=.false.;n=shape(temperature);nv=reactive_nvar(size(species));faces='periodic';reconstruction=2
    lower=0
    if(present(domain_lower)) lower=domain_lower
    if(present(order)) reconstruction=order
    if(present(boundary_faces)) faces=boundary_faces
    valid=all(n>=1).and.all(n<=huge(1)-6).and.all(shape(state)==[nv,n])
    valid=valid.and.all(shape(rate)==[nv,n]).and.size(boundary)==nv
    valid=valid.and.all(shape(fx)==[nv,n(1)+1,n(2),n(3)]).and.all(shape(fy)==[nv,n(1),n(2)+1,n(3)]).and. &
      all(shape(fz)==[nv,n(1),n(2),n(3)+1])
    valid=valid.and.all(ieee_is_finite([spacing,dt,time])).and.all(spacing>0).and.dt>0
    valid=valid.and.(reconstruction==1.or.reconstruction==2).and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.all(ieee_is_finite(lower))
    if(present(frame)) valid=valid.and.frame%is_valid()
    valid=valid.and.(present(internal_energy_rate).eqv.present(internal_energy_integrated_flux))
    if(present(internal_energy_rate)) valid=valid.and.all(shape(internal_energy_rate)==n)
    if(present(internal_energy_density)) then
      valid=valid.and.all(shape(internal_energy_density)==n)
      valid=valid.and.all(ieee_is_finite(internal_energy_density))
    endif
    if(.not.ready(valid)) return
    if(present(internal_energy_density)) allocate(held_ghost(n(1)+6,n(2)+6,n(3)+6))
    if(present(internal_energy_rate)) allocate(auxiliary_rate,mold=internal_energy_rate)
    allocate(ghost(nv,n(1)+6,n(2)+6,n(3)+6),gt(n(1)+6,n(2)+6,n(3)+6))
    if(present(fill_halos)) then
      call fill_halos(state,temperature,time,3,.false.,ghost,gt,valid, &
        internal_energy_density=internal_energy_density,internal_energy_ghost=held_ghost)
    else
      ! Local physical ghosts are built in ordinary relative energy. Prescribed
      ! inflow templates use that convention; published/callback states rotate.
      ordinary=state
      if(present(frame)) then
        call convert_regular_energy_3d(frame,lower,spacing,ordinary,0,.false.,valid)
        if(.not.ready(valid)) return
      endif
      if(present(inflow_mach)) then
        call fill_mach_inlet_boundary_ghosts_3d(species,ordinary,temperature,faces,inflow_mach,ghost,gt,valid, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          wall_velocity=wall_velocity,ghost_width=3, &
          internal_energy_density=internal_energy_density,internal_energy_ghost=held_ghost)
      else
        call fill_reactive_boundary_ghosts_3d(ordinary,temperature,faces,ghost,gt,valid, &
          inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=wall_velocity,ghost_width=3, &
          internal_energy_density=internal_energy_density,internal_energy_ghost=held_ghost)
      endif
      if(.not.ready(valid)) return
      call fill_pressure_outlet_ghosts_3d(species,n,faces,outlet_pressure,ghost,gt,valid,ghost_width=3, &
        internal_energy_ghost=held_ghost)
      if(.not.ready(valid)) return
      if(present(frame)) call convert_regular_energy_3d(frame,lower,spacing,ghost,3,.true.,valid)
    endif
    if(.not.ready(valid)) return
    allocate(private_rate,source=rate);allocate(px,source=fx);allocate(py,source=fy);allocate(pz,source=fz)
    if(present(internal_energy_rate)) then
      call compute_reactive_mol_halo_3d(species,ghost,gt,n,spacing,dt,private_rate,px,py,pz,private_boundary, &
        valid,reconstruction,lax_friedrichs,domain_lower=lower,frame=frame, &
        internal_energy_rate=auxiliary_rate,internal_energy_integrated_flux=auxiliary_flux, &
        internal_energy_ghost=held_ghost)
    else
      call compute_reactive_mol_halo_3d(species,ghost,gt,n,spacing,dt,private_rate,px,py,pz,private_boundary, &
        valid,reconstruction,lax_friedrichs,domain_lower=lower,frame=frame,internal_energy_ghost=held_ghost)
    endif
    if(.not.ready(valid)) return
    if(present(internal_energy_rate)) then
      internal_energy_rate=auxiliary_rate;internal_energy_integrated_flux=auxiliary_flux
    endif
    rate=private_rate;fx=px;fy=py;fz=pz;boundary=private_boundary;ok=.true.
  contains
    logical function ready(local_ok)
      logical, intent(in) :: local_ok
      ready=local_ok
      if(present(agree_stage)) call agree_stage(local_ok,ready)
    end function ready
  end subroutine compute_reactive_mol_3d
  ! Pinned MOL predictor/trapezoidal reaction iteration. All caller fields and
  ! ledgers are published only after the complete interval accepts. Callbacks
  ! own private physical side effects until this routine returns success.
  subroutine advance_reactive_mol_3d(species,state,temperature,spacing,dt,time,iterations,chemistry,rtol,atol,ok, &
      reactions,reaction_source,heat_release,transport,transport_flags,source_count,evaluate_sources,external_force, &
      source_change,boundary_change,chemical_change,integrated_flux,registered_flux,minimum_transport_theta, &
      order,lax_friedrichs,boundary_faces,inflow_state,inflow_temperature,inflow_mach,inflow_pressure, &
      outlet_pressure,wall_temperature,wall_velocity,constant_coefficients,fill_halos,agree_stage,evaluate_reaction, &
      physical_faces,domain_lower,frame,internal_energy_density,energy_correction,energy_policy, &
      internal_energy_integrated_flux,internal_energy_registered_flux,agree_flux_blend,minimum_flux_blend)
    procedure(reactive_mol_blend_agreement), optional :: agree_flux_blend
    real(dp), optional, intent(inout) :: minimum_flux_blend
    real(dp), optional, intent(inout) :: internal_energy_density(:,:,:),energy_correction
    type(eb_integrated_flux_3d), optional, intent(inout) :: internal_energy_integrated_flux,internal_energy_registered_flux
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    type(dual_energy_policy) :: policy
    real(dp), allocatable :: old_auxiliary(:,:,:),auxiliary(:,:,:),a0(:,:,:),a1(:,:,:), &
      predictor(:,:,:,:),predictor_auxiliary(:,:,:),predictor_reset(:,:,:),reset_delta(:,:,:)
    type(eb_integrated_flux_3d) :: af0,af1,auxiliary_physical_flux,auxiliary_registration
    real(dp) :: correction,initial_correction,predictor_correction,final_correction,cell_change,cell_potential
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt,time,rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry
    logical, intent(out) :: ok
    type(elementary_reaction), optional, intent(in) :: reactions(:)
    real(dp), optional, intent(inout) :: reaction_source(:,:,:,:),heat_release(:,:,:)
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: transport_flags(4),physical_faces(6)
    integer, optional, intent(in) :: source_count,order
    procedure(reactive_sdc_source_evaluator), optional :: evaluate_sources
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(reactive_sdc_halo_filler), optional :: fill_halos
    procedure(reactive_sdc_stage_agreement), optional :: agree_stage
    real(dp), optional, intent(in) :: external_force(3),inflow_state(:,:),inflow_temperature(:),inflow_mach(6)
    real(dp), optional, intent(in) :: inflow_pressure(6),outlet_pressure(6),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    logical, optional, intent(in) :: lax_friedrichs
    real(dp), optional, intent(inout) :: source_change(:),boundary_change(:),chemical_change(:),minimum_transport_theta
    type(eb_integrated_flux_3d), optional, intent(inout) :: integrated_flux,registered_flux
    real(dp), allocatable :: old(:,:,:,:),ot(:,:,:),gas(:,:,:,:),temp(:,:,:),recovered(:,:,:),ir(:,:,:,:),heat(:,:,:)
    real(dp), allocatable :: h0(:,:,:,:),h1(:,:,:,:),e0(:,:,:,:),e1(:,:,:,:),rates(:,:,:,:),families(:,:,:,:,:)
    type(eb_integrated_flux_3d) :: f0,f1,physical_flux,registration
    real(dp), allocatable :: b0(:),b1(:),boundary(:),source(:),chemical(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    real(dp), allocatable :: potential(:,:,:)
    real(dp) :: lower(3),position(3),force(3),theta,stage_theta,volume,accepted_blend
    integer :: n(3),nv,nfamilies,reconstruction,it,c,i,j,k
    logical :: valid,flags(4),physical(6),all_blend,all_no_blend
    character(len=32) :: faces(6),transport_faces(6)
    ok=.false.;n=shape(temperature);nv=reactive_nvar(size(species));nfamilies=0;reconstruction=2
    if(present(source_count)) nfamilies=source_count
    if(present(order)) reconstruction=order
    faces='periodic';force=0;flags=.false.;physical=.true.;lower=0
    if(present(domain_lower)) lower=domain_lower
    if(present(boundary_faces)) faces=boundary_faces
    if(present(external_force)) force=external_force
    if(present(transport_flags)) flags=transport_flags
    if(present(physical_faces)) physical=physical_faces
    valid=all(n>=1).and.all(n<=huge(1)-6).and.all(shape(state)==[nv,n]).and.iterations>=1.and.nfamilies>=0
    valid=valid.and.all(ieee_is_finite([spacing,dt,time,rtol,atol,force])).and.all(spacing>0).and.min(dt,rtol,atol)>0
    valid=valid.and.valid_reactive_boundaries_3d(faces).and.(reconstruction==1.or.reconstruction==2)
    valid=valid.and.(present(evaluate_sources).eqv.nfamilies>0)
    valid=valid.and..not.(any(flags).and..not.present(transport)).and..not.(flags(4).and..not.flags(3))
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    valid=valid.and.(all(physical).or.present(fill_halos))
    if(present(reaction_source)) valid=valid.and.all(shape(reaction_source)==[nv-4,n])
    if(present(heat_release)) valid=valid.and.all(shape(heat_release)==n)
    if(present(source_change)) valid=valid.and.size(source_change)==nv
    if(present(boundary_change)) valid=valid.and.size(boundary_change)==nv
    if(present(chemical_change)) valid=valid.and.size(chemical_change)==nv
    valid=valid.and.all(ieee_is_finite(lower))
    if(present(frame)) valid=valid.and.frame%is_valid()
    valid=valid.and.(present(internal_energy_density).eqv.present(energy_correction))
    valid=valid.and.(.not.present(energy_policy).or.present(internal_energy_density))
    valid=valid.and.(.not.(present(internal_energy_integrated_flux).or.present(internal_energy_registered_flux)).or. &
      present(internal_energy_density))
    if(present(internal_energy_density)) then
      valid=valid.and.all(shape(internal_energy_density)==n)
      valid=valid.and.all(ieee_is_finite(internal_energy_density)).and.ieee_is_finite(energy_correction)
    endif
    if(.not.ready(valid)) return
    if(present(agree_stage)) then
      call agree_stage(present(agree_flux_blend),all_blend)
      call agree_stage(.not.present(agree_flux_blend),all_no_blend)
      if(.not.all_blend.and..not.all_no_blend) return
    endif
    if(present(energy_policy)) policy=energy_policy
    if(present(internal_energy_density)) then
      allocate(old_auxiliary,source=internal_energy_density);allocate(auxiliary,source=internal_energy_density)
      allocate(a0,mold=auxiliary);allocate(a1,mold=auxiliary);allocate(reset_delta,mold=auxiliary)
    endif
    volume=product(spacing);correction=0;initial_correction=0;predictor_correction=0;final_correction=0
    if(chemistry) then
      if(.not.ready(present(reactions))) return
      if(.not.ready(size(reactions)>0)) return
    endif
    if(present(transport)) then
      if(.not.ready(size(transport)==size(species))) return
    endif
    if(present(constant_coefficients)) then
      if(.not.ready(present(transport))) return
      if(.not.ready(valid_constant_transport_coefficients(constant_coefficients))) return
    endif
    allocate(old,source=state);allocate(ot,source=temperature);allocate(gas,source=state);allocate(temp,source=temperature)
    allocate(recovered,source=temperature);allocate(heat,source=temperature);heat=0
    allocate(ir(nv-4,n(1),n(2),n(3)),source=0._dp)
    if(chemistry.and.present(reaction_source)) ir=reaction_source
    if(.not.ready(all(ieee_is_finite(ir)))) return
    if(present(frame)) then
      allocate(potential(n(1),n(2),n(3)))
      valid=.true.
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        position=lower+(real([i,j,k],dp)-.5_dp)*spacing
        call rotating_potential_3d(frame,position,potential(i,j,k),valid)
        if(.not.valid) exit
      enddo
        if(.not.valid) exit
      enddo
        if(.not.valid) exit
      enddo
      if(.not.ready(valid)) return
    endif
    call recover(old,ot,valid,old_auxiliary)
    if(.not.ready(valid)) return
    if(present(internal_energy_density)) initial_correction=volume*sum(reset_delta)
    allocate(h0,source=state);allocate(h1,source=state);allocate(e0,source=state);allocate(e1,source=state)
    allocate(rates,source=state);allocate(families(nv,n(1),n(2),n(3),nfamilies))
    allocate(b0(nv),b1(nv),boundary(nv),source(nv),chemical(nv))
    call allocate_flux(f0);call allocate_flux(f1)
    theta=1;accepted_blend=1;transport_faces=faces
    where(.not.physical) transport_faces='interior'
    call spatial(old,ot,time,h0,f0,b0,valid,a0,af0,old_auxiliary)
    if(.not.ready(valid)) return
    call external(old,ot,time,e0,valid)
    if(.not.ready(valid)) return
    gas=old+dt*(h0+e0)
    if(chemistry) then
      gas(6:,:,:,:)=gas(6:,:,:,:)+dt*ir(:nv-5,:,:,:)
      gas(5,:,:,:)=gas(5,:,:,:)+dt*ir(nv-4,:,:,:)
    endif
    if(present(internal_energy_density)) auxiliary=old_auxiliary+dt*a0
    call recover(gas,temp,valid,auxiliary)
    if(.not.ready(valid)) return
    if(present(internal_energy_density)) then
      predictor=gas;predictor_auxiliary=auxiliary;predictor_reset=reset_delta
      predictor_correction=volume*sum(reset_delta)
    endif
    do it=1,merge(iterations,1,chemistry)
      call spatial(gas,temp,time+dt,h1,f1,b1,valid,a1,af1,auxiliary)
      if(.not.ready(valid)) return
      e1=0
      if(it==1) then
        call external(gas,temp,time+dt,e1,valid)
        if(.not.ready(valid)) return
      endif
      rates=.5_dp*(h0+h1)+.5_dp*(e0+e1)
      if(present(internal_energy_density).and.it==1) rates(5,:,:,:)=rates(5,:,:,:)+.5_dp*predictor_reset/dt
      if(.not.ready(all(ieee_is_finite(rates)))) return
      if(chemistry) then
        call react_sdc_state_3d(species,reactions,old,ot,rates,dt,rtol,atol,gas,temp,ir,valid, &
          evaluate_reaction=evaluate_reaction,heat_release=heat,rotating_potential=potential, &
          internal_energy_density=auxiliary,old_internal_energy_density=old_auxiliary)
        valid=ready(valid)
        if(.not.valid) then
          ! A dry old cell can lose a tracer made by a completed new-time
          ! particle kick. Freeze the evaluated sources and repair the whole
          ! conservative spatial family; do not replay the particle callback.
          call repair_flux_window(valid)
          if(.not.ready(valid)) return
          rates=.5_dp*(h0+h1)+.5_dp*(e0+e1)
          if(present(internal_energy_density).and.it==1) rates(5,:,:,:)=rates(5,:,:,:)+.5_dp*predictor_reset/dt
          call react_sdc_state_3d(species,reactions,old,ot,rates,dt,rtol,atol,gas,temp,ir,valid, &
            evaluate_reaction=evaluate_reaction,heat_release=heat,rotating_potential=potential, &
            internal_energy_density=auxiliary,old_internal_energy_density=old_auxiliary)
          if(.not.ready(valid)) return
        endif
        if(present(internal_energy_density)) then
          ! The reactor adapter returns the forced internal-energy endpoint
          ! before adding kinetic/rotating energy (pinned React.cpp UEINT).
          call recover(gas,temp,valid,auxiliary)
          if(.not.ready(valid)) return
        endif
      else
        if(present(internal_energy_density)) then
          gas=.5_dp*(old+predictor)+.5_dp*dt*(h1+e1)
          auxiliary=.5_dp*(old_auxiliary+predictor_auxiliary)+.5_dp*dt*a1
        else
          gas=old+dt*rates
        endif
        call recover(gas,temp,valid,auxiliary)
        valid=ready(valid)
        if(.not.valid) then
          call repair_flux_window(valid)
          if(.not.ready(valid)) return
          rates=.5_dp*(h0+h1)+.5_dp*(e0+e1)
          if(present(internal_energy_density).and.it==1) rates(5,:,:,:)=rates(5,:,:,:)+.5_dp*predictor_reset/dt
          if(present(internal_energy_density)) then
            gas=.5_dp*(old+predictor)+.5_dp*dt*(h1+e1)
            auxiliary=.5_dp*(old_auxiliary+predictor_auxiliary)+.5_dp*dt*a1
          else
            gas=old+dt*rates
          endif
          call recover(gas,temp,valid,auxiliary)
          if(.not.ready(valid)) return
        endif
      endif
    enddo
    if(present(internal_energy_density)) then
      final_correction=volume*sum(reset_delta)
      correction=initial_correction+final_correction
      if(.not.chemistry.or.iterations==1) correction=correction+.5_dp*predictor_correction
    endif
    if(.not.ready(all(ieee_is_finite(gas)).and.all(ieee_is_finite(temp)).and.all(gas(6:,:,:,:)>=0))) return
    volume=product(spacing);source=0;chemical=0;boundary=.5_dp*(b0+b1)
    do c=1,nv
      source(c)=.5_dp*dt*volume*(sum(e0(c,:,:,:))+sum(e1(c,:,:,:)))
      chemical(c)=volume*sum(gas(c,:,:,:)-(old(c,:,:,:)+dt*rates(c,:,:,:)))
    enddo
    if(present(internal_energy_density)) then
      chemical(5)=chemical(5)-final_correction
      auxiliary_physical_flux=af0
      auxiliary_physical_flux%x=.5_dp*(af0%x+af1%x)
      auxiliary_physical_flux%y=.5_dp*(af0%y+af1%y)
      auxiliary_physical_flux%z=.5_dp*(af0%z+af1%z)
      auxiliary_registration=auxiliary_physical_flux
      if(.not.chemistry.and.iterations>1) then
        auxiliary_registration%x=.5_dp*af0%x
        auxiliary_registration%y=.5_dp*af0%y
        auxiliary_registration%z=.5_dp*af0%z
      endif
    endif
    physical_flux=f0
    physical_flux%x=.5_dp*(f0%x+f1%x);physical_flux%y=.5_dp*(f0%y+f1%y);physical_flux%z=.5_dp*(f0%z+f1%z)
    registration=physical_flux
    if(.not.chemistry.and.iterations>1) then
      registration%x=.5_dp*f0%x;registration%y=.5_dp*f0%y;registration%z=.5_dp*f0%z
    endif
    if(.not.ready(all(ieee_is_finite([boundary,source,chemical])).and.all(ieee_is_finite(heat)))) return
    if(present(internal_energy_density)) then
      if(.not.ready(ieee_is_finite(energy_correction+correction))) return
      internal_energy_density=auxiliary;energy_correction=energy_correction+correction
      if(present(internal_energy_integrated_flux)) internal_energy_integrated_flux=auxiliary_physical_flux
      if(present(internal_energy_registered_flux)) internal_energy_registered_flux=auxiliary_registration
    endif
    state=gas;temperature=temp
    if(present(reaction_source)) reaction_source=ir
    if(present(heat_release)) heat_release=heat
    if(present(source_change)) source_change=source
    if(present(boundary_change)) boundary_change=boundary
    if(present(chemical_change)) chemical_change=chemical
    if(present(integrated_flux)) integrated_flux=physical_flux
    if(present(registered_flux)) registered_flux=registration
    if(present(minimum_transport_theta)) minimum_transport_theta=theta
    if(present(minimum_flux_blend)) minimum_flux_blend=accepted_blend
    ok=.true.
  contains
    subroutine repair_flux_window(good)
      logical, intent(out) :: good
      real(dp) :: baseline(nv),candidate(nv),local_blend,blend,ratio,allowance
      integer :: x,y,z,component
      logical :: local_ok
      local_blend=1;good=.true.
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        baseline=old(:,x,y,z)+dt*(h0(:,x,y,z)+.5_dp*(e0(:,x,y,z)+e1(:,x,y,z)))
        candidate=old(:,x,y,z)+dt*rates(:,x,y,z)
        if(.not.all(ieee_is_finite([baseline,candidate]))) then
          good=.false.;cycle
        endif
        if(baseline(1)<=0.or.any(baseline(6:)<0)) then
          good=.false.;cycle
        endif
        do component=1,nv
          if(component>=2.and.component<=5) cycle
          if(candidate(component)>=0.and.(component/=1.or.candidate(component)>0)) cycle
          ! Evaluate base/(base-high) without overflowing the subtraction.
          allowance=0
          if(baseline(component)>0) then
            if(baseline(component)>-candidate(component)) then
              allowance=1/(1+(-candidate(component))/baseline(component))
            else
              ratio=baseline(component)/(-candidate(component))
              allowance=ratio/(1+ratio)
            endif
          endif
          local_blend=min(local_blend,allowance*(1-16*epsilon(1._dp)))
        enddo
      enddo
      enddo
      enddo
      good=ready(good)
      if(.not.good) return
      blend=local_blend
      if(present(agree_flux_blend)) then
        call agree_flux_blend(local_blend,blend,local_ok)
        good=local_ok
        if(good) good=ieee_is_finite(blend)
        if(good) good=blend>=0.and.blend<=local_blend
      else if(present(agree_stage)) then
        ! A distributed caller needs one coefficient for shared faces.
        good=local_blend>=1
      endif
      good=ready(good)
      if(.not.good) return
      good=blend<1
      if(.not.good) return
      if(blend==0) then
        h1=h0;f1=f0;b1=b0
        if(present(internal_energy_density)) then
          a1=a0;af1=af0
        endif
      else
        h1=h0+blend*(h1-h0);b1=b0+blend*(b1-b0)
        f1%x=f0%x+blend*(f1%x-f0%x)
        f1%y=f0%y+blend*(f1%y-f0%y)
        f1%z=f0%z+blend*(f1%z-f0%z)
        if(present(internal_energy_density)) then
          a1=a0+blend*(a1-a0)
          af1%x=af0%x+blend*(af1%x-af0%x)
          af1%y=af0%y+blend*(af1%y-af0%y)
          af1%z=af0%z+blend*(af1%z-af0%z)
        endif
      endif
      accepted_blend=min(accepted_blend,blend)
      good=all(ieee_is_finite(h1)).and.all(ieee_is_finite(b1))
    end subroutine repair_flux_window
    logical function ready(local_ok)
      logical, intent(in) :: local_ok
      ready=local_ok
      if(present(agree_stage)) call agree_stage(local_ok,ready)
    end function ready
    subroutine allocate_flux(f)
      type(eb_integrated_flux_3d), intent(out) :: f
      allocate(f%x(nv,0:n(1),n(2),n(3)),f%y(nv,n(1),0:n(2),n(3)),f%z(nv,n(1),n(2),0:n(3)))
      f%x=0;f%y=0;f%z=0
    end subroutine allocate_flux
    subroutine recover(u,t,good,aux)
      real(dp), intent(inout) :: u(:,:,:,:)
      real(dp), optional, intent(inout) :: aux(:,:,:)
      real(dp), intent(inout) :: t(:,:,:)
      logical, intent(out) :: good
      good=all(ieee_is_finite(u))
      if(.not.good) return
      good=all(u(1,:,:,:)>0).and.all(u(6:,:,:,:)>=0)
      if(.not.good) return
      if(present(aux)) then
        reset_delta=0
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          cell_change=0;cell_potential=0
          if(allocated(potential)) cell_potential=potential(i,j,k)
          call reset_nasa7_dual_energy_cell(species,u(:,i,j,k),aux(i,j,k),t(i,j,k),policy, &
            cell_change,good,cell_potential)
          if(.not.good) return
          reset_delta(i,j,k)=cell_change
        enddo
        enddo
        enddo
      endif
      call recover_reactive_temperatures_3d(species,u,t,n(1),n(2),n(3),recovered,good,potential,aux)
      if(good) t=recovered
    end subroutine recover
    subroutine external(u,t,stage_time,output,good)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time
      real(dp), intent(out) :: output(:,:,:,:)
      logical, intent(out) :: good
      real(dp) :: rotation_source(3),cell_position(3)
      integer :: x,y,z
      call compute_reactive_external_source_3d(u,force,output,good)
      good=ready(good)
      if(.not.good) return
      if(present(frame)) then
        do z=1,n(3)
        do y=1,n(2)
        do x=1,n(1)
          cell_position=lower+(real([x,y,z],dp)-.5_dp)*spacing
          call rotating_momentum_source_3d(frame,cell_position,u(1,x,y,z),u(2:4,x,y,z),rotation_source,good)
          if(.not.good) exit
          output(2:4,x,y,z)=output(2:4,x,y,z)+rotation_source
        enddo
          if(.not.good) exit
        enddo
          if(.not.good) exit
        enddo
        if(.not.ready(good)) return
      endif
      if(present(evaluate_sources)) then
        call evaluate_sources(u,t,stage_time,families,good)
        good=ready(good)
      if(.not.good) return
        output=output+sum(families,dim=5)
      endif
      good=all(ieee_is_finite(output))
    end subroutine external
    subroutine spatial(u,t,stage_time,output,f,b,good,aux_rate,aux_flux,stage_held)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stage_time
      real(dp), optional, intent(in) :: stage_held(:,:,:)
      real(dp), allocatable :: ghost_held(:,:,:)
      real(dp), optional, intent(inout) :: aux_rate(:,:,:)
      type(eb_integrated_flux_3d), optional, intent(inout) :: aux_flux
      real(dp), intent(inout) :: output(:,:,:,:)
      type(eb_integrated_flux_3d), intent(inout) :: f
      real(dp), intent(inout) :: b(:)
      logical, intent(out) :: good
      real(dp), allocatable :: tx(:,:,:,:),ty(:,:,:,:),tz(:,:,:,:),ghost(:,:,:,:),gt(:,:,:),ordinary(:,:,:,:)
      if(present(aux_rate)) then
        call compute_reactive_mol_3d(species,u,t,spacing,dt,stage_time,output,f%x,f%y,f%z,b,good, &
          reconstruction,lax_friedrichs,faces,inflow_state,inflow_temperature,inflow_mach,inflow_pressure, &
          outlet_pressure,wall_velocity,fill_halos,agree_stage,domain_lower=lower,frame=frame, &
          internal_energy_rate=aux_rate,internal_energy_integrated_flux=aux_flux,internal_energy_density=stage_held)
      else
        call compute_reactive_mol_3d(species,u,t,spacing,dt,stage_time,output,f%x,f%y,f%z,b,good, &
          reconstruction,lax_friedrichs,faces,inflow_state,inflow_temperature,inflow_mach,inflow_pressure, &
          outlet_pressure,wall_velocity,fill_halos,agree_stage,domain_lower=lower,frame=frame, &
          internal_energy_density=stage_held)
      endif
      good=ready(good)
      if(.not.good) return
      if(.not.present(transport)) return
      allocate(tx(nv,0:n(1),n(2),n(3)),ty(nv,n(1),0:n(2),n(3)),tz(nv,n(1),n(2),0:n(3)))
      ordinary=u
      if(present(frame)) then
        call convert_regular_energy_3d(frame,lower,spacing,ordinary,0,.false.,good)
        if(.not.ready(good)) return
      endif
      if(present(fill_halos)) then
        allocate(ghost(nv,n(1)+4,n(2)+4,n(3)+4),gt(n(1)+4,n(2)+4,n(3)+4))
        if(present(stage_held)) allocate(ghost_held(n(1)+4,n(2)+4,n(3)+4))
        call fill_halos(u,t,stage_time,2,.true.,ghost,gt,good, &
          internal_energy_density=stage_held,internal_energy_ghost=ghost_held)
        good=ready(good)
      if(.not.good) return
        if(present(frame)) then
          call convert_regular_energy_3d(frame,lower,spacing,ghost,2,.false.,good)
          if(.not.ready(good)) return
        endif
        call physical_transport_ghosted_fluxes_3d(species,transport,ordinary,t,ghost,gt,n,spacing,dt,transport_faces, &
          flags(1),flags(2),flags(3),flags(4),tx,ty,tz,stage_theta,good,wall_temperature,wall_velocity,constant_coefficients, &
          internal_energy_density=stage_held,internal_energy_ghost=ghost_held)
      else
        call physical_transport_fluxes_3d(species,transport,ordinary,t,n,spacing,dt,faces, &
          flags(1),flags(2),flags(3),flags(4),tx,ty,tz,stage_theta,good, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure, &
          constant_coefficients,internal_energy_density=stage_held)
      endif
      good=ready(good)
      if(.not.good) return
      theta=min(theta,stage_theta)
      f%x=f%x+dt*tx;f%y=f%y+dt*ty;f%z=f%z+dt*tz
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        output(:,i,j,k)=output(:,i,j,k)-(tx(:,i,j,k)-tx(:,i-1,j,k))/spacing(1)- &
          (ty(:,i,j,k)-ty(:,i,j-1,k))/spacing(2)-(tz(:,i,j,k)-tz(:,i,j,k-1))/spacing(3)
      enddo
      enddo
      enddo
      do c=1,nv
        b(c)=((sum(f%x(c,0,:,:))-sum(f%x(c,n(1),:,:)))*spacing(2)*spacing(3)+ &
          (sum(f%y(c,:,0,:))-sum(f%y(c,:,n(2),:)))*spacing(1)*spacing(3)+ &
          (sum(f%z(c,:,:,0))-sum(f%z(c,:,:,n(3))))*spacing(1)*spacing(2))
      enddo
      good=all(ieee_is_finite(output)).and.all(ieee_is_finite(b))
    end subroutine spatial
  end subroutine advance_reactive_mol_3d

  ! Convert private regular/halo fields at their actual cell positions. Keeping
  ! ordinary EOS ghosts separate avoids wrapping a spatial potential with a
  ! periodic state. Distributed halo callbacks own the rotating convention.
  subroutine convert_regular_energy_3d(frame,lower,spacing,u,width,to_rotating,ok)
    type(rotating_frame_3d), intent(in) :: frame
    real(dp), intent(in) :: lower(3),spacing(3)
    real(dp), intent(inout) :: u(:,:,:,:)
    integer, intent(in) :: width
    logical, intent(in) :: to_rotating
    logical, intent(out) :: ok
    real(dp), allocatable :: energy(:,:,:)
    real(dp) :: position(3)
    integer :: i,j,k
    ok=.false.
    if(size(u,1)<5.or.width<0) return
    if(.not.frame%is_valid().or..not.all(ieee_is_finite([lower,spacing])).or.any(spacing<=0)) return
    allocate(energy(size(u,2),size(u,3),size(u,4)))
    do k=1,size(u,4)
    do j=1,size(u,3)
    do i=1,size(u,2)
      position=lower+(real([i,j,k]-width,dp)-.5_dp)*spacing
      call convert_rotating_energy_3d(frame,position,u(1,i,j,k),u(5,i,j,k),to_rotating,energy(i,j,k),ok)
      if(.not.ok) return
    enddo
    enddo
    enddo
    u(5,:,:,:)=energy;ok=.true.
  end subroutine convert_regular_energy_3d
end module reactive_mol_3d_mod
