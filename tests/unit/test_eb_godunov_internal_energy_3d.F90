program test_eb_godunov_internal_energy_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved,reactive_pelec_primitive_flux_x
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_hydro_3d_mod, only: eb_plm_face_prediction_3d,compute_eb_godunov_pressure_work_3d, &
    compute_sdc_eb_hydro_3d,compute_reactive_eb_plm_ctu_fluxes_3d,compute_sdc_eb_boundary_hydro_3d, &
    correct_reactive_eb_spatial_hydro_flux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer, parameter :: n=10,m=n-6
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: g,physical
  type(eb_integrated_flux_3d) :: flux,base,hflux,hbase,kept_flux,kept_base,kept_hflux,kept_hbase
  type(eb_plm_face_prediction_3d) :: prediction(3)
  real(dp), allocatable :: ghost(:,:,:,:),temp(:,:,:),vf(:,:,:),held(:,:,:),flat(:,:,:),source(:,:,:,:)
  real(dp), allocatable :: rate(:,:,:,:),kept_rate(:,:,:,:),hrate(:,:,:),pwork(:,:,:),kept_hrate(:,:,:),kept_pwork(:,:,:)
  real(dp), allocatable :: q(:),ql(:),qr(:),riemann(:),sound(:,:,:),boundary(:),wall(:),kept_boundary(:),kept_wall(:)
  real(dp) :: nan_value,cs,swap,expected_h,inventory,ledger,coefficient,original,cellvol,maxres
  integer :: nv,i,j,k,c,mode,policy,axis,idx(3),cell(3),inside(3),dump_unit,ios
  logical :: ok
  character(len=1024) :: dump_path
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species));nan_value=ieee_value(0._dp,ieee_quiet_nan)
  allocate(ghost(nv,-2:n+3,-2:n+3,-2:n+3),temp(-2:n+3,-2:n+3,-2:n+3), &
    vf(-2:n+3,-2:n+3,-2:n+3),held(-2:n+3,-2:n+3,-2:n+3),sound(-2:n+3,-2:n+3,-2:n+3), &
    flat(0:n+1,0:n+1,0:n+1),source(nv,n,n,n),rate(nv,m,m,m),hrate(m,m,m),pwork(m,m,m), &
    q(nv),ql(nv),qr(nv),riemann(nv),boundary(nv),wall(nv))
  dump_unit=0
  if(command_argument_count()==1) then
    call get_command_argument(1,dump_path)
    open(newunit=dump_unit,file=trim(dump_path),status='replace',action='write',iostat=ios)
    call require(ios==0,'pressure-work reference fixture file')
    write(dump_unit,*) 3
  endif
  call check_pressure_work(0)
  call check_pressure_work(1)
  call check_pressure_work(2)
  if(dump_unit/=0) close(dump_unit)
  maxres=0
  do mode=0,3
    call build_axis_plane_eb_geometry_3d(n,n,n,0._dp,1._dp,0._dp,2._dp,0._dp,3._dp, &
      'x',merge(-.1_dp,.46_dp,mode==0),g,ok)
    if(mode==2) call build_axis_plane_eb_geometry_3d(n,n,n,0._dp,1._dp,0._dp,2._dp,0._dp,3._dp, &
      'y',.92_dp,g,ok)
    if(mode==3) call build_axis_plane_eb_geometry_3d(n,n,n,0._dp,1._dp,0._dp,2._dp,0._dp,3._dp, &
      'z',1.38_dp,g,ok)
    call require(ok,'regular/cut spatial geometry')
    call extract_eb_geometry_3d(g,[4,4,4],[n-3,n-3,n-3],physical,ok)
    call require(ok,'physical spatial geometry')
    cellvol=physical%dx*physical%dy*physical%dz
    ghost=nan_value;temp=nan_value;held=nan_value;sound=0;flat=1;source=0
    do k=-2,n+3
    do j=-2,n+3
    do i=-2,n+3
      inside=max(1,min(n,[i,j,k]));vf(i,j,k)=g%volume_fraction(inside(1),inside(2),inside(3))
      if(vf(i,j,k)<=0) cycle
      q=0;q(1)=1+.001_dp*sum(inside);q(2:4)=[4._dp,2._dp,-1._dp]-.1_dp*inside
      q(5)=150000+100*inside(2)
      do c=1,size(species)
        if(species(c)%name=='N2') q(c+5)=1
      enddo
      call reactive_primitive_to_conserved(species,q,ghost(:,i,j,k),temp(i,j,k),cs,ok)
      call require(ok,'actual gas halo state')
      sound(i,j,k)=cs
      ! Distinct held energy, deliberately inconsistent with UEDEN-kinetic.
      held(i,j,k)=ghost(5,i,j,k)-.5_dp*sum(ghost(2:4,i,j,k)**2)/ghost(1,i,j,k)+123*inside(1)
    enddo
    enddo
    enddo
    do policy=0,1
      coefficient=.1_dp*policy
      source=0;source(2,:,:,:)=.3_dp*policy;source(5,:,:,:)=10000._dp*policy
      call compute_sdc_eb_hydro_3d(species,g,ghost,temp,vf,1e-7_dp,flat,source,rate,flux,boundary,wall,ok, &
        difmag=coefficient,baseline_flux=base)
      call require(ok,'legacy gas spatial family')
      kept_rate=rate;kept_flux=flux;kept_base=base;kept_boundary=boundary;kept_wall=wall
      call evaluate_held()
      call require(ok,'paired Godunov spatial family')
      call check_gas_identity()
      call check_weights()
      inventory=1e-7_dp*cellvol*sum(physical%volume_fraction*(hrate+pwork))
      ledger=physical%dy*physical%dz*sum(physical%x_face_fraction(0,:,:)*hflux%x(1,0,:,:)- &
        physical%x_face_fraction(m,:,:)*hflux%x(1,m,:,:))+ &
        physical%dx*physical%dz*sum(physical%y_face_fraction(:,0,:)*hflux%y(1,:,0,:)- &
        physical%y_face_fraction(:,m,:)*hflux%y(1,:,m,:))+ &
        physical%dx*physical%dy*sum(physical%z_face_fraction(:,:,0)*hflux%z(1,:,:,0)- &
        physical%z_face_fraction(:,:,m)*hflux%z(1,:,:,m))
      maxres=max(maxres,abs(inventory-ledger))
      call require(abs(inventory-ledger)<1e-15_dp,'independent held advective inventory plus pressure work')
      call require(maxval(abs(hflux%x))+maxval(abs(hflux%y))+maxval(abs(hflux%z))>0,'nonzero actual held advection')
      call require(maxval(abs(pwork))>1,'nonzero Godunov pressure work')
      call require(all(pack(hrate,physical%volume_fraction<=0)==0),'covered held rate masked')
      if(mode==0.and.policy==1) then
        kept_hflux=hflux;kept_pwork=pwork
        held=held+17._dp*spread(spread([(real(i,dp),i=-2,n+3)],2,n+6),3,n+6)
        call evaluate_held()
        call require(ok,'changed actual held state')
        call check_gas_identity()
        call require(all(pwork==kept_pwork),'held viscosity leaves pressure work unchanged')
        ! Uniform compression: divergence=-.1*(1/dx+1/dy+1/dz).
        expected_h=1e-7_dp*g%dx*coefficient*(-.1_dp*(1/g%dx+1/g%dy+1/g%dz))*17
        call require(maxval(abs((hflux%x-kept_hflux%x)-expected_h))<1e-15_dp, &
          'viscosity consumes independent held jump and exact gas divergence')
        call require(all(hflux%y==kept_hflux%y).and.all(hflux%z==kept_hflux%z), &
          'independent x held gradient leaves other fluxes unchanged')
      endif
      kept_hflux=hflux;kept_hbase=hbase;kept_hrate=hrate;kept_pwork=pwork
      original=held(n-3,n-3,n-3);held(n-3,n-3,n-3)=nan_value
      call evaluate_held()
      call require(.not.ok,'late active held NaN refusal')
      held(n-3,n-3,n-3)=original
      call check_all_unchanged()
      call compute_sdc_eb_hydro_3d(species,g,ghost,temp,vf,1e-7_dp,flat,source,rate,flux,boundary,wall,ok, &
        internal_energy_state=held,internal_energy_rate=hrate)
      call require(.not.ok,'incomplete held family refusal')
      call check_all_unchanged()
    enddo
    call compute_reactive_eb_plm_ctu_fluxes_3d(species,g,ghost,temp,vf,1e-7_dp,flat,prediction,ok, &
      conserved_source=source,capture_internal_energy_flux=.true.)
    call require(ok,'captured actual corrected Riemann flux')
    do axis=1,3
      idx=[3,3,3];idx(axis)=2
      if(mode>0) idx(mode)=3
      ql=prediction(axis)%left(:,idx(1),idx(2),idx(3));qr=prediction(axis)%right(:,idx(1),idx(2),idx(3))
      if(ql(1)<=0.or.qr(1)<=0) cycle
      swap=ql(2);ql(2)=ql(axis+1);ql(axis+1)=swap
      swap=qr(2);qr(2)=qr(axis+1);qr(axis+1)=swap
      cell=idx+2;cell(axis)=cell(axis)+1
      cs=.5_dp*(sound(idx(1)+2,idx(2)+2,idx(3)+2)+sound(cell(1),cell(2),cell(3)))
      call reactive_pelec_primitive_flux_x(species,ql,qr,cs,riemann,ok,passive_upwind=.true., &
        internal_energy_flux=expected_h)
      call require(ok,'independent final Riemann UEINT')
      call require(prediction(axis)%internal_energy_flux(idx(1),idx(2),idx(3))==expected_h, &
        'normal and transverse corrected stages retain actual Riemann UEINT')
    enddo
  enddo
  call check_boundary_family()
  call check_paired_flux_correction()
  print *, 'PASS: Godunov UEINT advection, pdivu, actual-held viscosity, gas identity and atomic refusal; residual=',maxres
contains
  subroutine check_paired_flux_correction()
    type(eb_geometry_3d) :: shape
    type(eb_integrated_flux_3d) :: raw,low,corrected,reference,aux,aux_low,kept_aux
    real(dp) :: u(nv,m,m,m),forcing(nv,m,m,m),rhs(nv,m,m,m),raw_rhs(nv,m,m,m),expected(nv,m,m,m)
    real(dp) :: h_rhs(m,m,m),kept_h(m,m,m),theta,reference_h,divergence,dt_local
    integer :: ii,jj,kk
    logical :: good
    dt_local=.001_dp
    call build_axis_plane_eb_geometry_3d(m,m,m,0._dp,1._dp,0._dp,2._dp,0._dp,3._dp, &
      'x',-.1_dp,shape,good)
    call require(good,'paired spatial correction geometry')
    u=0;u(1,:,:,:)=1;u(nv,:,:,:)=1;u(5,:,:,:)=200000;forcing=0
    allocate(raw%x(nv,0:m,m,m),raw%y(nv,m,0:m,m),raw%z(nv,m,m,0:m))
    raw%x=0;raw%y=0;raw%z=0;low=raw
    raw%x([1,nv],2,2,2)=1
    allocate(aux%x(1,0:m,m,m),aux%y(1,m,0:m,m),aux%z(1,m,m,0:m))
    aux%x=0;aux%y=0;aux%z=0;aux_low=aux
    aux%x(1,2,2,2)=13;aux_low%x(1,2,2,2)=5
    raw_rhs=0
    raw_rhs([1,nv],2,2,2)=-1/(shape%dx*dt_local)
    raw_rhs([1,nv],3,2,2)=1/(shape%dx*dt_local)
    h_rhs=-42
    h_rhs(2,2,2)=h_rhs(2,2,2)-13/(shape%dx*dt_local)
    h_rhs(3,2,2)=h_rhs(3,2,2)+13/(shape%dx*dt_local)
    rhs=raw_rhs;corrected=raw
    call correct_reactive_eb_spatial_hydro_flux_3d(shape,u,forcing,dt_local,rhs,corrected,low,good)
    call require(good,'legacy nonzero spatial correction')
    expected=rhs;reference=corrected;rhs=raw_rhs;corrected=raw
    call correct_reactive_eb_spatial_hydro_flux_3d(shape,u,forcing,dt_local,rhs,corrected,low,good, &
      internal_energy_rate=h_rhs,internal_energy_flux=aux,internal_energy_baseline=aux_low)
    call require(good,'paired nonzero spatial correction')
    call require(all(rhs==expected).and.all(corrected%x==reference%x).and.all(corrected%y==reference%y).and. &
      all(corrected%z==reference%z),'held spatial correction preserves exact gas limiter')
    theta=corrected%x(1,2,2,2)
    call require(theta>0.and.theta<1,'fixture activates shared face correction')
    reference_h=5+theta*(13-5)
    call require(abs(aux%x(1,2,2,2)-reference_h)<1e-14_dp,'same gas face coefficient scales independent held family')
    do kk=1,m
    do jj=1,m
    do ii=1,m
      divergence=(aux%x(1,ii-1,jj,kk)-aux%x(1,ii,jj,kk))/(shape%dx*dt_local)
      call require(abs(h_rhs(ii,jj,kk)-(divergence-42))<1e-10_dp, &
        'corrected held rate preserves pressure-work residual')
    enddo
    enddo
    enddo
    kept_aux=aux;kept_h=h_rhs;aux_low%x(1,m,m,m)=nan_value
    call correct_reactive_eb_spatial_hydro_flux_3d(shape,u,forcing,dt_local,rhs,corrected,low,good, &
      internal_energy_rate=h_rhs,internal_energy_flux=aux,internal_energy_baseline=aux_low)
    call require(.not.good.and.all(rhs==expected).and.all(h_rhs==kept_h), &
      'late invalid held correction baseline preserves both rates')
    call require(all(aux%x==kept_aux%x).and.all(aux%y==kept_aux%y).and.all(aux%z==kept_aux%z).and. &
      all(corrected%x==reference%x).and.all(corrected%y==reference%y).and.all(corrected%z==reference%z), &
      'late invalid held correction baseline preserves both face families')
  end subroutine check_paired_flux_correction
  subroutine check_boundary_family()
    type(eb_geometry_3d) :: extended
    type(eb_integrated_flux_3d) :: gasface,gasbase,holdface,holdbase,referenceface,referencebase
    real(dp) :: initial(nv,m,m,m),initial_t(m,m,m),h(m,m,m),forcing(nv,m,m,m),r(nv,m,m,m),reference(nv,m,m,m)
    real(dp) :: inlet(nv,6),inlet_t(6),inlet_h(6),b(nv),w(nv),rb(nv),rw(nv),hr(m,m,m),pw(m,m,m),rhsave(m,m,m)
    real(dp) :: before(m,m,m),sound_value,inventory_change,boundary_change
    integer :: ii,jj,kk,species_index,boundary_mode
    logical :: good
    character(len=24) :: bc(6)
    character(len=24), parameter :: modes(6)=[character(len=24):: &
      'outflow','periodic','slip_wall','no_slip_wall','isothermal_wall','inflow']
    call build_axis_plane_eb_geometry_3d(m+12,m+12,m+12,-1.5_dp,2.5_dp,-3._dp,5._dp,-4.5_dp,7.5_dp, &
      'x',-2._dp,extended,good)
    call require(good,'boundary extended regular geometry')
    do kk=1,m
    do jj=1,m
    do ii=1,m
      q=0;q(1)=1;q(2:4)=[4._dp,2._dp,-1._dp]-.1_dp*[ii,jj,kk];q(5)=150000
      do species_index=1,size(species)
        if(species(species_index)%name=='N2') q(species_index+5)=1
      enddo
      call reactive_primitive_to_conserved(species,q,initial(:,ii,jj,kk),initial_t(ii,jj,kk),sound_value,good)
      call require(good,'boundary actual input')
      h(ii,jj,kk)=initial(5,ii,jj,kk)-.5_dp*sum(initial(2:4,ii,jj,kk)**2)+37*ii
    enddo
    enddo
    enddo
    forcing=0
    do ii=1,6
      inlet(:,ii)=initial(:,1,1,1);inlet_t(ii)=initial_t(1,1,1);inlet_h(ii)=h(1,1,1)+713
    enddo
    do boundary_mode=1,size(modes)
      bc=modes(boundary_mode)
      if(boundary_mode==6) then
        bc='outflow';bc(1)='inflow'
      endif
      call compute_sdc_eb_boundary_hydro_3d(species,initial,initial_t,extended,1e-7_dp,forcing,r,gasface,b,w,good,bc, &
        shock_flattening=.false.,inflow_state=inlet,inflow_temperature=inlet_t,baseline_flux=gasbase)
      if(.not.good) print *, 'Rejected boundary mode: ',trim(modes(boundary_mode))
      call require(good,'legacy physical boundary spatial family')
      reference=r;referenceface=gasface;referencebase=gasbase;rb=b;rw=w
      call compute_sdc_eb_boundary_hydro_3d(species,initial,initial_t,extended,1e-7_dp,forcing,r,gasface,b,w,good,bc, &
        shock_flattening=.false.,inflow_state=inlet,inflow_temperature=inlet_t,baseline_flux=gasbase, &
        internal_energy_state=h,internal_energy_rate=hr,internal_energy_integrated_flux=holdface, &
        internal_energy_baseline_flux=holdbase,internal_energy_pressure_work=pw,inflow_internal_energy=inlet_h)
      call require(good,'actual held physical boundary spatial family')
      call require(all(r==reference).and.all(b==rb).and.all(w==rw), &
        'physical boundary held request preserves gas rates and ledgers')
      call require(all(gasface%x==referenceface%x).and.all(gasface%y==referenceface%y).and. &
        all(gasface%z==referenceface%z),'physical boundary held request preserves gas faces')
      call require(all(gasbase%x==referencebase%x).and.all(gasbase%y==referencebase%y).and. &
        all(gasbase%z==referencebase%z),'physical boundary held request preserves gas baseline')
      inventory_change=1e-7_dp*extended%dx*extended%dy*extended%dz*sum(hr+pw)
      boundary_change=extended%dy*extended%dz*sum(holdface%x(1,0,:,:)-holdface%x(1,m,:,:))+ &
        extended%dx*extended%dz*sum(holdface%y(1,:,0,:)-holdface%y(1,:,m,:))+ &
        extended%dx*extended%dy*sum(holdface%z(1,:,:,0)-holdface%z(1,:,:,m))
      call require(abs(inventory_change-boundary_change)<1e-15_dp,'physical boundary held advective inventory')
      if(boundary_mode==2) call require(abs(boundary_change)<1e-15_dp,'periodic held face cancellation')
      rhsave=hr;before=pw;inlet_h(1)=nan_value
      call compute_sdc_eb_boundary_hydro_3d(species,initial,initial_t,extended,1e-7_dp,forcing,r,gasface,b,w,good,bc, &
        shock_flattening=.false.,inflow_state=inlet,inflow_temperature=inlet_t,baseline_flux=gasbase, &
        internal_energy_state=h,internal_energy_rate=hr,internal_energy_integrated_flux=holdface, &
        internal_energy_baseline_flux=holdbase,internal_energy_pressure_work=pw,inflow_internal_energy=inlet_h)
      call require(.not.good.and.all(hr==rhsave).and.all(pw==before).and.all(r==reference), &
        'invalid independent inlet held energy preserves complete spatial family')
      inlet_h(1)=h(1,1,1)+713
    enddo
  end subroutine check_boundary_family
  subroutine check_weights()
    type(eb_integrated_flux_3d) :: high,high_h,low,low_h,expected,expected_h
    real(dp), allocatable :: high_rate(:,:,:,:),high_rhs(:,:,:),high_work(:,:,:),saved_b(:),saved_w(:)
    real(dp) :: weights(3)=[0._dp,.375_dp,1._dp],weight,advection,residual
    integer :: test
    high=flux;high_h=hflux;low=base;low_h=hbase
    high_rate=rate;high_rhs=hrate;high_work=pwork;saved_b=boundary;saved_w=wall
    do test=1,size(weights)
      weight=weights(test)
      call compute_sdc_eb_hydro_3d(species,g,ghost,temp,vf,1e-7_dp,flat,source,rate,flux,boundary,wall,ok, &
        difmag=coefficient,baseline_flux=base,internal_energy_state=held,internal_energy_rate=hrate, &
        internal_energy_integrated_flux=hflux,internal_energy_baseline_flux=hbase, &
        internal_energy_pressure_work=pwork,spatial_flux_weight=weight)
      call require(ok,'native paired Godunov interval face weight')
      call require(all(pwork==high_work),'interval weighting preserves actual native pressure work')
      call require(all(wall==saved_w),'interval weighting preserves physical embedded-wall volume work')
      if(weight==1) then
        call require(all(rate==high_rate).and.all(hrate==high_rhs).and.all(boundary==saved_b), &
          'weight-one preserves gas/UEINT rates and volume ledgers exactly')
        expected=high;expected_h=high_h
      else
        expected=low;expected_h=low_h
        if(weight>0) then
          expected%x=(1-weight)*low%x+weight*high%x
          expected%y=(1-weight)*low%y+weight*high%y
          expected%z=(1-weight)*low%z+weight*high%z
          expected_h%x=(1-weight)*low_h%x+weight*high_h%x
          expected_h%y=(1-weight)*low_h%y+weight*high_h%y
          expected_h%z=(1-weight)*low_h%z+weight*high_h%z
        endif
      endif
      call require(all(flux%x==expected%x).and.all(flux%y==expected%y).and.all(flux%z==expected%z).and. &
        all(hflux%x==expected_h%x).and.all(hflux%y==expected_h%y).and.all(hflux%z==expected_h%z), &
        'gas and UEINT use the same coefficient and exact endpoint/convex integrated faces')
      ! Independent scalar conservation, before StateRedist. Pressure work is
      ! a native volume term, so add it back to test only face advection.
      inventory=1e-7_dp*cellvol*sum(physical%volume_fraction*(hrate+pwork))
      advection=physical%dy*physical%dz*sum(physical%x_face_fraction(0,:,:)*hflux%x(1,0,:,:)- &
        physical%x_face_fraction(m,:,:)*hflux%x(1,m,:,:))+ &
        physical%dx*physical%dz*sum(physical%y_face_fraction(:,0,:)*hflux%y(1,:,0,:)- &
        physical%y_face_fraction(:,m,:)*hflux%y(1,:,m,:))+ &
        physical%dx*physical%dy*sum(physical%z_face_fraction(:,:,0)*hflux%z(1,:,:,0)- &
        physical%z_face_fraction(:,:,m)*hflux%z(1,:,:,m))
      residual=abs(inventory-advection);maxres=max(maxres,residual)
      call require(residual<1e-15_dp,'weighted independent UEINT face inventory and fixed pressure-work residual')
    enddo
    do test=1,3
      weight=nan_value
      if(test==2) weight=-.1_dp
      if(test==3) weight=1.1_dp
      call compute_sdc_eb_hydro_3d(species,g,ghost,temp,vf,1e-7_dp,flat,source,rate,flux,boundary,wall,ok, &
        difmag=coefficient,baseline_flux=base,internal_energy_state=held,internal_energy_rate=hrate, &
        internal_energy_integrated_flux=hflux,internal_energy_baseline_flux=hbase, &
        internal_energy_pressure_work=pwork,spatial_flux_weight=weight)
      call require(.not.ok,'invalid native interval coefficient refuses')
      call require(all(rate==high_rate).and.all(hrate==high_rhs).and.all(pwork==high_work).and. &
        all(boundary==saved_b).and.all(wall==saved_w),'invalid weight retains prior rates and native volume terms')
      call require(all(flux%x==high%x).and.all(flux%y==high%y).and.all(flux%z==high%z).and. &
        all(hflux%x==high_h%x).and.all(hflux%y==high_h%y).and.all(hflux%z==high_h%z), &
        'invalid weight retains both prior face families')
    enddo
  end subroutine check_weights
  subroutine evaluate_held()
    call compute_sdc_eb_hydro_3d(species,g,ghost,temp,vf,1e-7_dp,flat,source,rate,flux,boundary,wall,ok, &
      difmag=coefficient,baseline_flux=base,internal_energy_state=held,internal_energy_rate=hrate, &
      internal_energy_integrated_flux=hflux,internal_energy_baseline_flux=hbase,internal_energy_pressure_work=pwork)
  end subroutine evaluate_held
  subroutine check_gas_identity()
    call require(all(rate==kept_rate).and.all(boundary==kept_boundary).and.all(wall==kept_wall), &
      'paired request preserves gas rate and ledgers exactly')
    call require(all(flux%x==kept_flux%x).and.all(flux%y==kept_flux%y).and.all(flux%z==kept_flux%z), &
      'paired request preserves gas faces exactly')
    call require(all(base%x==kept_base%x).and.all(base%y==kept_base%y).and.all(base%z==kept_base%z), &
      'paired request preserves gas baseline exactly')
  end subroutine check_gas_identity
  subroutine check_all_unchanged()
    call check_gas_identity()
    call require(all(hrate==kept_hrate).and.all(pwork==kept_pwork),'failure preserves held rate and pressure work')
    call require(all(hflux%x==kept_hflux%x).and.all(hflux%y==kept_hflux%y).and.all(hflux%z==kept_hflux%z), &
      'failure preserves held faces')
    call require(all(hbase%x==kept_hbase%x).and.all(hbase%y==kept_hbase%y).and.all(hbase%z==kept_hbase%z), &
      'failure preserves held baseline')
  end subroutine check_all_unchanged
  subroutine check_pressure_work(cut)
    integer, intent(in) :: cut
    integer, parameter :: d=6
    type(eb_geometry_3d) :: geometry
    type(eb_plm_face_prediction_3d) :: faces(3)
    real(dp) :: phi(0:d,0:d,0:d),pressure(d,d,d),work(d-2,d-2,d-2),saved(d-2,d-2,d-2)
    real(dp) :: analytic,previous
    integer :: ii,jj,kk,a,lo(3),hi(3)
    logical :: good
    if(cut==0) then
      call build_axis_plane_eb_geometry_3d(d,d,d,0._dp,1._dp,0._dp,2._dp,0._dp,3._dp, &
        'x',-.1_dp,geometry,good)
    else
      do kk=0,d
      do jj=0,d
      do ii=0,d
        phi(ii,jj,kk)=real(ii,dp)/d-.31_dp*real(jj,dp)/d+.17_dp*real(kk,dp)/d-.41_dp
        if(cut==2) phi(ii,jj,kk)=sum(([real(ii,dp)/d,real(jj,dp)/d,real(kk,dp)/d]- &
          [.45_dp,.55_dp,.5_dp])**2)-.4_dp**2
      enddo
      enddo
      enddo
      call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,2._dp,3._dp],geometry,good)
    endif
    call require(good,'manufactured pressure-work geometry')
    pressure=170003._dp
    do a=1,3
      lo=1;lo(a)=0;hi=d
      allocate(faces(a)%pressure(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        faces(a)%velocity(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      do kk=lo(3),hi(3)
      do jj=lo(2),hi(2)
      do ii=lo(1),hi(1)
        faces(a)%pressure(ii,jj,kk)=150001._dp
        faces(a)%velocity(ii,jj,kk)=.1_dp*a*ii+.2_dp*a*jj-.3_dp*a*kk+.013_dp*ii*jj*kk
      enddo
      enddo
      enddo
    enddo
    call compute_eb_godunov_pressure_work_3d(geometry,faces,pressure,work,good)
    call require(good,'manufactured pressure work')
    if(cut==0) then
      do kk=2,d-1
      do jj=2,d-1
      do ii=2,d-1
        analytic=150001._dp*((.1_dp+.013_dp*jj*kk)/geometry%dx+ &
          (.4_dp+.013_dp*ii*kk)/geometry%dy+(-.9_dp+.013_dp*ii*jj)/geometry%dz)
        call require(abs(work(ii-1,jj-1,kk-1)-analytic)<2e-9_dp,'independent regular pressure work')
      enddo
      enddo
      enddo
    else
      call require(any(geometry%x_face_fraction>0.and.geometry%x_face_fraction<1), &
        'nontrivial cut centroid pressure-work fixture')
      call require(any(geometry%volume_fraction(2:d-1,2:d-1,2:d-1)==0),'covered pressure-work fixture')
    endif
    if(dump_unit/=0) then
      write(dump_unit,*) d,d,d
      write(dump_unit,'(3(es26.17,1x))') geometry%dx,geometry%dy,geometry%dz
      write(dump_unit,*) geometry%cell_type
      write(dump_unit,'(es26.17)') geometry%volume_fraction,pressure
      write(dump_unit,'(es26.17)') geometry%x_face_fraction,geometry%y_face_fraction,geometry%z_face_fraction
      write(dump_unit,'(es26.17)') geometry%x_face_centroid_y,geometry%x_face_centroid_z, &
        geometry%y_face_centroid_x,geometry%y_face_centroid_z,geometry%z_face_centroid_x,geometry%z_face_centroid_y
      do a=1,3
        write(dump_unit,'(es26.17)') faces(a)%pressure,faces(a)%velocity
      enddo
      write(dump_unit,'(es26.17)') work
    endif
    saved=work;previous=faces(3)%velocity(d-1,d-1,d-1)
    faces(3)%velocity(d-1,d-1,d-1)=nan_value
    call compute_eb_godunov_pressure_work_3d(geometry,faces,pressure,work,good)
    call require(.not.good.and.all(work==saved),'late pressure-work refusal is atomic')
    faces(3)%velocity(d-1,d-1,d-1)=previous
    if(cut==0) then
      pressure=nan_value
      call compute_eb_godunov_pressure_work_3d(geometry,faces,pressure,work,good)
      call require(good.and.all(work==saved),'regular stencil does not consume cell pressure')
    endif
  end subroutine check_pressure_work
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL: ',label
    error stop 1
  end subroutine require
end program test_eb_godunov_internal_energy_3d
