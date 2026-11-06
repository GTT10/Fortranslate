program test_eb_mol_hydro_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  use mixture_thermo_mod, only: mixture_mass_properties
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d,limit_mol_centroid_species_fluxes_3d, &
    eb_mol_diagnostics_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: g,physical
  type(eb_integrated_flux_3d) :: flux,saved,limiter_raw,limiter_flux,aux_flux,kept_aux_flux,aux_raw,patch_raw
  type(eb_mol_diagnostics_3d) :: diagnostic,kept_diagnostic
  real(dp) :: balance(3),spacing(3),total_boundary(12),total_wall(12)
  real(dp), allocatable :: irrelevant_held(:,:,:)
  type(rotating_frame_3d), allocatable :: frame
  real(dp), allocatable :: potential,q_back(:),ordinary(:),zero_energy(:)
  real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),q(:),rates(:,:,:,:),boundary(:),wall(:),inventory(:)
  real(dp), allocatable :: aux_rate(:,:,:),kept_rate(:,:,:,:),kept_aux_rate(:,:,:)
  real(dp), allocatable :: limiter_state(:,:,:,:),limiter_alpha(:,:,:,:)
  real(dp) :: phi(0:10,0:10,0:10),sound,t,residual,dt,position(3),t_back,c_back,omega,donor_mass,face_mass
  real(dp) :: aux_boundary,aux_inventory,held_value,original_energy
  integer :: i,j,k,c,nv,mode,order,rotation,solver,cases,face_cell(3),left(3),right(3)
  logical :: ok,found
  character(len=64) :: argument
  omega=400._dp
  if(command_argument_count()>0) then
    call get_command_argument(1,argument)
    read(argument,*) omega
  end if
  call load_h2o2_elementary_thermo(species,ok)
  if(.not.ok) error stop 'thermo'
  nv=reactive_nvar(size(species));dt=1e-7_dp
  call check_held_input()
  call check_spatial_halo()
  allocate(state(nv,10,10,10),temp(10,10,10),q(reactive_nprim(size(species))), &
    rates(nv,4,4,4),boundary(nv),wall(nv),inventory(nv),q_back(reactive_nprim(size(species))),ordinary(nv),zero_energy(nv))
  do k=0,10
  do j=0,10
  do i=0,10
    phi(i,j,k)=(-.75_dp+.25_dp*i)-.3_dp*(-.75_dp+.25_dp*j)-.44_dp
  end do
  end do
  end do
  call build_level_set_eb_geometry_3d(phi,[-.75_dp,-.75_dp,-.75_dp],[1.75_dp,1.75_dp,1.75_dp],g,ok)
  if(.not.ok) error stop 'geometry'
  call extract_eb_geometry_3d(g,[4,4,4],[7,7,7],physical,ok)
  if(.not.ok) error stop 'physical geometry'
  allocate(aux_rate(4,4,4),kept_rate(nv,4,4,4))
  cases=0
  do rotation=0,1
  if(rotation==1) then
    allocate(frame,potential)
    frame%angular_speed=omega;frame%axis=3;frame%axis_point=[.13_dp,-.07_dp,.2_dp]
  end if
  do mode=1,2
    do k=1,10
    do j=1,10
    do i=1,10
      q=0;q(1)=1;q(5)=150000
      if(mode==2) then
        q(1)=1+.01_dp*i;q(2:4)=[.1_dp*i,.2_dp*j,-.1_dp*k];q(5)=150000+100*j
      end if
      do c=1,size(species)
        if(species(c)%name=='N2') q(c+5)=1
      end do
      if(allocated(frame)) then
        position=[g%x_lower+(i-.5_dp)*g%dx,g%y_lower+(j-.5_dp)*g%dy,g%z_lower+(k-.5_dp)*g%dz]
        call rotating_potential_3d(frame,position,potential,ok)
        if(.not.ok) error stop 'cell potential'
      end if
      call reactive_primitive_to_conserved(species,q,ordinary,t,sound,ok)
      if(.not.ok) error stop 'ordinary state'
      call reactive_primitive_to_conserved(species,q,zero_energy,t,sound,ok,0._dp)
      if(.not.ok.or.any(zero_energy/=ordinary)) error stop 'zero potential identity'
      call reactive_primitive_to_conserved(species,q,state(:,i,j,k),t,sound,ok,potential)
      if(.not.ok) error stop 'initial state'
      if(allocated(potential)) then
        if(state(5,i,j,k)/=ordinary(5)-q(1)*potential) error stop 'rotating energy definition'
        if(any(state(1:4,i,j,k)/=ordinary(1:4)).or.any(state(6:,i,j,k)/=ordinary(6:))) &
          error stop 'rotating unaffected components'
      end if
      call reactive_conserved_to_primitive(species,state(:,i,j,k),t,q_back,t_back,c_back,ok,potential)
      if(.not.ok) error stop 'rotating EOS recovery'
      if(maxval(abs(q_back-q)/max(1._dp,abs(q)))>1e-12_dp.or.abs(t_back/t-1)>1e-12_dp.or. &
        abs(c_back/sound-1)>1e-12_dp) error stop 'EOS roundtrip'
      temp(i,j,k)=t
      if(g%volume_fraction(i,j,k)<=0) state(:,i,j,k)=ieee_value(0._dp,ieee_quiet_nan)
    end do
    end do
    end do
    if(rotation==1.and.minval(state(5,:,:,:),mask=g%volume_fraction>0)>=0) &
      error stop 'negative rotating energy fixture'
    found=.false.
    do k=1,10
    do j=1,10
    do i=1,10
      if(g%volume_fraction(i,j,k)<=0) cycle
      found=.true.;exit
    enddo
      if(found) exit
    enddo
      if(found) exit
    enddo
    if(.not.found) error stop 'retained EOS physical cell fixture'
    position=[g%x_lower+(i-.5_dp)*g%dx,g%y_lower+(j-.5_dp)*g%dy,g%z_lower+(k-.5_dp)*g%dz]
    if(allocated(frame)) call rotating_potential_3d(frame,position,potential,ok)
    original_energy=state(5,i,j,k)
    held_value=original_energy-.5_dp*sum(state(2:4,i,j,k)**2)/state(1,i,j,k)
    if(allocated(potential)) held_value=held_value+state(1,i,j,k)*potential
    held_value=held_value+1000._dp
    call recover_reactive_eb_cell_3d(species,g,[i,j,k],state(:,i,j,k),temp(i,j,k), &
      q_back,t_back,c_back,ok,frame,held_value)
    if(.not.ok.or.t_back<=temp(i,j,k).or.state(5,i,j,k)/=original_energy) &
      error stop 'retained EB EOS uses auxiliary without changing total energy'
    do solver=0,1
    do order=1,2
      call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok,lax_friedrichs=solver==1,frame=frame)
      if(.not.ok) error stop 'MOL spatial rate'
      kept_rate=rates;saved=flux
      call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok, &
        lax_friedrichs=solver==1,frame=frame,diagnostics=diagnostic)
      if(.not.ok) error stop 'diagnostic spatial operator'
      if(any(rates/=kept_rate).or.any(flux%x/=saved%x).or.any(flux%y/=saved%y).or.any(flux%z/=saved%z)) &
        error stop 'diagnostic preserves operator exactly'
      if(any(diagnostic%rate/=rates)) error stop 'diagnostic actual rates'
      if(any(dt*diagnostic%centroid_flux%x/=flux%x).or.any(dt*diagnostic%centroid_flux%y/=flux%y).or. &
        any(dt*diagnostic%centroid_flux%z/=flux%z)) error stop 'diagnostic instantaneous flux units'
      spacing=[physical%dx,physical%dy,physical%dz]
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical%volume_fraction(i,j,k)<=0) cycle
        do c=2,4
          balance(1)=(physical%x_face_fraction(i-1,j,k)*diagnostic%centroid_flux%x(c,i-1,j,k)- &
            physical%x_face_fraction(i,j,k)*diagnostic%centroid_flux%x(c,i,j,k))/spacing(1)
          balance(2)=(physical%y_face_fraction(i,j-1,k)*diagnostic%centroid_flux%y(c,i,j-1,k)- &
            physical%y_face_fraction(i,j,k)*diagnostic%centroid_flux%y(c,i,j,k))/spacing(2)
          balance(3)=(physical%z_face_fraction(i,j,k-1)*diagnostic%centroid_flux%z(c,i,j,k-1)- &
            physical%z_face_fraction(i,j,k)*diagnostic%centroid_flux%z(c,i,j,k))/spacing(3)
          residual=sum(balance)/physical%volume_fraction(i,j,k)+diagnostic%wall_rate(c,i,j,k)-rates(c,i,j,k)
          if(abs(residual)>1e-10_dp*max(1._dp,sum(abs(balance))/physical%volume_fraction(i,j,k), &
            abs(diagnostic%wall_rate(c,i,j,k)))) error stop 'diagnostic directional momentum closure'
        enddo
      enddo
      enddo
      enddo
      call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok, &
        lax_friedrichs=solver==1,frame=frame,internal_energy_rate=aux_rate,internal_energy_integrated_flux=aux_flux)
      if(.not.ok) error stop 'MOL auxiliary spatial rate'
      if(any(rates/=kept_rate).or.any(flux%x/=saved%x).or.any(flux%y/=saved%y).or.any(flux%z/=saved%z)) &
        error stop 'auxiliary output preserves conserved operator'
      ! A deliberately different held scalar must not alter either face or wall
      ! pressure in explicit total-energy mode; auxiliary transport still closes.
      if(.not.allocated(irrelevant_held)) allocate(irrelevant_held(10,10,10))
      irrelevant_held=1.e6_dp
      kept_aux_flux=aux_flux;kept_aux_rate=aux_rate
      total_boundary=boundary;total_wall=wall
      call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok, &
        lax_friedrichs=solver==1,frame=frame,internal_energy_rate=aux_rate, &
        internal_energy_integrated_flux=aux_flux,internal_energy_density=irrelevant_held, &
        primitive_from_total_energy=.true.)
      if(.not.ok) error stop 'explicit total-energy MOL admission'
      if(any(rates/=kept_rate).or.any(flux%x/=saved%x).or.any(flux%y/=saved%y).or.any(flux%z/=saved%z)) &
        error stop 'total-energy face and wall operator identity'
      if(any(boundary/=total_boundary).or.any(wall/=total_wall)) error stop 'total-energy wall budget identity'
      if(any(aux_rate/=kept_aux_rate).or.any(aux_flux%x/=kept_aux_flux%x).or. &
        any(aux_flux%y/=kept_aux_flux%y).or.any(aux_flux%z/=kept_aux_flux%z)) &
        error stop 'total-energy mode preserves auxiliary transport'
      if(solver==1) then
        kept_aux_flux=aux_flux;kept_aux_rate=aux_rate
        call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok, &
          lax_friedrichs=.true.,frame=frame,raw_face_flux=patch_raw,limiter_alpha=limiter_alpha, &
          internal_energy_rate=aux_rate,internal_energy_integrated_flux=aux_flux, &
          internal_energy_raw_face_flux=aux_raw)
        if(.not.ok) error stop 'MOL auxiliary raw flux'
        if(any(aux_rate/=kept_aux_rate).or.any(aux_flux%x/=kept_aux_flux%x).or. &
          any(aux_flux%y/=kept_aux_flux%y).or.any(aux_flux%z/=kept_aux_flux%z)) &
          error stop 'raw auxiliary output preserves accepted operator'
        if(any(shape(aux_raw%x)/=shape(aux_flux%x)).or.any(shape(aux_raw%y)/=shape(aux_flux%y)).or. &
          any(shape(aux_raw%z)/=shape(aux_flux%z))) error stop 'raw auxiliary face layout'
      endif
      aux_inventory=dt*physical%dx*physical%dy*physical%dz*sum(physical%volume_fraction*aux_rate)
      aux_boundary=physical%dy*physical%dz*sum( &
        physical%x_face_fraction(0,:,:)*aux_flux%x(1,0,:,:)- &
        physical%x_face_fraction(4,:,:)*aux_flux%x(1,4,:,:))+ &
        physical%dx*physical%dz*sum(physical%y_face_fraction(:,0,:)*aux_flux%y(1,:,0,:)- &
        physical%y_face_fraction(:,4,:)*aux_flux%y(1,:,4,:))+ &
        physical%dx*physical%dy*sum(physical%z_face_fraction(:,:,0)*aux_flux%z(1,:,:,0)- &
        physical%z_face_fraction(:,:,4)*aux_flux%z(1,:,:,4))
      if(abs(aux_inventory-aux_boundary)>1e-15_dp) error stop 'auxiliary volume-weighted flux conservation'
      do c=1,nv
        inventory(c)=dt*physical%dx*physical%dy*physical%dz*sum(physical%volume_fraction*rates(c,:,:,:))
      end do
      residual=maxval(abs(inventory-boundary-wall))
      if(residual>1e-15_dp) then
        print *, 'FAILED conservation: omega, rotation, mode, solver, order, residual:', &
          omega,rotation,mode,solver,order,residual
        error stop 'boundary and wall conservation'
      end if
      if(mode==1.and.(rotation==0.or.solver==0).and.maxval(abs(rates))>1e-7_dp) error stop 'stationary inclined wall balance'
      if(mode==2.and.maxval(abs(rates))<1) error stop 'nontrivial spatial rate'
      cases=cases+1
      print *, 'rotation, mode, solver, order, conservation residual:',rotation,mode,solver,order,residual
    end do
    end do
  end do
  end do
  ! A centroid correction that would remove a species absent from its donor
  ! must scale the entire shared face vector. A positive donor limits the
  ! correction by its cell mass without changing the baseline face flux.
  found=.false.;face_cell=0
  face_search: do k=1,physical%nz
  do j=1,physical%ny
  do i=1,physical%nx-1
    if(physical%volume_fraction(i,j,k)<=0.or.physical%volume_fraction(i+1,j,k)<=0) cycle
    if(physical%x_face_fraction(i,j,k)<=0) cycle
    face_cell=[i,j,k];found=.true.;exit face_search
  end do
  end do
  end do face_search
  if(.not.found) error stop 'centroid limiter face fixture'
  left=face_cell;right=left;right(1)=right(1)+1
  allocate(limiter_state(nv,physical%nx,physical%ny,physical%nz))
  limiter_state=0;limiter_state(1,:,:,:)=1
  do c=1,size(species)
    if(species(c)%name=='N2') limiter_state(c+5,:,:,:)=1
  end do
  allocate(limiter_raw%x(nv+1,0:physical%nx+2,physical%ny+2,physical%nz+2), &
    limiter_raw%y(nv+1,physical%nx+2,0:physical%ny+2,physical%nz+2), &
    limiter_raw%z(nv+1,physical%nx+2,physical%ny+2,0:physical%nz+2))
  allocate(limiter_flux%x(nv+1,0:physical%nx,physical%ny,physical%nz), &
    limiter_flux%y(nv+1,physical%nx,0:physical%ny,physical%nz), &
    limiter_flux%z(nv+1,physical%nx,physical%ny,0:physical%nz))
  limiter_raw%x=0;limiter_raw%y=0;limiter_raw%z=0
  limiter_flux%x=0;limiter_flux%y=0;limiter_flux%z=0
  limiter_flux%x([1,5,6,nv+1],i,j,k)=[1._dp,7._dp,1._dp,11._dp]
  call limit_mol_centroid_species_fluxes_3d(physical,limiter_state,limiter_raw,limiter_flux,dt,ok)
  if(.not.ok.or.any(limiter_flux%x(:,i,j,k)/=0)) error stop 'empty donor centroid correction'
  if(lbound(limiter_flux%x,2)/=0.or.lbound(limiter_flux%y,3)/=0.or. &
    lbound(limiter_flux%z,4)/=0) error stop 'centroid limiter preserves face indexing'
  face_mass=dt*physical%dy*physical%dz*physical%x_face_fraction(i,j,k)
  donor_mass=.5_dp*face_mass
  limiter_state(6,i,j,k)=donor_mass/(physical%volume_fraction(i,j,k)* &
    physical%dx*physical%dy*physical%dz)
  limiter_flux%x([1,5,6,nv+1],i,j,k)=[1._dp,7._dp,1._dp,11._dp]
  call limit_mol_centroid_species_fluxes_3d(physical,limiter_state,limiter_raw,limiter_flux,dt,ok)
  if(.not.ok.or.abs(limiter_flux%x(6,i,j,k)-.5_dp)>1e-13_dp) error stop 'finite donor centroid bound'
  if(abs(limiter_flux%x(1,i,j,k)-limiter_flux%x(6,i,j,k))>1e-15_dp.or. &
    abs(limiter_flux%x(5,i,j,k)-7*limiter_flux%x(6,i,j,k))>1e-14_dp) &
    error stop 'centroid correction scales entire conserved face'
  if(abs(limiter_flux%x(nv+1,i,j,k)-11*limiter_flux%x(6,i,j,k))>1e-14_dp) &
    error stop 'centroid correction scales auxiliary flux with species allowance'
  saved=limiter_flux
  limiter_state(6,right(1),right(2),right(3))=ieee_value(0._dp,ieee_quiet_nan)
  call limit_mol_centroid_species_fluxes_3d(physical,limiter_state,limiter_raw,limiter_flux,dt,ok)
  if(ok.or.any(limiter_flux%x/=saved%x).or.any(limiter_flux%y/=saved%y).or. &
    any(limiter_flux%z/=saved%z)) error stop 'centroid limiter rollback'
  limiter_state(6,right(1),right(2),right(3))=0
  solver=0
  saved=flux;rates=17;boundary=19;wall=23
  kept_aux_flux=aux_flux;kept_aux_rate=aux_rate
  kept_diagnostic=diagnostic
  state(1,10,10,10)=ieee_value(0._dp,ieee_quiet_nan)
  call compute_mol_eb_hydro_3d(species,g,state,temp,dt,2,rates,flux,boundary,wall,ok,lax_friedrichs=solver==1,frame=frame, &
    internal_energy_rate=aux_rate,internal_energy_integrated_flux=aux_flux,diagnostics=diagnostic)
  if(any(diagnostic%rate/=kept_diagnostic%rate).or. &
    any(diagnostic%pressure%x/=kept_diagnostic%pressure%x)) error stop "diagnostic rollback"
  if(any(aux_rate/=kept_aux_rate).or.any(aux_flux%x/=kept_aux_flux%x).or. &
    any(aux_flux%y/=kept_aux_flux%y).or.any(aux_flux%z/=kept_aux_flux%z)) error stop 'auxiliary rollback'
  if(ok.or.any(rates/=17).or.any(boundary/=19).or.any(wall/=23)) error stop 'rate rollback'
  if(any(flux%x/=saved%x).or.any(flux%y/=saved%y).or.any(flux%z/=saved%z)) error stop 'flux rollback'
  state(1,10,10,10)=q(1)
  frame%axis=0
  call compute_mol_eb_hydro_3d(species,g,state,temp,dt,2,rates,flux,boundary,wall,ok,frame=frame)
  if(ok.or.any(rates/=17).or.any(boundary/=19).or.any(wall/=23)) error stop 'invalid frame rollback'
  if(any(flux%x/=saved%x).or.any(flux%y/=saved%y).or.any(flux%z/=saved%z)) error stop 'invalid frame flux rollback'
  call reactive_primitive_to_conserved(species,q,ordinary,t,sound,ok,-1._dp)
  if(ok.or.any(ordinary/=0).or.t/=0.or.sound/=0) error stop 'negative potential rejection'
  call reactive_conserved_to_primitive(species,state(:,10,10,10),temp(10,10,10),q_back,t,sound,ok, &
    ieee_value(0._dp,ieee_quiet_nan))
  if(ok.or.any(q_back/=0).or.t/=0.or.sound/=0) error stop 'NaN potential rejection'
  print *, 'PASS: MOL inclined-wall spatial rate, EOS roundtrip and rollback; cases:',cases
contains
  subroutine check_held_input()
    type(eb_geometry_3d) :: grid
    type(rotating_frame_3d) :: spin
    type(eb_integrated_flux_3d) :: f,hf,kept_f,kept_hf
    real(dp) :: field_phi(0:10,0:10,0:10),cell(12),q(12),y(7),rho,t,cs,p,position(3),phi
    real(dp) :: mw,r,cp,cv,gamma,h,e,entropy
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),stored(:,:,:),rates(:,:,:,:),aux(:,:,:)
    real(dp) :: b(12),w(12),expected(12),kept_rates(12,4,4,4),kept_aux(4,4,4)
    integer :: i,j,k
    logical :: good
    field_phi=1
    call build_level_set_eb_geometry_3d(field_phi,[-.75_dp,-.75_dp,-.75_dp],[1.75_dp,1.75_dp,1.75_dp],grid,good)
    if(.not.good) error stop 'held MOL regular geometry'
    allocate(gas(12,10,10,10),temp(10,10,10),stored(10,10,10),rates(12,4,4,4),aux(4,4,4))
    rho=.31_dp;t=650._dp;y=0;y(7)=1
    call mixture_mass_properties(species,y,t,mw,r,cp,cv,gamma,h,e,entropy,good)
    if(.not.good) error stop 'held MOL independent NASA properties'
    p=rho*universal_gas_constant/species(7)%molecular_weight*t
    q=0;q(1)=rho;q(5)=p;q(12)=1
    call reactive_primitive_to_conserved(species,q,cell,t,cs,good)
    if(.not.good) error stop 'held MOL seed'
    spin%angular_speed=1.e10_dp;spin%axis=3
    do k=1,10
    do j=1,10
    do i=1,10
      position=[grid%x_lower,grid%y_lower,grid%z_lower]+ &
        (real([i,j,k],dp)-.5_dp)*[grid%dx,grid%dy,grid%dz]
      call rotating_potential_3d(spin,position,phi,good)
      if(.not.good) error stop 'held MOL potential'
      gas(:,i,j,k)=cell;gas(5,i,j,k)=cell(5)-rho*phi
      temp(i,j,k)=t;stored(i,j,k)=rho*e
    enddo
    enddo
    enddo
    call compute_mol_eb_hydro_3d(species,grid,gas,temp,1.e-25_dp,2,rates,f,b,w,good,frame=spin, &
      internal_energy_rate=aux,internal_energy_integrated_flux=hf,internal_energy_density=stored)
    if(.not.good) error stop 'actual held scalar admitted by rotating MOL kernel'
    if(maxval(abs(f%x(2,:,:,:)/1.e-25_dp-p))>3.e-10_dp*p) error stop 'independent MOL x pressure'
    if(maxval(abs(f%y(3,:,:,:)/1.e-25_dp-p))>3.e-10_dp*p) error stop 'independent MOL y pressure'
    if(maxval(abs(f%z(4,:,:,:)/1.e-25_dp-p))>3.e-10_dp*p) error stop 'independent MOL z pressure'
    if(maxval(abs(rates))>3.e-10_dp*p/grid%dx.or.any(aux/=0)) error stop 'uniform MOL gas/scalar equilibrium'
    kept_f=f;kept_hf=hf
    rates=17;aux=29;b=19;w=23;kept_rates=rates;kept_aux=aux
    stored(10,10,10)=1.e30_dp
    call compute_mol_eb_hydro_3d(species,grid,gas,temp,1.e-25_dp,2,rates,f,b,w,good,frame=spin, &
      internal_energy_rate=aux,internal_energy_integrated_flux=hf,internal_energy_density=stored)
    if(good.or.any(rates/=kept_rates).or.any(aux/=kept_aux).or.any(b/=19).or.any(w/=23)) &
      error stop 'late MOL held EOS refusal rollback'
    if(any(f%x/=kept_f%x).or.any(f%y/=kept_f%y).or.any(f%z/=kept_f%z).or. &
      any(hf%x/=kept_hf%x).or.any(hf%y/=kept_hf%y).or.any(hf%z/=kept_hf%z)) &
      error stop 'late MOL held paired flux rollback'
  end subroutine check_held_input
  subroutine check_spatial_halo()
    type(eb_mol_diagnostics_3d) :: halo_diagnostic
    type(eb_geometry_3d) :: wide,center,native
    type(eb_integrated_flux_3d) :: grown_flux,valid_flux,native_flux,kept_flux,raw_flux, &
      grown_aux_flux,valid_aux_flux,kept_aux_flux
    type(rotating_frame_3d), allocatable :: spin
    real(dp), allocatable :: u(:,:,:,:),temperature(:,:,:),grown_rate(:,:,:,:),valid_rate(:,:,:,:), &
      native_rate(:,:,:,:),halo(:,:,:,:),kept_halo(:,:,:,:),alpha(:,:,:,:),phi_wide(:,:,:),pot, &
      aux_halo(:,:,:),kept_aux_halo(:,:,:)
    real(dp) :: grown_aux(10,10,10),valid_aux(4,4,4)
    real(dp) :: prim(reactive_nprim(size(species))),b(nv),w(nv),bn(nv),wn(nv),bg(nv),wg(nv),tt,cs,pos(3),balance
    integer :: ii,jj,kk,cc,rotation_id,solver_id,order_id
    logical :: good
    allocate(phi_wide(0:16,0:16,0:16),u(nv,16,16,16),temperature(16,16,16), &
      grown_rate(nv,10,10,10),valid_rate(nv,4,4,4),native_rate(nv,4,4,4))
    do kk=0,16
    do jj=0,16
    do ii=0,16
      phi_wide(ii,jj,kk)=(-1.5_dp+.25_dp*ii)-.3_dp*(-1.5_dp+.25_dp*jj)-.44_dp
    enddo
    enddo
    enddo
    call build_level_set_eb_geometry_3d(phi_wide,[-1.5_dp,-1.5_dp,-1.5_dp],[2.5_dp,2.5_dp,2.5_dp],wide,good)
    if(.not.good) error stop 'spatial halo geometry'
    call extract_eb_geometry_3d(wide,[7,7,7],[10,10,10],center,good)
    if(.not.good) error stop 'spatial halo valid geometry'
    call extract_eb_geometry_3d(wide,[4,4,4],[13,13,13],native,good)
    if(.not.good) error stop 'spatial halo native geometry'
    do rotation_id=0,1
      if(rotation_id==1) then
        allocate(spin,pot)
        spin%angular_speed=omega;spin%axis=3;spin%axis_point=[.13_dp,-.07_dp,.2_dp]
      endif
      do kk=1,16
      do jj=1,16
      do ii=1,16
        pos=[wide%x_lower+(ii-.5_dp)*wide%dx,wide%y_lower+(jj-.5_dp)*wide%dy,wide%z_lower+(kk-.5_dp)*wide%dz]
        prim=0;prim(1)=1+.01_dp*pos(1);prim(2:4)=.1_dp*pos;prim(5)=150000+100*pos(2)
        do cc=1,size(species)
          if(species(cc)%name=='N2') prim(cc+5)=1
        enddo
        if(allocated(spin)) then
          call rotating_potential_3d(spin,pos,pot,good)
          if(.not.good) error stop 'spatial halo potential'
        endif
        call reactive_primitive_to_conserved(species,prim,u(:,ii,jj,kk),tt,cs,good,pot)
        if(.not.good) error stop 'spatial halo state'
        temperature(ii,jj,kk)=tt
        if(wide%volume_fraction(ii,jj,kk)<=0) u(:,ii,jj,kk)=ieee_value(0._dp,ieee_quiet_nan)
      enddo
      enddo
      enddo
      do solver_id=0,1
      do order_id=1,2
        call compute_mol_eb_hydro_3d(species,wide,u,temperature,dt,order_id,grown_rate,grown_flux,bg,wg,good, &
          lax_friedrichs=solver_id==1,frame=spin,internal_energy_rate=grown_aux, &
          internal_energy_integrated_flux=grown_aux_flux)
        if(.not.good) error stop 'spatial halo grown control'
        call compute_mol_eb_hydro_3d(species,wide,u,temperature,dt,order_id,valid_rate,valid_flux,b,w,good, &
          lax_friedrichs=solver_id==1,frame=spin,spatial_rate_halo=halo, &
          internal_energy_rate=valid_aux,internal_energy_integrated_flux=valid_aux_flux, &
          internal_energy_spatial_rate_halo=aux_halo,diagnostics=halo_diagnostic)
        if(.not.good) error stop 'spatial halo evaluation'
        if(any(halo_diagnostic%rate/=grown_rate)) error stop 'diagnostic halo rate indexing'
        if(any(shape(halo_diagnostic%cell_pressure)/=[16,16,16])) error stop 'diagnostic input pressure shape'
        if(any(dt*halo_diagnostic%centroid_flux%x(1:nv,:,:,:)/=grown_flux%x)) error stop 'diagnostic halo flux indexing'
        if(any(shape(halo)/=[nv,16,16,16])) error stop 'spatial halo shape'
        if(any(halo(:,4:13,4:13,4:13)/=grown_rate)) error stop 'spatial halo equals grown operator'
        if(any(shape(aux_halo)/=[16,16,16]).or.any(aux_halo(4:13,4:13,4:13)/=grown_aux)) &
          error stop 'retained spatial halo equals grown operator'
        if(any(valid_aux/=grown_aux(4:7,4:7,4:7))) error stop 'retained spatial halo valid crop'
        if(any(valid_aux_flux%x/=grown_aux_flux%x(:,3:7,4:7,4:7)).or. &
          any(valid_aux_flux%y/=grown_aux_flux%y(:,4:7,3:7,4:7)).or. &
          any(valid_aux_flux%z/=grown_aux_flux%z(:,4:7,4:7,3:7))) &
          error stop 'retained spatial halo face crop'
        if(any(valid_rate/=grown_rate(:,4:7,4:7,4:7))) error stop 'spatial halo valid crop'
        if(any(valid_flux%x/=grown_flux%x(:,3:7,4:7,4:7)).or. &
          any(valid_flux%y/=grown_flux%y(:,4:7,3:7,4:7)).or. &
          any(valid_flux%z/=grown_flux%z(:,4:7,4:7,3:7))) error stop 'spatial halo face crop'
        if(any(halo(:,:3,:,:)/=0).or.any(halo(:,14:,:,:)/=0).or.any(halo(:,:,:3,:)/=0).or. &
          any(halo(:,:,14:,:)/=0).or.any(halo(:,:,:,:3)/=0).or.any(halo(:,:,:,14:)/=0)) &
          error stop 'spatial halo unused layers zero'
        balance=0
        do cc=1,nv
          balance=max(balance,abs(dt*center%dx*center%dy*center%dz*sum(center%volume_fraction*valid_rate(cc,:,:,:))-b(cc)-w(cc)))
        enddo
        if(balance>1.e-15_dp) error stop 'spatial halo valid conservation'
        if(solver_id==0) then
          call compute_mol_eb_hydro_3d(species,native,u(:,4:13,4:13,4:13),temperature(4:13,4:13,4:13), &
            dt,order_id,native_rate,native_flux,bn,wn,good,frame=spin)
          if(.not.good) error stop 'spatial halo independent native control'
          if(maxval(abs(native_rate-valid_rate)/max(1._dp,abs(native_rate)))>1.e-12_dp) &
            error stop 'spatial halo native valid equality'
        else
          call compute_mol_eb_hydro_3d(species,wide,u,temperature,dt,order_id,valid_rate,valid_flux,b,w,good, &
            lax_friedrichs=.true.,frame=spin,raw_face_flux=raw_flux,limiter_alpha=alpha, &
            include_boundary_correction=.true.,spatial_rate_halo=halo)
          if(.not.good.or.any(shape(alpha)/=[nv-5,4,4,4])) error stop 'spatial halo limiter output'
          if(any(shape(raw_flux%x)/=[nv,5,4,4]).or.lbound(raw_flux%x,2)/=0) error stop 'spatial halo raw flux bounds'
        endif
      enddo
      enddo
    enddo
    kept_flux=valid_flux;kept_halo=halo;kept_aux_halo=aux_halo;kept_aux_flux=valid_aux_flux
    valid_rate=17;valid_aux=29;b=19;w=23
    u(1,10,10,10)=ieee_value(0._dp,ieee_quiet_nan)
    call compute_mol_eb_hydro_3d(species,wide,u,temperature,dt,2,valid_rate,valid_flux,b,w,good,frame=spin, &
      spatial_rate_halo=halo,internal_energy_rate=valid_aux,internal_energy_integrated_flux=valid_aux_flux, &
      internal_energy_spatial_rate_halo=aux_halo)
    if(good.or.any(valid_rate/=17).or.any(b/=19).or.any(w/=23).or.any(halo/=kept_halo)) error stop 'spatial halo joint rollback'
    if(any(valid_flux%x/=kept_flux%x).or.any(valid_flux%y/=kept_flux%y).or.any(valid_flux%z/=kept_flux%z)) &
      error stop 'spatial halo flux rollback'
    if(any(valid_aux/=29).or.any(aux_halo/=kept_aux_halo).or. &
      any(valid_aux_flux%x/=kept_aux_flux%x).or.any(valid_aux_flux%y/=kept_aux_flux%y).or. &
      any(valid_aux_flux%z/=kept_aux_flux%z)) error stop 'retained spatial halo rollback'
  end subroutine check_spatial_halo
end program
