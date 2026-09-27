program test_eb_mol_hydro_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: g,physical
  type(eb_integrated_flux_3d) :: flux,saved
  type(rotating_frame_3d), allocatable :: frame
  real(dp), allocatable :: potential,q_back(:),ordinary(:),zero_energy(:)
  real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),q(:),rates(:,:,:,:),boundary(:),wall(:),inventory(:)
  real(dp) :: phi(0:10,0:10,0:10),sound,t,residual,dt,position(3),t_back,c_back,omega
  integer :: i,j,k,c,nv,mode,order,rotation,solver,cases
  logical :: ok
  character(len=64) :: argument
  omega=400._dp
  if(command_argument_count()>0) then
    call get_command_argument(1,argument)
    read(argument,*) omega
  end if
  call load_h2o2_elementary_thermo(species,ok)
  if(.not.ok) error stop 'thermo'
  nv=reactive_nvar(size(species));dt=1e-7_dp
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
    do solver=0,1
    do order=1,2
      call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok,lax_friedrichs=solver==1,frame=frame)
      if(.not.ok) error stop 'MOL spatial rate'
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
  solver=0
  saved=flux;rates=17;boundary=19;wall=23
  state(1,10,10,10)=ieee_value(0._dp,ieee_quiet_nan)
  call compute_mol_eb_hydro_3d(species,g,state,temp,dt,2,rates,flux,boundary,wall,ok,lax_friedrichs=solver==1,frame=frame)
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
end program
