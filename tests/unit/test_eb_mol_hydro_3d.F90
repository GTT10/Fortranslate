program test_eb_mol_hydro_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: g,physical
  type(eb_integrated_flux_3d) :: flux,saved
  real(dp), allocatable :: state(:,:,:,:),temp(:,:,:),q(:),rates(:,:,:,:),boundary(:),wall(:),inventory(:)
  real(dp) :: phi(0:10,0:10,0:10),sound,t,residual,dt
  integer :: i,j,k,c,nv,mode,order
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  if(.not.ok) error stop 'thermo'
  nv=reactive_nvar(size(species));dt=1e-7_dp
  allocate(state(nv,10,10,10),temp(10,10,10),q(reactive_nprim(size(species))), &
    rates(nv,4,4,4),boundary(nv),wall(nv),inventory(nv))
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
      call reactive_primitive_to_conserved(species,q,state(:,i,j,k),t,sound,ok)
      if(.not.ok) error stop 'initial state'
      temp(i,j,k)=t
      if(g%volume_fraction(i,j,k)<=0) state(:,i,j,k)=ieee_value(0._dp,ieee_quiet_nan)
    end do
    end do
    end do
    do order=1,2
      call compute_mol_eb_hydro_3d(species,g,state,temp,dt,order,rates,flux,boundary,wall,ok)
      if(.not.ok) error stop 'MOL spatial rate'
      do c=1,nv
        inventory(c)=dt*physical%dx*physical%dy*physical%dz*sum(physical%volume_fraction*rates(c,:,:,:))
      end do
      residual=maxval(abs(inventory-boundary-wall))
      if(residual>1e-15_dp) error stop 'boundary and wall conservation'
      if(mode==1.and.maxval(abs(rates))>1e-7_dp) error stop 'stationary inclined wall balance'
      if(mode==2.and.maxval(abs(rates))<1) error stop 'nontrivial spatial rate'
      print *, 'mode, order, conservation residual:',mode,order,residual
    end do
  end do
  saved=flux;rates=17;boundary=19;wall=23
  state(1,10,10,10)=ieee_value(0._dp,ieee_quiet_nan)
  call compute_mol_eb_hydro_3d(species,g,state,temp,dt,2,rates,flux,boundary,wall,ok)
  if(ok.or.any(rates/=17).or.any(boundary/=19).or.any(wall/=23)) error stop 'rate rollback'
  if(any(flux%x/=saved%x).or.any(flux%y/=saved%y).or.any(flux%z/=saved%z)) error stop 'flux rollback'
  print *, 'PASS: MOL inclined-wall spatial rate and rollback'
end program
