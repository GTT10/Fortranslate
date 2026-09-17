program test_eb_coarse_only_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,advance_reactive_eb_multipatch_full_3d
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,reactive_eb_integrals_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(eb_geometry_3d) :: geometry
  type(eb_reactive_patch_3d) :: empty(0)
  type(eb_integrated_flux_3d) :: flux
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: initial(:,:,:,:),initial_t(:,:,:),u(:,:,:,:),t(:,:,:),reference(:,:,:,:),reference_t(:,:,:)
  real(dp), allocatable :: stage(:,:,:,:),stage_t(:,:,:),b(:),w(:),chem(:),before(:),after(:)
  real(dp) :: theta,half_theta,expected_theta,dt
  integer :: lower(3,0),nv,mode,integrator_index
  logical :: ok,chemistry,transport_on,subcycle
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nv=reactive_nvar(size(species));dt=1.e-8_dp
  call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,geometry,ok)
  call require(ok,'geometry')
  allocate(initial(nv,4,4,4),initial_t(4,4,4),b(nv),w(nv),chem(nv),before(nv),after(nv))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,geometry,initial,initial_t,ok)
  call require(ok,'initial state')
  call reactive_eb_integrals_3d(initial,geometry,before,ok)
  call require(ok,'initial inventory')
  do integrator_index=1,2
  do mode=0,7
    chemistry=btest(mode,0);transport_on=btest(mode,1);subcycle=btest(mode,2)
    u=initial;t=initial_t;reference=initial;reference_t=initial_t;expected_theta=1._dp
    if (chemistry) call reaction()
    if (transport_on) call diffuse()
    stage=reference;stage_t=reference_t
    call advance_reactive_eb_euler_3d(species,stage,stage_t,geometry,'rusanov',dt,reference,reference_t,ok)
    call require(ok,'independent hydro')
    if (transport_on) call diffuse()
    if (chemistry) call reaction()
    call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,geometry,lower,3,'rusanov',dt, &
      chemistry,2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,u,t,empty,b,w,chem,theta,ok, &
      integrators(integrator_index),subcycle)
    call require(ok,'zero-patch composition')
    call require(all(u==reference).and.all(t==reference_t),'exact independent split parity')
    call require(theta==expected_theta,'transport limiter parity')
    call reactive_eb_integrals_3d(u,geometry,after,ok)
    call require(ok,'final inventory')
    call require(all(abs(after-before-b-w-chem)<=2.e-14_dp*max(1._dp,abs(before))),'component budget')
    if (chemistry) call require(any(abs(chem(6:))>1.e-15_dp),'active chemical species budget')
    u=initial;t=initial_t;b=17._dp;w=19._dp;chem=23._dp;theta=.75_dp
    call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,geometry,lower,3,'invalid',dt, &
      chemistry,2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,u,t,empty,b,w,chem,theta,ok, &
      integrators(integrator_index),subcycle)
    call require(.not.ok,'rejected hydro after optional reaction and transport')
    call require(all(u==initial).and.all(t==initial_t).and.all(b==17._dp).and.all(w==19._dp).and. &
      all(chem==23._dp).and.theta==.75_dp,'state and diagnostic rollback')
  end do
  end do
  print *, 'test_eb_coarse_only_3d: PASS'
contains
  subroutine reaction()
    call advance_reactive_chemistry_3d(species,reactions,reference,reference_t,4,4,4,.5_dp*dt, &
      2.e-7_dp,1.e-12_dp,ok,geometry%volume_fraction>0._dp,integrators(integrator_index))
    call require(ok,'independent reaction half')
  end subroutine reaction
  subroutine diffuse()
    call advance_reactive_eb_transport_patch_euler_3d(species,transport,geometry,.5_dp*dt,.true.,.true., &
      .true.,.true.,reference,reference_t,flux,half_theta,ok)
    call require(ok,'independent transport half')
    expected_theta=min(expected_theta,half_theta)
  end subroutine diffuse
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_coarse_only_3d
