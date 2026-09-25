program test_eb_multipatch_full_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_species_component
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,advance_reactive_eb_multipatch_full_3d, &
    advance_reactive_eb_multipatch_euler_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(eb_geometry_3d) :: full,coarse
  type(eb_reactive_patch_3d) :: patches(2),initial(2),reference(2),swapped(2)
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),ic(:,:,:,:),it(:,:,:),rc(:,:,:,:),rt(:,:,:)
  real(dp), allocatable :: b(:),w(:),chem(:),rb(:),rw(:),rchem(:),before(:),after(:)
  real(dp) :: theta,rtheta,atoms(3),elements(3),scale(3)
  integer :: lower(3,2),upper(3,2),p,nvar,mode,s,component
  integer :: time_mode
  character(len=32) :: argument
  real(dp) :: coefficient=0
  logical :: subcycling
  logical :: valid,active(4,4,4)
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call get_command_argument(1,argument)
  if(argument=="viscosity") coefficient=.037_dp
  call load_h2o2_elementary_thermo(species,valid)
  call require(valid,'thermo')
  call load_h2o2_elementary_transport(transport,valid)
  call require(valid,'transport')
  call load_h2o2_elementary_mechanism(reactions,valid)
  call require(valid,'chemistry')
  nvar=reactive_nvar(size(species))
  allocate(uc(nvar,4,4,4),tc(4,4,4),b(nvar),w(nvar),chem(nvar),rb(nvar),rw(nvar),rchem(nvar),before(nvar),after(nvar))
  call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,full,valid)
  call require(valid,'full geometry')
  call coarsen_eb_geometry_3d(full,2,coarse,valid)
  call require(valid,'coarse geometry')
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,valid)
  call require(valid,'coarse initialization')
  lower(:,1)=[2,2,2];lower(:,2)=[3,2,2];upper(:,1)=[2,3,3];upper(:,2)=[3,3,3]
  active=coarse%volume_fraction>0._dp;active(2:3,2:3,2:3)=.false.
  do p=1,2
    call extract_eb_geometry_3d(full,(lower(:,p)-1)*2+1,upper(:,p)*2,patches(p)%geometry,valid)
    call require(valid,'fine extraction')
    allocate(patches(p)%state(nvar,2,4,4),patches(p)%temperature(2,4,4))
    config%nx=2
    call initialize_reactive_eb_density_sheet_3d(species,config,patches(p)%geometry, &
      patches(p)%state,patches(p)%temperature,valid)
    call require(valid,'fine initialization')
  end do
  ! These parent cells must never own reaction solves; restriction replaces them.
  uc(:,2:3,2:3,2:3)=-99._dp;tc(2:3,2:3,2:3)=-99._dp
  ic=uc;it=tc;initial=patches
  call inventory(before)
  do time_mode=1,2
  subcycling=time_mode==2
  do mode=1,2
    call reset()
    call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,coarse,lower,2,'rusanov',1.e-8_dp, &
      .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,uc,tc,patches,b,w,chem,theta,valid, &
        integrators(mode),subcycling,difmag=coefficient)
    call require(valid,'reacting transport multiple-patch composition')
    call inventory(after)
    call require(all(abs(after-before-b-w-chem)<=3.e-12_dp*max(1._dp,abs(before))),'full composite inventories')
    call require(all(abs(chem(1:5))<=3.e-12_dp*max(1._dp,abs(before(1:5)))),'reaction mass momentum energy')
    call require(theta>=0._dp.and.theta<=1._dp,'composition limiter')
    call require(maxval(abs(chem(6:)))>1.e-12_dp,'reaction species response')
    elements=0._dp;scale=0._dp
    do s=1,size(species)
      select case(trim(species(s)%name))
      case('H2');atoms=[2._dp,0._dp,0._dp]
      case('H');atoms=[1._dp,0._dp,0._dp]
      case('O');atoms=[0._dp,1._dp,0._dp]
      case('O2');atoms=[0._dp,2._dp,0._dp]
      case('OH');atoms=[1._dp,1._dp,0._dp]
      case('H2O');atoms=[2._dp,1._dp,0._dp]
      case('N2');atoms=[0._dp,0._dp,2._dp]
      case default;error stop 'unexpected species'
      end select
      component=reactive_species_component(s)
      elements=elements+atoms*chem(component)/species(s)%molecular_weight
      scale=scale+atoms*abs(before(component))/species(s)%molecular_weight
    end do
    call require(all(abs(elements)<=2.e-10_dp*max(1.e-30_dp,scale)),'reaction element inventories')
    rc=uc;rt=tc;reference=patches;rb=b;rw=w;rchem=chem;rtheta=theta
    call reset()
    swapped(1)=patches(2);swapped(2)=patches(1)
    call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,coarse,lower(:,2:1:-1),2,'rusanov',1.e-8_dp, &
      .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,uc,tc,swapped,b,w,chem,theta,valid, &
        integrators(mode),subcycling,difmag=coefficient)
    call require(valid.and.all(uc==rc).and.all(tc==rt).and.all(b==rb).and.all(w==rw).and. &
      all(chem==rchem).and.theta==rtheta,'full composition permutation invariance')
    do p=1,2
      call require(all(swapped(3-p)%state==reference(p)%state).and. &
        all(swapped(3-p)%temperature==reference(p)%temperature),'fine reacting permutation invariance')
    end do
  end do
  call reset()
  call advance_reactive_eb_multipatch_euler_3d(species,coarse,lower,2,'rusanov',1.e-8_dp,uc,tc,patches,rb,rw,valid, &
        subcycling,difmag=coefficient)
  call require(valid,'hydro reference')
  rc=uc;rt=tc;reference=patches
  call reset()
  call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,coarse,lower,2,'rusanov',1.e-8_dp, &
    .false.,0._dp,0._dp,.false.,.true.,.true.,.true.,.true.,uc,tc,patches,b,w,chem,theta,valid,subcycle=subcycling, &
        difmag=coefficient)
  call require(valid.and.all(uc==rc).and.all(tc==rt).and.all(b==rb).and.all(w==rw).and. &
    all(chem==0._dp).and.theta==1._dp,'disabled operators exact hydro parity')
  do p=1,2
    call require(all(patches(p)%state==reference(p)%state),'disabled operators fine parity')
  end do
  call reset()
  call advance_reactive_eb_multipatch_full_3d(species,reactions,transport,coarse,lower,2,'invalid',1.e-8_dp, &
    .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,uc,tc,patches,b,w,chem,theta,valid,'explicit', &
        subcycling,difmag=coefficient)
  call require(.not.valid,'hydro failure after accepted reaction and transport halves')
  call require(all(uc==ic).and.all(tc==it).and.all(b==17._dp).and.all(w==19._dp).and. &
    all(chem==23._dp).and.theta==29._dp,'all coarse and diagnostic outputs roll back')
  do p=1,2
    call require(all(patches(p)%state==initial(p)%state).and.all(patches(p)%temperature==initial(p)%temperature), &
      'all fine patches roll back')
  end do
  end do
  write(*,'(a)') 'test_eb_multipatch_full_3d: PASS'
contains
  subroutine reset()
    uc=ic;tc=it;patches=initial;b=17._dp;w=19._dp;chem=23._dp;theta=29._dp
  end subroutine reset
  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: c,q
    do c=1,nvar
      values(c)=coarse%dx*coarse%dy*coarse%dz*sum(uc(c,:,:,:)*coarse%volume_fraction,mask=active)
      do q=1,2
        values(c)=values(c)+patches(q)%geometry%dx*patches(q)%geometry%dy*patches(q)%geometry%dz* &
          sum(patches(q)%state(c,:,:,:)*patches(q)%geometry%volume_fraction)
      end do
    end do
  end subroutine inventory
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_multipatch_full_3d
