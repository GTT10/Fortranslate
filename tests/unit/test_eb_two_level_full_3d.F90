program test_eb_two_level_full_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_species_component
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: advance_reactive_eb_two_level_euler_3d, &
    advance_reactive_eb_two_level_strang_3d,advance_reactive_eb_two_level_full_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  real(dp) :: theta
  real(dp), allocatable :: reference_chem(:)
  type(eb_geometry_3d) :: full,coarse,fine
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),uf(:,:,:,:),tf(:,:,:),initial_c(:,:,:,:),initial_tc(:,:,:)
  real(dp), allocatable :: initial_f(:,:,:,:),initial_tf(:,:,:),hydro_c(:,:,:,:),hydro_tc(:,:,:)
  real(dp), allocatable :: hydro_f(:,:,:,:),hydro_tf(:,:,:),b(:),w(:),chem(:),hb(:),hw(:),before(:),after(:)
  real(dp) :: atoms(3),element_change(3),element_scale(3)
  integer :: nvar,c,mode,s,component
  logical :: ok,refined(4,4,4)
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nvar=reactive_nvar(size(species))
  allocate(uc(nvar,4,4,4),uf(nvar,4,4,4),tc(4,4,4),tf(4,4,4))
  allocate(reference_chem(nvar))
  allocate(b(nvar),w(nvar),chem(nvar),hb(nvar),hw(nvar),before(nvar),after(nvar))
  call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,full,ok)
  call require(ok,'fine geometry')
  call coarsen_eb_geometry_3d(full,2,coarse,ok)
  call require(ok,'coarse geometry')
  call extract_eb_geometry_3d(full,[3,3,3],[6,6,6],fine,ok)
  call require(ok,'fine patch')
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,ok)
  call require(ok,'coarse initial')
  call initialize_reactive_eb_density_sheet_3d(species,config,fine,uf,tf,ok)
  call require(ok,'fine initial')
  initial_c=uc;initial_tc=tc;initial_f=uf;initial_tf=tf
  refined=.false.;refined(2:3,2:3,2:3)=.true.
  call inventory(before)
  call advance_reactive_eb_two_level_euler_3d(species,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
    uc,tc,uf,tf,hb,hw,ok)
  call require(ok,'hydro reference')
  hydro_c=uc;hydro_tc=tc;hydro_f=uf;hydro_tf=tf
  call reset()
  call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
    .false.,0._dp,0._dp,.false.,.true.,.true.,.true.,.true.,uc,tc,uf,tf,b,w,chem,theta,ok)
  call require(ok.and.all(uc==hydro_c).and.all(tc==hydro_tc).and.all(uf==hydro_f).and.all(tf==hydro_tf).and. &
    all(b==hb).and.all(w==hw).and.all(chem==0._dp),'disabled chemistry exact hydro parity')
  do mode=1,2
    call reset()
    ! Refined parent storage does not own chemical solves and is overwritten by
    ! validated restriction. An erroneous coarse solve here would fail on rho.
    uc(:,2:3,2:3,2:3)=-99._dp;tc(2:3,2:3,2:3)=-99._dp
    call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
      .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,uc,tc,uf,tf,b,w,chem,theta,ok,integrators(mode))
    call require(ok,'full reacting transport hierarchy')
    call require(theta>=0._dp.and.theta<=1._dp,'full composition limiter')
    call inventory(after)
    do c=1,nvar
      call require(abs(after(c)-before(c)-b(c)-w(c)-chem(c))<=3.e-12_dp*max(1._dp,abs(before(c))), &
        'composite boundary wall chemistry ledger')
    end do
    call require(maxval(abs(uf-hydro_f))>1.e-10_dp,'chemistry response nonzero')
    call require(maxval(abs(chem(1:5)))<=3.e-12_dp*max(1._dp,maxval(abs(before(1:5)))), &
      'chemistry preserves mass momentum energy')
    element_change=0._dp;element_scale=0._dp
    do s=1,size(species)
      select case(trim(species(s)%name))
      case('H2');atoms=[2._dp,0._dp,0._dp]
      case('H');atoms=[1._dp,0._dp,0._dp]
      case('O');atoms=[0._dp,1._dp,0._dp]
      case('O2');atoms=[0._dp,2._dp,0._dp]
      case('OH');atoms=[1._dp,1._dp,0._dp]
      case('H2O');atoms=[2._dp,1._dp,0._dp]
      case('N2');atoms=[0._dp,0._dp,2._dp]
      case default;error stop 'unexpected fixed species'
      end select
      component=reactive_species_component(s)
      element_change=element_change+atoms*chem(component)/species(s)%molecular_weight
      element_scale=element_scale+atoms*abs(before(component))/species(s)%molecular_weight
    end do
    call require(all(abs(element_change)<=2.e-10_dp*max(1.e-30_dp,element_scale)),'reaction element conservation')
    call require(all(tc>0._dp).and.all(tf>0._dp),'positive reacting temperatures')
    call require(all(uc(:,:,:,1)==initial_c(:,:,:,1)).and.all(tc(:,:,1)==initial_tc(:,:,1)), &
      'covered storage preserved')
    hydro_c=initial_c;hydro_tc=initial_tc;hydro_f=initial_f;hydro_tf=initial_tf
    call advance_reactive_eb_two_level_strang_3d(species,reactions,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
      .true.,2.e-7_dp,1.e-12_dp,hydro_c,hydro_tc,hydro_f,hydro_tf,hb,hw,reference_chem,ok,integrators(mode))
    call require(ok,'reaction hydro reference')
    call require(maxval(abs(uf-hydro_f))>1.e-12_dp,'full composition transport response')
    call reset()
    call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
      .true.,2.e-7_dp,1.e-12_dp,.false.,.true.,.true.,.true.,.true.,uc,tc,uf,tf,b,w,chem,theta,ok,integrators(mode))
    call require(ok.and.all(uc==hydro_c).and.all(tc==hydro_tc).and.all(uf==hydro_f).and.all(tf==hydro_tf).and. &
      all(b==hb).and.all(w==hw).and.all(chem==reference_chem).and.theta==1._dp, &
      'disabled transport exact reaction hydro parity')

  end do
  call reset()
  call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
    .false.,0._dp,0._dp,.true.,.true.,.true.,.true.,.true.,uc,tc,uf,tf,b,w,chem,theta,ok)
  call require(ok.and.all(chem==0._dp),'transport hydro with chemistry disabled')
  call inventory(after)
  do c=1,nvar
    call require(abs(after(c)-before(c)-b(c)-w(c))<=3.e-12_dp*max(1._dp,abs(before(c))), &
      'nonreacting transport hydro inventory')
  end do
  call reset()
  call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,[2,2,2],2,'rusanov',1.e-8_dp, &
    .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,uc,tc,uf,tf,b,w,chem,theta,ok,'invalid')
  call require(.not.ok,'bad chemistry integrator rejected')
  call unchanged()
  call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,[2,2,2],2,'invalid',1.e-8_dp, &
    .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true.,uc,tc,uf,tf,b,w,chem,theta,ok,'explicit')
  call require(.not.ok,'hydro failure after first reaction rejected')
  call unchanged()
  write(*,'(a)') 'test_eb_two_level_full_3d: PASS'
contains
  subroutine reset()
    uc=initial_c;tc=initial_tc;uf=initial_f;tf=initial_tf;b=17._dp;w=19._dp;chem=23._dp;theta=29._dp
  end subroutine reset
  subroutine unchanged()
    call require(all(uc==initial_c).and.all(tc==initial_tc).and.all(uf==initial_f).and.all(tf==initial_tf).and. &
      all(b==17._dp).and.all(w==19._dp).and.all(chem==23._dp).and.theta==29._dp,'reaction hierarchy transaction rollback')
  end subroutine unchanged
  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: component
    do component=1,size(values)
      values(component)=coarse%dx*coarse%dy*coarse%dz* &
        sum(uc(component,:,:,:)*coarse%volume_fraction,mask=.not.refined)+ &
        fine%dx*fine%dy*fine%dz*sum(uf(component,:,:,:)*fine%volume_fraction)
    end do
  end subroutine inventory
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) error stop message
  end subroutine require
end program test_eb_two_level_full_3d
