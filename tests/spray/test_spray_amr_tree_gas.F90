program test_spray_amr_tree_gas
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d,average_down_amr_patch_tree_3d
  use spray_amr_tree_coupling_mod, only: deposit_amr_tree_spray_gas,deposit_amr_tree_spray_field
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: gas,temp,saved,saved_t,expected,averaged
  real(dp), allocatable :: root(:,:,:,:),primitive(:),cell(:),y(:),feedback(:,:),before(:),after(:),scale(:)
  real(dp) :: positions(3,3),spacing(3)=.001_dp,zero(3)=0._dp,t,sound,r
  integer :: nv,node,n(3),i,j,k,p,water,nitrogen,mode
  logical :: ok,periodic(3)=[.true.,.true.,.true.]
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermodynamics')
  nv=reactive_nvar(size(species))
  allocate(root(nv,4,4,4),primitive(nv),cell(nv),y(size(species)),feedback(nv,3),before(nv),after(nv),scale(nv))
  water=0;nitrogen=0
  do p=1,size(species)
    if (species(p)%name=='H2O') water=p
    if (species(p)%name=='N2') nitrogen=p
  end do
  call require(water>0.and.nitrogen>0,'species identity')
  y=0._dp;y(nitrogen)=1._dp
  r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0._dp;primitive(1)=1._dp;primitive(5)=r*1000._dp;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  do k=1,4
  do j=1,4
  do i=1,4
    root(:,i,j,k)=cell
  end do
  end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  call prolong_amr_patch_tree_3d(tree,root,gas,ok)
  call require(ok,'gas fields')
  allocate(temp%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(temp%fields(node)%values(1,n(1),n(2),n(3)))
    temp%fields(node)%values=1000._dp
  end do
  y=0._dp;y(water)=1._dp;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'vapor gas constant')
  primitive(1)=1._dp;primitive(2:4)=[10._dp,-20._dp,30._dp];primitive(5)=r*1500._dp;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'vapor feedback state')
  positions(:,1)=spacing*[1.9_dp,1.1_dp,1.1_dp]
  positions(:,2)=spacing*[.9_dp,.9_dp,.9_dp]
  positions(:,3)=spacing*[3.9_dp,1.1_dp,1.1_dp]
  do p=1,3
    feedback(:,p)=real(p,dp)*1.e-11_dp*cell
  end do
  call composite_integrals_amr_patch_tree_3d(tree,gas,spacing,before,ok)
  call require(ok,'initial inventory')
  expected=gas
  call deposit_amr_tree_spray_field(tree,expected,positions,feedback,zero,spacing,periodic,ok)
  call require(ok,'field reference')
  call deposit_amr_tree_spray_gas(species,tree,gas,temp,positions,feedback,zero,spacing,periodic,ok)
  call require(ok,'thermodynamic feedback')
  averaged=gas
  call average_down_amr_patch_tree_3d(tree,averaged,ok)
  call require(ok,'parent consistency')
  do node=1,3
    call require(all(gas%fields(node)%values==expected%fields(node)%values),'conserved payload unmodified by EOS')
    call require(all(gas%fields(node)%values==averaged%fields(node)%values),'covered-parent consistency')
    n=tree%nodes(node)%shape
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      call reactive_conserved_to_primitive(species,gas%fields(node)%values(:,i,j,k),1000._dp,primitive,t,sound,ok)
      call require(ok.and.t==temp%fields(node)%values(1,i,j,k),'independent cell recovery')
    end do
    end do
    end do
  end do
  call require(maxval(temp%fields(3)%values)>1000._dp,'vapor heating')
  call composite_integrals_amr_patch_tree_3d(tree,gas,spacing,after,ok)
  scale=max(abs(before),sum(abs(feedback),dim=2),tiny(1._dp))
  call require(ok.and.all(abs(after-before-sum(feedback,dim=2))<=3.e-13_dp*scale),'mass momentum energy species budget')
  saved=gas;saved_t=temp
  do mode=1,3
    feedback=0._dp
    select case(mode)
    case(1)
      feedback(1,3)=-1._dp ! Negative receiving density.
    case(2)
      feedback(5,3)=-1._dp ! Energy below the thermodynamic range.
    case(3)
      feedback(5+nitrogen,3)=-1._dp ! Negative receiving species mass.
    end select
    call deposit_amr_tree_spray_gas(species,tree,gas,temp,positions,feedback,zero,spacing,periodic,ok)
    call require(.not.ok,'inadmissible feedback rejection')
    call unchanged()
  end do
  call deposit_amr_tree_spray_gas(species,tree,gas,temp,positions(:,:0),feedback(:,:0),zero,spacing,periodic,ok)
  call require(ok,'empty feedback');call unchanged()
  gas%fields(3)%values(1,1,1,1)=-1._dp;saved=gas
  feedback=0._dp
  call deposit_amr_tree_spray_gas(species,tree,gas,temp,positions,feedback,zero,spacing,periodic,ok)
  call require(.not.ok,'invalid initial gas rejected');call unchanged()
  print *, 'PASS: transactional AMR spray gas feedback'
contains
  subroutine unchanged()
    integer :: q
    do q=1,3
      call require(all(gas%fields(q)%values==saved%fields(q)%values),'conserved rollback')
      call require(all(temp%fields(q)%values==saved_t%fields(q)%values),'temperature rollback')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
