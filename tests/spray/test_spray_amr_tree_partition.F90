program test_spray_amr_tree_partition
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,initialize_amr_patch_tree_3d
  use spray_parcel_mod, only: spray_parcel
  use spray_amr_tree_partition_mod
  implicit none
  type(amr_patch_3d) :: patches(3)
  type(amr_patch_tree_3d) :: tree,coarse
  type(spray_parcel), allocatable :: parcels(:),original(:),saved(:)
  real(dp), allocatable :: mass(:,:),original_mass(:,:),saved_mass(:,:)
  real(dp) :: positions(3,6),zero(3)=0._dp,spacing(3)=1._dp
  integer, allocatable :: locations(:,:),counts(:),offsets(:),previous(:,:)
  integer :: p,q,owners(4)=[2,1,0,3]
  logical :: ok
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'left patch')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,patches(2),ok)
  call require(ok,'right patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,4,1,4,2,patches(3),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patches,tree,ok)
  call require(ok,'tree')
  positions(:,1)=[.25_dp,.25_dp,.25_dp]
  positions(:,2)=[1._dp,.25_dp,.25_dp]
  positions(:,3)=[2._dp,.25_dp,.25_dp]
  positions(:,4)=[4.25_dp,.25_dp,.25_dp]
  positions(:,5)=[-.25_dp,.25_dp,.25_dp]
  positions(:,6)=[.25_dp,2._dp,.25_dp]
  call locate_amr_tree_spray_cells(tree,positions,zero,spacing,[.true.,.false.,.false.],locations,ok)
  call require(ok,'periodic and face lookup')
  call require(all(locations(1,:)==[4,2,3,4,3,2]),'finest half-open node')
  call require(all(locations(:,1)==[4,2,2,2]),'nested local cell')
  call require(all(locations(:,2)==[2,3,1,1]),'nested upper face')
  call require(all(locations(:,3)==[3,1,1,1]),'peer face')
  call require(all(locations(:,5)==[3,4,1,1]),'negative periodic wrap')
  previous=locations
  call locate_amr_tree_spray_cells(tree,positions,zero,spacing,[.false.,.false.,.false.],locations,ok)
  call require(.not.ok.and.all(locations==previous),'nonperiodic out-of-domain rollback')
  positions(:,4)=[4._dp,.25_dp,.25_dp]
  call locate_amr_tree_spray_cells(tree,positions(:,:4),zero,spacing,[.false.,.false.,.false.],locations,ok)
  call require(.not.ok.and.all(locations==previous),'upper domain face rejection')
  positions(:,4)=[4.25_dp,.25_dp,.25_dp]
  allocate(parcels(6),mass(2,6))
  do p=1,6
    parcels(p)=spray_parcel(int(p,int64),positions(:,p),[1._dp,2._dp,3._dp], &
      1.e-9_dp*p,350._dp,100._dp)
    parcels(p)%parent_id=p-1;parcels(p)%material_id=1+mod(p,2)
    parcels(p)%deformation=.1_dp*p;parcels(p)%deformation_rate=.2_dp*p
    parcels(p)%rt_time=.01_dp*p;parcels(p)%khrt_initial_count=100._dp
    parcels(p)%khrt_shed_volume=1.e-15_dp*p
    mass(:,p)=parcels(p)%mass*[.25_dp,.75_dp]
  end do
  original=parcels;original_mass=mass
  call partition_amr_tree_spray_parcels(tree,owners,5,parcels,zero,spacing,[.true.,.false.,.false.], &
    counts,offsets,locations,ok,mass)
  call require(ok,'stable owner partition')
  call require(all(counts==[2,2,0,2,0]).and.all(offsets==[0,2,4,4,6]),'empty owners and offsets')
  call require(all(parcels%id==[3,5,2,6,1,4]),'stable rank order')
  do p=1,6
    q=int(parcels(p)%id)
    call require(all(mass(:,p)==original_mass(:,q)),'component association')
    call require(all(parcels(p)%position==original(q)%position).and. &
      all(parcels(p)%velocity==original(q)%velocity),'kinematics unchanged')
    call require(parcels(p)%parent_id==original(q)%parent_id.and. &
      parcels(p)%material_id==original(q)%material_id,'identity unchanged')
    call require(all([parcels(p)%mass,parcels(p)%temperature,parcels(p)%multiplicity, &
      parcels(p)%deformation,parcels(p)%deformation_rate,parcels(p)%rt_time, &
      parcels(p)%khrt_initial_count,parcels(p)%khrt_shed_volume]== &
      [original(q)%mass,original(q)%temperature,original(q)%multiplicity, &
      original(q)%deformation,original(q)%deformation_rate,original(q)%rt_time, &
      original(q)%khrt_initial_count,original(q)%khrt_shed_volume]),'thermodynamics and breakup unchanged')
  end do
  parcels=original;mass=original_mass;mass(1,2)=-1._dp
  saved=parcels;saved_mass=mass;previous=locations
  call partition_amr_tree_spray_parcels(tree,owners,5,parcels,zero,spacing,[.true.,.false.,.false.], &
    counts,offsets,locations,ok,mass)
  call require(.not.ok.and.all(parcels%id==saved%id).and.all(mass==saved_mass).and. &
    all(locations==previous).and.all(counts==[2,2,0,2,0]),'invalid inventory rollback')
  mass=original_mass;parcels(2)%id=parcels(1)%id
  call partition_amr_tree_spray_parcels(tree,owners,5,parcels,zero,spacing,[.true.,.false.,.false.], &
    counts,offsets,locations,ok,mass)
  call require(.not.ok,'duplicate ID rejection')
  parcels=original;owners(4)=5
  call partition_amr_tree_spray_parcels(tree,owners,5,parcels,zero,spacing,[.true.,.false.,.false.], &
    counts,offsets,locations,ok,mass)
  call require(.not.ok.and.all(locations==previous),'invalid owner rollback')
  positions(1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call locate_amr_tree_spray_cells(tree,positions,zero,spacing,[.true.,.false.,.false.],locations,ok)
  call require(.not.ok.and.all(locations==previous),'nonfinite lookup rollback')
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),coarse,ok)
  call require(ok,'coarse regrid')
  call partition_amr_tree_spray_parcels(coarse,[1],3,parcels,zero,spacing,[.true.,.false.,.false.], &
    counts,offsets,locations,ok,mass)
  call require(ok.and.all(locations(1,:)==1).and.all(counts==[0,6,0]),'coarsened topology reassignment')
  deallocate(parcels,mass);allocate(parcels(0),mass(2,0))
  call partition_amr_tree_spray_parcels(coarse,[0],2,parcels,zero,spacing,[.false.,.false.,.false.], &
    counts,offsets,locations,ok,mass)
  call require(ok.and.all(counts==0).and.size(locations,2)==0,'empty population')
  print *, 'PASS: AMR spray ownership and partition'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
