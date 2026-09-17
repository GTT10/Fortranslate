program test_eb_remap_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,remap_reactive_eb_multipatch_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_reactive_patch_3d) :: old(2),target(2),saved(2),empty(0),permuted(2)
  type(eb_geometry_3d) :: full,coarse
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),su(:,:,:,:),st(:,:,:),before(:),after(:)
  integer :: bad_lower(3,2),ol(3,2),nl(3,2),empty_lower(3,0),r,nv,p,i,j,k,axis,new_r
  logical :: ok
  character(len=1), parameter :: axes(3)=['x','y','z']
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(u(nv,4,4,4),t(4,4,4),before(nv),after(nv))
  do axis=1,3
  do r=2,3
    call build_axis_plane_eb_geometry_3d(4*r,4*r,4*r,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      axes(axis),.31_dp,full,ok)
    call require(ok,'fine geometry')
    call coarsen_eb_geometry_3d(full,r,coarse,ok)
    call require(ok,'coarse geometry')
    config%nx=4;config%ny=4;config%nz=4
    call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
    call require(ok,'initialize')
    ol(:,1)=[1,1,1];ol(:,2)=[3,3,3]
    nl(:,1)=[2,1,1];nl(:,2)=[3,3,2]
    do p=1,2
      call setup(old(p),ol(:,p),r)
      config%nx=2*r;config%ny=2*r;config%nz=2*r
      call initialize_reactive_eb_density_sheet_3d(species,config,old(p)%geometry,old(p)%state,old(p)%temperature,ok)
      call require(ok,'fine initialize')
      do k=1,2*r
      do j=1,2*r
      do i=1,2*r
        old(p)%state(:,i,j,k)=old(p)%state(:,i,j,k)*(1._dp+.01_dp*(i+3*j+5*k+p))
      end do
      end do
      end do
      call setup(target(p),nl(:,p),r)
    end do
    su=u;st=t
    call inventory(old,ol,r,before)
    call remap_reactive_eb_multipatch_3d(species,coarse,ol,r,old,nl,r,target,u,t,ok)
    call require(ok,'move patches')
    call inventory(target,nl,r,after)
    call require(maxval(abs(after-before)/max(1._dp,abs(before)))<2.e-14_dp,'move conservation')
    call require(all(target(1)%state(:,1:r,:,:)==old(1)%state(:,r+1:2*r,:,:)),'exact overlap state')
    call require(all(target(1)%temperature(1:r,:,:)==old(1)%temperature(r+1:2*r,:,:)),'exact overlap temperature')
    saved=target
    u=su;t=st
    permuted=target(2:1:-1)
    call remap_reactive_eb_multipatch_3d(species,coarse,ol(:,2:1:-1),r,old(2:1:-1), &
      nl(:,2:1:-1),r,permuted,u,t,ok)
    call require(ok,'permuted remap')
    call require(all(permuted(2)%state==saved(1)%state).and.all(permuted(1)%state==saved(2)%state), &
      'permutation exactness')
    su=u;st=t
    old(1)%state(1,2*r,2*r,2*r)=-1._dp
    call remap_reactive_eb_multipatch_3d(species,coarse,ol,r,old,nl,r,target,u,t,ok)
    call require(.not.ok,'invalid source rejection')
    call require(all(u==su).and.all(t==st).and.all(target(1)%state==saved(1)%state).and. &
      all(target(2)%state==saved(2)%state),'failed remap rollback')
    bad_lower=nl;bad_lower(:,2)=nl(:,1)
    call remap_reactive_eb_multipatch_3d(species,coarse,nl,r,saved,bad_lower,r,target,u,t,ok)
    call require(.not.ok,'overlapping destination rejection')
    call require(all(u==su).and.all(t==st).and.all(target(1)%state==saved(1)%state).and. &
      all(target(2)%state==saved(2)%state),'topology rejection rollback')
    u(1,4,2,4)=-1._dp
    call remap_reactive_eb_multipatch_3d(species,coarse,nl,r,saved,nl,r,target,u,t,ok)
    call require(.not.ok,'invalid exposed coarse state rejection')
    call require(u(1,4,2,4)==-1._dp.and.all(t==st).and.all(target(1)%state==saved(1)%state), &
      'invalid coarse state rollback')
    u=su
    ! Remove every fine patch, then introduce refinement with a changed ratio.
    call remap_reactive_eb_multipatch_3d(species,coarse,nl,r,saved,empty_lower,r,empty,u,t,ok)
    call require(ok,'complete coarsening')
    call inventory(empty,empty_lower,r,after)
    call require(maxval(abs(after-before)/max(1._dp,abs(before)))<2.e-14_dp,'coarsening conservation')
    new_r=5-r
    call build_axis_plane_eb_geometry_3d(4*new_r,4*new_r,4*new_r,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      axes(axis),.31_dp,full,ok)
    call require(ok,'changed refinement geometry')
    do p=1,2
      call setup(target(p),ol(:,p),new_r)
    end do
    call remap_reactive_eb_multipatch_3d(species,coarse,empty_lower,r,empty,ol,new_r,target,u,t,ok)
    call require(ok,'refinement from coarse only')
    call inventory(target,ol,new_r,after)
    call require(maxval(abs(after-before)/max(1._dp,abs(before)))<2.e-14_dp,'refinement conservation')
    ! Direct ratio change with existing fine patches also conserves their
    ! subcell-varying inventories, using the restricted coarse donor.
    call remap_reactive_eb_multipatch_3d(species,coarse,nl,r,saved,ol,new_r,target,u,t,ok)
    call require(ok,'direct refinement ratio change')
    call inventory(target,ol,new_r,after)
    call require(maxval(abs(after-before)/max(1._dp,abs(before)))<2.e-14_dp,'ratio-change conservation')
  end do
  end do
  print *, 'test_eb_remap_3d: PASS'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
  subroutine setup(patch,lower,ratio)
    type(eb_reactive_patch_3d), intent(inout) :: patch
    integer, intent(in) :: lower(3),ratio
    call extract_eb_geometry_3d(full,(lower-1)*ratio+1,(lower+1)*ratio,patch%geometry,ok)
    call require(ok,'extract patch')
    if (allocated(patch%state)) deallocate(patch%state,patch%temperature)
    allocate(patch%state(nv,2*ratio,2*ratio,2*ratio),patch%temperature(2*ratio,2*ratio,2*ratio))
    patch%state=-17._dp;patch%temperature=-19._dp
  end subroutine setup
  subroutine inventory(patches,lower,ratio,amount)
    type(eb_reactive_patch_3d), intent(in) :: patches(:)
    integer, intent(in) :: lower(:,:),ratio
    real(dp), intent(out) :: amount(:)
    logical :: refined(4,4,4)
    integer :: a,b,c,q,d(3)
    real(dp) :: volume
    refined=.false.;amount=0._dp
    do q=1,size(patches)
      d=[patches(q)%geometry%nx,patches(q)%geometry%ny,patches(q)%geometry%nz]
      refined(lower(1,q):lower(1,q)+d(1)/ratio-1,lower(2,q):lower(2,q)+d(2)/ratio-1, &
        lower(3,q):lower(3,q)+d(3)/ratio-1)=.true.
      volume=patches(q)%geometry%dx*patches(q)%geometry%dy*patches(q)%geometry%dz
      do c=1,d(3)
      do b=1,d(2)
      do a=1,d(1)
        amount=amount+volume*patches(q)%geometry%volume_fraction(a,b,c)*patches(q)%state(:,a,b,c)
      end do
      end do
      end do
    end do
    volume=coarse%dx*coarse%dy*coarse%dz
    do c=1,4
    do b=1,4
    do a=1,4
      if (.not.refined(a,b,c)) amount=amount+volume*coarse%volume_fraction(a,b,c)*u(:,a,b,c)
    end do
    end do
    end do
  end subroutine inventory
end program test_eb_remap_3d
