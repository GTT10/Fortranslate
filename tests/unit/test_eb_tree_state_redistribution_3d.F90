program test_eb_tree_state_redistribution_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,redistribute_eb_tree_state_3d, &
    redistribute_sdc_eb_tree_state_level_rate_3d,redistribute_sdc_eb_tree_state_rate_3d
  use eb_reactive_redistribution_3d_mod, only: reactive_eb_weighted_state_redistribute_3d
#ifdef PELEF_TEST_MPI_REDIST
  use mpi_f08
  use mpi_eb_tree_redistribution_3d_mod, only: redistribute_mpi_eb_tree_state_3d, &
    redistribute_sdc_mpi_eb_tree_state_level_rate_3d
#endif
  implicit none
  type(eb_reactive_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,coarse
  type(amr_patch_3d) :: boxes(2)
  type(amr_tree_state_3d) :: raw,result,saved,empty,uniform
  real(dp), allocatable :: flat(:,:,:,:),expected(:,:,:,:)
  real(dp) :: before(2),after(2),error,x,y,z,weight,minimum(2),maximum(2)
  integer :: layout,p,n(3),i,j,k,offset,children
#ifdef PELEF_TEST_MPI_REDIST
  integer :: rank,nranks,status
#endif
  logical :: ok
#ifdef PELEF_TEST_MPI_REDIST
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
#endif
  call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    'x',.49_dp,finest,ok)
  call require(ok,'fine geometry')
  call coarsen_eb_geometry_3d(finest,2,coarse,ok)
  call require(ok,'coarse geometry')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),ok)
  call require(ok,'left patch')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,boxes(2),ok)
  call require(ok,'right patch')
  allocate(flat(2,8,8,8),expected(2,8,8,8))
  do k=1,8
  do j=1,8
  do i=1,8
    flat(:,i,j,k)=value((i-.5_dp)/8,(j-.5_dp)/8,(k-.5_dp)/8)
    if(finest%volume_fraction(i,j,k)<=0) flat(:,i,j,k)=0
  end do
  end do
  end do
  call reactive_eb_weighted_state_redistribute_3d(finest,flat,expected,ok)
  call require(ok,'single grid reference')
  do layout=0,3
    children=min(layout,2)
    if(layout==3) children=1
    if(layout==3) then
      call initialize_amr_patch_tree_3d([4,4,4],[1],boxes(2:2),tree%topology,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[(1,p=1,children)],boxes(:children),tree%topology,ok)
    end if
    call require(ok,'topology')
    allocate(tree%nodes(children+1),raw%fields(children+1))
    tree%nodes(1)%geometry=coarse
    do p=2,children+1
      offset=4*(p-2)
      if(layout==3) offset=4
      call extract_eb_geometry_3d(finest,[offset+1,1,1],[offset+4,8,8],tree%nodes(p)%geometry,ok)
      call require(ok,'child geometry')
    end do
    before=0;minimum=huge(1._dp);maximum=-huge(1._dp)
    do p=1,children+1
      n=tree%topology%nodes(p)%shape
      allocate(tree%nodes(p)%state(2,n(1),n(2),n(3)),source=0._dp)
      allocate(tree%nodes(p)%temperature(n(1),n(2),n(3)),source=300._dp)
      allocate(raw%fields(p)%values(2,n(1),n(2),n(3)))
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(p==1) then
          x=(i-.5_dp)/4;y=(j-.5_dp)/4;z=(k-.5_dp)/4
        else
          offset=4*(p-2)
          if(layout==3) offset=4
          x=(i+offset-.5_dp)/8;y=(j-.5_dp)/8;z=(k-.5_dp)/8
        end if
        raw%fields(p)%values(:,i,j,k)=value(x,y,z)
        weight=fluid_volume(p,i,j,k)
        before=before+weight*raw%fields(p)%values(:,i,j,k)
        if(weight>0) then
          minimum=min(minimum,raw%fields(p)%values(:,i,j,k))
          maximum=max(maximum,raw%fields(p)%values(:,i,j,k))
        end if
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) &
          raw%fields(p)%values(:,i,j,k)=ieee_value(0._dp,ieee_quiet_nan)
      end do
      end do
      end do
    end do
    call redistribute_eb_tree_state_3d(tree,raw,result,ok)
    call require(ok,'composite redistribution')
    call check_decomposition(raw,result)
#ifdef PELEF_TEST_MPI_REDIST
    call check_mpi(raw,result)
#endif
    after=0
    do p=1,children+1
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        after=after+fluid_volume(p,i,j,k)*result%fields(p)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    error=maxval(abs(after-before)/max(1._dp,abs(before)))
    print *, 'layout, composite residual: ',layout,error
    call require(error<3.e-14_dp,'signed component conservation')
    if(layout==0) then
      block
        real(dp) :: ref(2,4,4,4),input(2,4,4,4)
        input=raw%fields(1)%values
        do k=1,4
        do j=1,4
        do i=1,4
          if(coarse%volume_fraction(i,j,k)<=0) input(:,i,j,k)=0
        end do
        end do
        end do
        call reactive_eb_weighted_state_redistribute_3d(coarse,input,ref,ok)
        call require(ok,'root reference')
        call require(maxval(abs(result%fields(1)%values-ref))<2.e-14_dp,'single patch equivalence')
      end block
    else if(layout==2) then
      do p=2,3
        offset=4*(p-2)
        call require(maxval(abs(result%fields(p)%values-expected(:,offset+1:offset+4,:,:)))<2.e-14_dp, &
          'peer split equals single grid')
      end do
    else
      call require(any(result%fields(1)%values(:,3,2,2)/=raw%fields(1)%values(:,3,2,2)).or. &
        any(result%fields(1)%values(:,2,2,2)/=raw%fields(1)%values(:,2,2,2)), &
        'coarse/fine redistribution changes coarse state')
    end if
    do p=1,children+1
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(fluid_volume(p,i,j,k)<=0) cycle
        call require(all(result%fields(p)%values(:,i,j,k)>=minimum-2.e-14_dp).and. &
          all(result%fields(p)%values(:,i,j,k)<=maximum+2.e-14_dp),'state remains within input bounds')
      end do
      end do
      end do
    end do
    uniform=raw
    do p=1,children+1
      uniform%fields(p)%values=2.5_dp
    end do
    call redistribute_eb_tree_state_3d(tree,uniform,result,ok)
    call require(ok,'uniform composite state')
#ifdef PELEF_TEST_MPI_REDIST
    call check_mpi(uniform,result)
#endif
    do p=1,children+1
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(fluid_volume(p,i,j,k)<=0) cycle
        call require(maxval(abs(result%fields(p)%values(:,i,j,k)-2.5_dp))<2.e-14_dp,'constant preservation')
      end do
      end do
      end do
    end do
    saved=result
    call redistribute_eb_tree_state_3d(tree,uniform,result,ok,target_volume_fraction=0._dp)
    call require(.not.ok,'invalid target rejected')
    raw%fields(1)%values(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
    call redistribute_eb_tree_state_3d(tree,raw,result,ok)
    call require(.not.ok,'nonfinite active rate rejected')
    do p=1,children+1
      call require(all(result%fields(p)%values==saved%fields(p)%values),'failure retains output')
    end do
    call redistribute_eb_tree_state_3d(tree,empty,result,ok)
    call require(.not.ok,'absent fields rejected')
    deallocate(tree%nodes,raw%fields)
  end do
  print *, 'test_eb_tree_state_redistribution_3d: PASS'
#ifdef PELEF_TEST_MPI_REDIST
  call MPI_Finalize(status)
#endif
contains
  subroutine check_decomposition(input,reference)
    type(amr_tree_state_3d), intent(in) :: input,reference
    type(amr_tree_state_3d) :: delta,total,rate,transfer,combined,full
    real(dp) :: inventory(2),balance(2),h
    integer :: depth,q,a,b,c,d(3)
    logical :: good
#ifdef PELEF_TEST_MPI_REDIST
    type(amr_tree_state_3d) :: local,local_state,mrate,mtransfer
    type(eb_geometry_3d) :: metrics(size(tree%nodes))
    integer :: owners(size(tree%nodes)),distribution
#endif
    total=input;combined=input
    do q=1,size(tree%nodes)
      total%fields(q)%values=0;combined%fields(q)%values=0
      tree%nodes(q)%state=input%fields(q)%values
    end do
    h=.25_dp
    call redistribute_sdc_eb_tree_state_rate_3d(tree,input,h,full,good)
    call require(good,'full affine state rate')
    do depth=0,tree%topology%level_count()-1
      call redistribute_eb_tree_state_3d(tree,input,delta,good,correction_depth=depth)
      call require(good,'selected conservative neighborhood correction')
#ifdef PELEF_TEST_MPI_REDIST
      call check_mpi(input,delta,depth=depth)
#endif
      call redistribute_sdc_eb_tree_state_level_rate_3d(tree,input,depth,h,rate,transfer,good)
      call require(good,'selected affine rate and integrated transfer')
      inventory=0;balance=0
      do q=1,size(tree%nodes)
        total%fields(q)%values=total%fields(q)%values+delta%fields(q)%values
        combined%fields(q)%values=combined%fields(q)%values+rate%fields(q)%values+transfer%fields(q)%values/h
        if(tree%topology%nodes(q)%level==depth) then
          call require(all(transfer%fields(q)%values==0),'selected transfer zero')
        else
          call require(all(rate%fields(q)%values==0),'inactive rate zero')
        end if
        d=tree%topology%nodes(q)%shape
        do c=1,d(3)
        do b=1,d(2)
        do a=1,d(1)
          if(fluid_volume(q,a,b,c)<=0) cycle
          inventory=inventory+fluid_volume(q,a,b,c)*delta%fields(q)%values(:,a,b,c)
          balance=balance+fluid_volume(q,a,b,c)*(h*rate%fields(q)%values(:,a,b,c)+transfer%fields(q)%values(:,a,b,c))
          if(tree%topology%nodes(q)%level==depth) &
            balance=balance-h*fluid_volume(q,a,b,c)*input%fields(q)%values(:,a,b,c)
        end do
        end do
        end do
      end do
      call require(maxval(abs(inventory))<3.e-14_dp,'each level baseline correction conserves')
      call require(maxval(abs(balance))<3.e-14_dp,'each level affine rate conserves with transfers')
#ifdef PELEF_TEST_MPI_REDIST
      do distribution=1,2
        owners=[(mod(q-1,nranks),q=1,size(tree%nodes))]
        if(distribution==2) owners=0
        allocate(local%fields(size(tree%nodes)),local_state%fields(size(tree%nodes)))
        do q=1,size(tree%nodes)
          metrics(q)=tree%nodes(q)%geometry
          if(owners(q)/=rank) cycle
          local%fields(q)%values=input%fields(q)%values
          local_state%fields(q)%values=tree%nodes(q)%state
        end do
        call redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree%topology,owners,metrics,local_state,local, &
          2,depth,h,mrate,mtransfer,MPI_COMM_WORLD,good)
        call require(good,'MPI level state rate')
        do q=1,size(tree%nodes)
          call require(allocated(mrate%fields(q)%values).eqv.(owners(q)==rank),'owner-local level rate')
          call require(allocated(mtransfer%fields(q)%values).eqv.(owners(q)==rank),'owner-local level transfer')
          if(owners(q)/=rank) cycle
          call require(all(mrate%fields(q)%values==rate%fields(q)%values),'exact serial MPI affine level rate')
          call require(all(mtransfer%fields(q)%values==transfer%fields(q)%values),'exact serial MPI level transfer')
        end do
        deallocate(local%fields,local_state%fields)
      end do
#endif
    end do
    do q=1,size(tree%nodes)
      d=tree%topology%nodes(q)%shape
      call require(maxval(abs(combined%fields(q)%values-full%fields(q)%values))<3.e-14_dp, &
        'sum of level rates and transfers equals full affine operator')
      do c=1,d(3)
      do b=1,d(2)
      do a=1,d(1)
        if(tree%nodes(q)%geometry%volume_fraction(a,b,c)<=0) cycle
        call require(maxval(abs(total%fields(q)%values(:,a,b,c)- &
          (reference%fields(q)%values(:,a,b,c)-input%fields(q)%values(:,a,b,c))))<3.e-14_dp, &
          'sum of group corrections reconstructs full state displacement')
      end do
      end do
      end do
    end do
  end subroutine check_decomposition
#ifdef PELEF_TEST_MPI_REDIST
  subroutine check_mpi(input,reference,labels,depth)
    type(amr_tree_state_3d), intent(in) :: input,reference
    character(len=*), optional, intent(in) :: labels(6)
    integer, optional, intent(in) :: depth
    type(amr_tree_state_3d) :: local,output,previous
    type(eb_geometry_3d) :: geometry(size(tree%nodes))
    integer :: owners(size(tree%nodes)),q,distribution
    character(len=16) :: faces(6),bad_faces(6)
    logical :: good
    real(dp) :: target
    faces='outflow'
    if(present(labels)) faces=labels
    do q=1,size(tree%nodes)
      geometry(q)=tree%nodes(q)%geometry
    end do
    do distribution=1,2
      owners=[(mod(q-1,nranks),q=1,size(tree%nodes))]
      if(distribution==2) owners=0
      allocate(local%fields(size(tree%nodes)))
      do q=1,size(tree%nodes)
        if(owners(q)==rank) local%fields(q)%values=input%fields(q)%values
      end do
      call redistribute_mpi_eb_tree_state_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,boundary_faces=faces,correction_depth=depth)
      call require(good,'MPI composite redistribution')
      do q=1,size(tree%nodes)
        call require(allocated(output%fields(q)%values).eqv.(owners(q)==rank),'MPI owner-local result')
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==reference%fields(q)%values),'exact serial/MPI redistribution')
      end do
      previous=output;target=.5_dp
      if(rank==nranks-1) target=.75_dp
      call redistribute_mpi_eb_tree_state_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,target_volume_fraction=target)
      call require(.not.good,'rank-dependent target rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==previous%fields(q)%values),'MPI target rollback')
      end do
      bad_faces=faces
      if(rank==nranks-1) bad_faces(1:2)='slip_wall'
      call redistribute_mpi_eb_tree_state_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,boundary_faces=bad_faces)
      call require(.not.good,'rank-dependent redistribution boundaries rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==previous%fields(q)%values),'MPI boundary rollback')
      end do
      if(owners(1)==rank) local%fields(1)%values(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
      call redistribute_mpi_eb_tree_state_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,boundary_faces=faces,correction_depth=depth)
      call require(.not.good,'one-owner nonfinite rate rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==previous%fields(q)%values),'MPI source rollback')
      end do
      deallocate(local%fields)
    end do
  end subroutine check_mpi
#endif
  function value(x,y,z) result(v)
    real(dp), intent(in) :: x,y,z
    real(dp) :: v(2)
    v=[x*x+3*y-z,2+x-y+z*z]
  end function value
  real(dp) function fluid_volume(p,i,j,k) result(v)
    integer, intent(in) :: p,i,j,k
    v=0
    if(p==1.and.layout==1.and.i<=2) return
    if(p==1.and.layout==2) return
    if(p==1.and.layout==3.and.i>=3) return
    associate(g=>tree%nodes(p)%geometry)
      v=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
    end associate
  end function fluid_volume
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(condition) return
    print *, 'FAIL: ',message
#ifdef PELEF_TEST_MPI_REDIST
    call MPI_Abort(MPI_COMM_WORLD,1,status)
#endif
    error stop 1
  end subroutine require
end program test_eb_tree_state_redistribution_3d
