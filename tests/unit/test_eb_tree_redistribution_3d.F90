program test_eb_tree_redistribution_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,redistribute_sdc_eb_tree_flux_rate_3d, &
    redistribute_sdc_eb_tree_level_rate_3d
  use eb_reactive_redistribution_3d_mod, only: reactive_eb_flux_redistribute_3d
#ifdef PELEF_TEST_MPI_REDIST
  use mpi_f08
  use mpi_eb_tree_redistribution_3d_mod, only: redistribute_sdc_mpi_eb_tree_flux_rate_3d, &
    redistribute_sdc_mpi_eb_tree_level_rate_3d
#endif
  implicit none
  type(eb_reactive_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,coarse
  type(amr_patch_3d) :: boxes(2)
  type(amr_tree_state_3d) :: raw,result,saved,empty
  real(dp), allocatable :: flat(:,:,:,:),expected(:,:,:,:)
  real(dp) :: before(2),after(2),error,x,y,z,weight
  integer :: layout,p,n(3),i,j,k,offset,direction
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
    'y',.31_dp,finest,ok)
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
  call reactive_eb_flux_redistribute_3d(finest,flat,expected,ok)
  call require(ok,'single grid reference')
  do layout=0,2
    call initialize_amr_patch_tree_3d([4,4,4],[(1,p=1,layout)],boxes(:layout),tree%topology,ok)
    call require(ok,'topology')
    allocate(tree%nodes(layout+1),raw%fields(layout+1))
    tree%nodes(1)%geometry=coarse
    do p=2,layout+1
      offset=4*(p-2)
      call extract_eb_geometry_3d(finest,[offset+1,1,1],[offset+4,8,8],tree%nodes(p)%geometry,ok)
      call require(ok,'child geometry')
    end do
    before=0
    do p=1,layout+1
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
          x=(i+4*(p-2)-.5_dp)/8;y=(j-.5_dp)/8;z=(k-.5_dp)/8
        end if
        raw%fields(p)%values(:,i,j,k)=value(x,y,z)
        weight=fluid_volume(p,i,j,k)
        before=before+weight*raw%fields(p)%values(:,i,j,k)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) &
          raw%fields(p)%values(:,i,j,k)=ieee_value(0._dp,ieee_quiet_nan)
      end do
      end do
      end do
    end do
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,raw,result,ok)
    call require(ok,'composite redistribution')
    call check_levels(raw,result)
#ifdef PELEF_TEST_MPI_REDIST
    call check_mpi(raw,result)
#endif
    after=0
    do p=1,layout+1
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
        call reactive_eb_flux_redistribute_3d(coarse,input,ref,ok)
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
      call require(any(result%fields(1)%values(:,3,2,2)/=raw%fields(1)%values(:,3,2,2)), &
        'coarse interface receives redistribution')
      call require(any(result%fields(2)%values(:,4,3,3)/=raw%fields(2)%values(:,4,3,3)), &
        'fine interface receives redistribution')
    end if
    if(layout==1) then
      do direction=1,2
        do p=1,2
          raw%fields(p)%values=0
        end do
        if(direction==1) then
          raw%fields(1)%values(:,3,2,2)=[1._dp,-2._dp]
        else
          raw%fields(2)%values(:,4,3,3)=[1._dp,-2._dp]
        end if
        call redistribute_sdc_eb_tree_flux_rate_3d(tree,raw,result,ok)
        call require(ok,'coarse/fine impulse')
        if(direction==1) then
          call require(result%fields(2)%values(1,4,3,3)>0,'coarse impulse reaches fine leaf')
        else
          call require(result%fields(1)%values(1,3,2,2)>0,'fine impulse reaches coarse leaf')
        end if
        before=0;after=0
        do p=1,2
          n=tree%topology%nodes(p)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            weight=fluid_volume(p,i,j,k)
            before=before+weight*raw%fields(p)%values(:,i,j,k)
            after=after+weight*result%fields(p)%values(:,i,j,k)
          end do
          end do
          end do
        end do
        call require(maxval(abs(after-before))<3.e-14_dp,'coarse/fine impulse conservation')
      end do
    end if
    call check_periodic()
    saved=result
    raw%fields(1)%values(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,raw,result,ok)
    call require(.not.ok,'nonfinite active rate rejected')
    do p=1,layout+1
      call require(all(result%fields(p)%values==saved%fields(p)%values),'failure retains output')
    end do
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,empty,result,ok)
    call require(.not.ok,'absent fields rejected')
    deallocate(tree%nodes,raw%fields)
  end do
  print *, 'test_eb_tree_redistribution_3d: PASS'
#ifdef PELEF_TEST_MPI_REDIST
  call MPI_Finalize(status)
#endif
contains
  subroutine check_levels(input,reference,labels)
    type(amr_tree_state_3d), intent(in) :: input,reference
    character(len=*), optional, intent(in) :: labels(6)
    type(amr_tree_state_3d) :: rate,transfer,reconstructed,masked,previous_rate,previous_transfer
    real(dp) :: dt,weight,original(2),changed(2)
    integer :: depth,q,d(3),x,y,z
    logical :: good
    reconstructed=reference
    do q=1,size(tree%nodes)
      reconstructed%fields(q)%values=0
    end do
    do depth=0,tree%topology%level_count()-1
      dt=.25_dp/real(2**depth,dp)
      masked=input
      do q=1,size(tree%nodes)
        if(tree%topology%nodes(q)%level/=depth) &
          masked%fields(q)%values=ieee_value(0._dp,ieee_quiet_nan)
      end do
      call redistribute_sdc_eb_tree_level_rate_3d(tree,masked,depth,dt,rate,transfer,good,labels)
      call require(good,'source-level redistribution ignores other-level rates')
      original=0;changed=0
      do q=1,size(tree%nodes)
        reconstructed%fields(q)%values=reconstructed%fields(q)%values+ &
          rate%fields(q)%values+transfer%fields(q)%values/dt
        if(tree%topology%nodes(q)%level==depth) then
          call require(all(transfer%fields(q)%values==0),'no same-level deferred transfer')
        else
          call require(all(rate%fields(q)%values==0),'no other-level immediate rate')
        end if
        d=tree%topology%nodes(q)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          weight=fluid_volume(q,x,y,z)
          if(weight<=0) cycle
          if(tree%topology%nodes(q)%level==depth) original=original+weight*dt*input%fields(q)%values(:,x,y,z)
          changed=changed+weight*(dt*rate%fields(q)%values(:,x,y,z)+transfer%fields(q)%values(:,x,y,z))
        end do
        end do
        end do
      end do
      call require(maxval(abs(changed-original))<3.e-14_dp,'each source-level transfer balances')
#ifdef PELEF_TEST_MPI_REDIST
      call check_level_mpi(masked,depth,dt,rate,transfer,labels)
#endif
      previous_rate=rate;previous_transfer=transfer
      call redistribute_sdc_eb_tree_level_rate_3d(tree,masked,depth,-dt,rate,transfer,good,labels)
      call require(.not.good,'negative transfer interval rejected')
      do q=1,size(tree%nodes)
        call require(all(rate%fields(q)%values==previous_rate%fields(q)%values),'level rate rollback')
        call require(all(transfer%fields(q)%values==previous_transfer%fields(q)%values),'level transfer rollback')
      end do
    end do
    do q=1,size(tree%nodes)
      call require(maxval(abs(reconstructed%fields(q)%values-reference%fields(q)%values))<2.e-14_dp, &
        'sum of source-level contributions recovers composite operator')
    end do
  end subroutine check_levels
#ifdef PELEF_TEST_MPI_REDIST
  subroutine check_level_mpi(input,depth,dt,reference_rate,reference_transfer,labels)
    type(amr_tree_state_3d), intent(in) :: input,reference_rate,reference_transfer
    integer, intent(in) :: depth
    real(dp), intent(in) :: dt
    character(len=*), optional, intent(in) :: labels(6)
    type(amr_tree_state_3d) :: local,rate,transfer
    type(eb_geometry_3d) :: geometry(size(tree%nodes))
    integer :: owners(size(tree%nodes)),distribution,q,bad_depth
    real(dp) :: bad_dt
    logical :: good
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
      call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree%topology,owners,geometry,local,2,depth,dt, &
        rate,transfer,MPI_COMM_WORLD,good,labels)
      call require(good,'MPI source-level redistribution')
      do q=1,size(tree%nodes)
        call require(allocated(rate%fields(q)%values).eqv.(owners(q)==rank),'owner-local level rate')
        call require(allocated(transfer%fields(q)%values).eqv.(owners(q)==rank),'owner-local transfer')
        if(owners(q)/=rank) cycle
        call require(all(rate%fields(q)%values==reference_rate%fields(q)%values),'exact serial/MPI level rate')
        call require(all(transfer%fields(q)%values==reference_transfer%fields(q)%values),'exact serial/MPI transfer')
      end do
      bad_depth=depth;bad_dt=dt
      if(rank==nranks-1) bad_depth=depth+1
      call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree%topology,owners,geometry,local,2,bad_depth,dt, &
        rate,transfer,MPI_COMM_WORLD,good,labels)
      call require(.not.good,'rank-dependent source depth rejected')
      if(rank==nranks-1) bad_dt=2*dt
      call redistribute_sdc_mpi_eb_tree_level_rate_3d(tree%topology,owners,geometry,local,2,depth,bad_dt, &
        rate,transfer,MPI_COMM_WORLD,good,labels)
      call require(.not.good,'rank-dependent transfer interval rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(rate%fields(q)%values==reference_rate%fields(q)%values),'MPI level rate rollback')
        call require(all(transfer%fields(q)%values==reference_transfer%fields(q)%values),'MPI transfer rollback')
      end do
      deallocate(local%fields)
    end do
  end subroutine check_level_mpi
#endif
#ifdef PELEF_TEST_MPI_REDIST
  subroutine check_mpi(input,reference,labels)
    type(amr_tree_state_3d), intent(in) :: input,reference
    character(len=*), optional, intent(in) :: labels(6)
    type(amr_tree_state_3d) :: local,output,previous
    type(eb_geometry_3d) :: geometry(size(tree%nodes))
    integer :: owners(size(tree%nodes)),q,distribution
    character(len=16) :: faces(6),bad_faces(6)
    logical :: good
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
      call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,faces)
      call require(good,'MPI composite redistribution')
      do q=1,size(tree%nodes)
        call require(allocated(output%fields(q)%values).eqv.(owners(q)==rank),'MPI owner-local result')
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==reference%fields(q)%values),'exact serial/MPI redistribution')
      end do
      previous=output;bad_faces=faces
      if(rank==nranks-1) bad_faces(1:2)='slip_wall'
      call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,bad_faces)
      call require(.not.good,'rank-dependent redistribution boundaries rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==previous%fields(q)%values),'MPI boundary rollback')
      end do
      if(owners(1)==rank) local%fields(1)%values(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
      call redistribute_sdc_mpi_eb_tree_flux_rate_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,faces)
      call require(.not.good,'one-owner nonfinite rate rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==previous%fields(q)%values),'MPI source rollback')
      end do
      deallocate(local%fields)
    end do
  end subroutine check_mpi
#endif
  subroutine check_periodic()
    character(len=16) :: faces(6)
    type(amr_tree_state_3d) :: periodic,shifted,input
    real(dp), allocatable :: packed(:,:,:,:),reference(:,:,:,:)
    integer :: m,q,side
    faces='outflow';faces(1:2)='periodic'
    input=raw
    if(layout==0.or.layout==2) then
      m=4
      if(layout==2) m=8
      allocate(packed(2,m,m,m),reference(2,m,m,m))
      if(layout==0) then
        packed=input%fields(1)%values
      else
        packed(:,:4,:,:)=input%fields(2)%values
        packed(:,5:,:,:)=input%fields(3)%values
      end if
      call redistribute_sdc_eb_tree_flux_rate_3d(tree,input,periodic,ok,faces)
      call require(ok,'periodic redistribution')
      call check_levels(input,periodic,faces)
#ifdef PELEF_TEST_MPI_REDIST
      call check_mpi(input,periodic,faces)
#endif
      if(layout==0) then
        reference=periodic%fields(1)%values
      else
        reference(:,:4,:,:)=periodic%fields(2)%values
        reference(:,5:,:,:)=periodic%fields(3)%values
      end if
      packed=cshift(packed,1,dim=2);reference=cshift(reference,1,dim=2)
      if(layout==0) then
        input%fields(1)%values=packed
      else
        input%fields(2)%values=packed(:,:4,:,:)
        input%fields(3)%values=packed(:,5:,:,:)
      end if
      call redistribute_sdc_eb_tree_flux_rate_3d(tree,input,shifted,ok,faces)
      call require(ok,'translated periodic redistribution')
      if(layout==0) then
        packed=shifted%fields(1)%values
      else
        packed(:,:4,:,:)=shifted%fields(2)%values
        packed(:,5:,:,:)=shifted%fields(3)%values
      end if
      call require(maxval(abs(packed-reference))<2.e-14_dp,'periodic translation crosses patch and domain seams')
    else
      do side=1,2
        do q=1,2
          input%fields(q)%values=0
        end do
        if(side==1) then
          input%fields(1)%values(:,4,2,2)=[1._dp,-2._dp]
        else
          input%fields(2)%values(:,1,3,3)=[1._dp,-2._dp]
        end if
        call redistribute_sdc_eb_tree_flux_rate_3d(tree,input,periodic,ok,faces)
        call require(ok,'periodic coarse/fine impulse')
        call check_levels(input,periodic,faces)
#ifdef PELEF_TEST_MPI_REDIST
        call check_mpi(input,periodic,faces)
#endif
        if(side==1) then
          call require(periodic%fields(2)%values(1,1,3,3)>0,'periodic coarse-to-fine transfer')
        else
          call require(periodic%fields(1)%values(1,4,2,2)>0,'periodic fine-to-coarse transfer')
        end if
        before=0;after=0
        do q=1,2
          n=tree%topology%nodes(q)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            weight=fluid_volume(q,i,j,k)
            before=before+weight*input%fields(q)%values(:,i,j,k)
            after=after+weight*periodic%fields(q)%values(:,i,j,k)
          end do
          end do
          end do
        end do
        call require(maxval(abs(after-before))<3.e-14_dp,'periodic coarse/fine conservation')
      end do
    end if
    shifted=periodic
    faces(2)='outflow'
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,input,shifted,ok,faces)
    call require(.not.ok,'unpaired periodic faces rejected')
    do q=1,layout+1
      call require(all(shifted%fields(q)%values==periodic%fields(q)%values),'periodic control rollback')
    end do
    faces='outflow';faces(3:4)='periodic'
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,input,shifted,ok,faces)
    call require(.not.ok,'mismatched periodic geometry rejected')
  end subroutine check_periodic
  function value(x,y,z) result(v)
    real(dp), intent(in) :: x,y,z
    real(dp) :: v(2)
    v=[x*x+3*y-z,-2*x+y*y+z*z]
  end function value
  real(dp) function fluid_volume(p,i,j,k) result(v)
    integer, intent(in) :: p,i,j,k
    v=0
    if(p==1.and.layout>=1.and.i<=2) return
    if(p==1.and.layout==2) return
    associate(g=>tree%nodes(p)%geometry)
      v=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
    end associate
  end function fluid_volume
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(condition) return
    print *, 'FAIL: ',message
    error stop 1
  end subroutine require
end program test_eb_tree_redistribution_3d
