program test_eb_tree_state_redistribution_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,redistribute_eb_tree_state_3d, &
    redistribute_sdc_eb_tree_state_level_rate_3d,redistribute_sdc_eb_tree_state_rate_3d
  use eb_reactive_redistribution_3d_mod, only: reactive_eb_weighted_state_redistribute_3d, &
    state_redist_register_3d,compute_state_redist_level_stage_3d
#ifdef PELEF_TEST_MPI_REDIST
  use mpi_f08
  use mpi_eb_tree_redistribution_3d_mod, only: compute_state_redist_mpi_level_stage_3d,redistribute_mpi_eb_tree_state_3d, &
    redistribute_sdc_mpi_eb_tree_state_level_rate_3d,redistribute_sdc_mpi_eb_tree_state_rate_3d
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
  call check_prepared_level_stages()
  call check_boundary_reconstruction()
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
    call check_high_order_tree(raw)
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
  subroutine check_prepared_level_stages()
    type(amr_patch_tree_3d) :: topology
    type(amr_patch_3d) :: patches(2)
    type(eb_geometry_3d) :: wide(3)
    type(amr_tree_state_3d) :: old,rhs,reference,scales,local_scales,local_old,local_rhs,local_rate,kept_rate
    type(state_redist_register_3d), allocatable :: expected_regs(:),regs(:),kept_regs(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: dx(3),lo(3),hi(3),local_dt
    integer :: node,nn(3),dims(3),ii,jj,kk,mode,owners(3),myrank,number
    logical :: good
#ifdef PELEF_TEST_MPI_REDIST
    myrank=rank;number=nranks
#else
    myrank=0;number=1
#endif
    call initialize_amr_patch_3d(6,4,4,1,2,2,3,2,3,2,patches(1),good)
    call require(good,'prepared stage left patch')
    call initialize_amr_patch_3d(6,4,4,3,4,2,3,2,3,2,patches(2),good)
    call require(good,'prepared stage right patch')
    call initialize_amr_patch_tree_3d([6,4,4],[1,1],patches,topology,good)
    call require(good,'prepared stage tree')
    call amr_tree_integer_layout_3d(topology,origin,domain,ratios,good)
    call require(good,'prepared stage integer layout')
    allocate(old%fields(3),rhs%fields(3),reference%fields(3),scales%fields(3))
    do node=1,3
      nn=topology%nodes(node)%shape;dims=nn+12
      dx=1._dp/real(domain(:,topology%nodes(node)%level+1),dp)
      lo=(real(origin(:,node),dp)-6._dp)*dx;hi=lo+real(dims,dp)*dx
      call build_axis_plane_eb_geometry_3d(dims(1),dims(2),dims(3),lo(1),hi(1),lo(2),hi(2),lo(3),hi(3), &
        'y',.71_dp,wide(node),good)
      call require(good,'prepared stage wide geometry')
      allocate(old%fields(node)%values(2,dims(1),dims(2),dims(3)), &
        rhs%fields(node)%values(2,dims(1),dims(2),dims(3)),reference%fields(node)%values(2,nn(1),nn(2),nn(3)))
      allocate(scales%fields(node)%values(1,dims(1),dims(2),dims(3)))
      do kk=1,dims(3)
      do jj=1,dims(2)
      do ii=1,dims(1)
        scales%fields(node)%values(1,ii,jj,kk)=.5_dp+.125_dp*modulo(sum(origin(:,node))+ii+jj+kk,5_int64)
        old%fields(node)%values(:,ii,jj,kk)=value(lo(1)+(ii-.5_dp)*dx(1),lo(2)+(jj-.5_dp)*dx(2), &
          lo(3)+(kk-.5_dp)*dx(3))
        rhs%fields(node)%values(:,ii,jj,kk)=[.1_dp,-.2_dp]+.01_dp*(lo(2)+(jj-.5_dp)*dx(2))
      enddo
      enddo
      enddo
      reference%fields(node)%values=-8._dp
    enddo
    call compute_state_redist_level_stage_3d(topology,wide,old,rhs,2,1,.01_dp,reference,expected_regs,good, &
      max_order=4,register_scale=.5_dp,update_scales=scales)
    call require(good,'prepared stage serial control')
    call require(maxval(abs(expected_regs(2)%fine_delta))>1.e-8_dp,'prepared stage exercises fine ghost exchange')
    do mode=0,2
      owners=0
      if(mode==1) owners=[0,modulo(1,number),modulo(2,number)]
      if(mode==2) owners=[modulo(1,number),0,modulo(1,number)]
      allocate(local_old%fields(3),local_rhs%fields(3),local_rate%fields(3),local_scales%fields(3))
      do node=1,3
        if(owners(node)/=myrank) cycle
        local_old%fields(node)%values=old%fields(node)%values
        local_rhs%fields(node)%values=rhs%fields(node)%values
        local_scales%fields(node)%values=scales%fields(node)%values
        nn=topology%nodes(node)%shape
        allocate(local_rate%fields(node)%values(2,nn(1),nn(2),nn(3)),source=-8._dp)
      enddo
#ifdef PELEF_TEST_MPI_REDIST
      call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
        local_rate,regs,MPI_COMM_WORLD,good,max_order=4,register_scale=.5_dp,update_scales=local_scales)
#else
      call compute_state_redist_level_stage_3d(topology,wide,local_old,local_rhs,2,1,.01_dp, &
        local_rate,regs,good,max_order=4,register_scale=.5_dp,owned=owners==myrank,update_scales=local_scales)
#endif
      call require(good,'prepared stage owner-local update')
      do node=1,3
        if(owners(node)/=myrank) then
          call require(.not.allocated(local_rate%fields(node)%values),'prepared stage no remote gas')
          call require(.not.allocated(regs(node)%fine_delta),'prepared stage no remote ledger')
          cycle
        endif
        call require(all(local_rate%fields(node)%values==reference%fields(node)%values),'prepared stage serial rate equality')
        if(node==1) cycle
        call require(all(regs(node)%fine_delta==expected_regs(node)%fine_delta),'prepared stage serial ledger equality')
      enddo
      kept_rate=local_rate;kept_regs=regs
      if(owners(3)==myrank) local_rhs%fields(3)%values(:,7,10,7)=ieee_value(1._dp,ieee_quiet_nan)
#ifdef PELEF_TEST_MPI_REDIST
      call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
        local_rate,regs,MPI_COMM_WORLD,good,max_order=4,register_scale=.5_dp,update_scales=local_scales)
#else
      call compute_state_redist_level_stage_3d(topology,wide,local_old,local_rhs,2,1,.01_dp, &
        local_rate,regs,good,max_order=4,register_scale=.5_dp,owned=owners==myrank,update_scales=local_scales)
#endif
      call require(.not.good,'prepared stage remote failure collectively rejected')
      do node=1,3
        if(owners(node)/=myrank) cycle
        call require(all(local_rate%fields(node)%values==kept_rate%fields(node)%values),'prepared stage collective rate rollback')
        if(node==1) cycle
        call require(all(regs(node)%fine_delta==kept_regs(node)%fine_delta),'prepared stage collective ledger rollback')
      enddo
#ifdef PELEF_TEST_MPI_REDIST
      if(owners(3)==myrank) then
        local_rhs%fields(3)%values=rhs%fields(3)%values
        local_scales%fields(3)%values(1,7,10,7)=0
      endif
      call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
        local_rate,regs,MPI_COMM_WORLD,good,max_order=4,register_scale=.5_dp,update_scales=local_scales)
      call require(.not.good,'prepared stage invalid owner scale rejected collectively')
      do node=2,3
        if(owners(node)/=myrank) cycle
        call require(all(local_rate%fields(node)%values==kept_rate%fields(node)%values),'prepared scale rate rollback')
        call require(all(regs(node)%fine_delta==kept_regs(node)%fine_delta),'prepared scale ledger rollback')
      enddo
      if(owners(3)==myrank) local_scales%fields(3)%values=scales%fields(3)%values
      if(number>1) then
        local_dt=.01_dp
        if(myrank==0) local_dt=.02_dp
        call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,local_dt, &
          local_rate,regs,MPI_COMM_WORLD,good,max_order=4,register_scale=.5_dp,update_scales=local_scales)
        call require(.not.good,'prepared stage inconsistent dt rejected')
        if(myrank==0) then
          call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
            local_rate,regs,MPI_COMM_WORLD,good,max_order=0)
        else
          call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
            local_rate,regs,MPI_COMM_WORLD,good)
        endif
        call require(.not.good,'prepared stage inconsistent optional presence rejected')
        if(myrank==0) then
          call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
            local_rate,regs,MPI_COMM_WORLD,good,update_scales=local_scales)
        else
          call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,local_old,local_rhs,2,1,.01_dp, &
            local_rate,regs,MPI_COMM_WORLD,good)
        endif
        call require(.not.good,'prepared stage inconsistent update scale presence rejected')
      endif
#endif
      deallocate(local_old%fields,local_rhs%fields,local_rate%fields,local_scales%fields,regs)
    enddo
  end subroutine check_prepared_level_stages
  subroutine check_boundary_reconstruction()
    type(eb_reactive_tree_3d) :: local_tree
    type(amr_tree_state_3d) :: input,output,halos
    type(eb_geometry_3d) :: hg(1)
    type(amr_patch_3d) :: no_children(0)
    real(dp) :: ref(1,8,8,8),x,y,z,before,after
    integer :: a,b,c,da,db,dc,order
    character(len=16) :: faces(6)
#ifdef PELEF_TEST_MPI_REDIST
    type(amr_tree_state_3d) :: local_input,local_halos,mpi_output
    integer :: owners(1),distribution
#endif
    logical :: good
    call initialize_amr_patch_tree_3d([8,8,8],[integer::],no_children,local_tree%topology,good)
    call require(good,'boundary tree topology')
    allocate(local_tree%nodes(1),input%fields(1),halos%fields(1))
    call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',.37_dp,local_tree%nodes(1)%geometry,good)
    call require(good,'boundary tree geometry')
    call build_axis_plane_eb_geometry_3d(12,12,12,-.25_dp,1.25_dp,-.25_dp,1.25_dp,-.25_dp,1.25_dp, &
      'x',.37_dp,hg(1),good)
    call require(good,'boundary tree halo geometry')
    allocate(local_tree%nodes(1)%state(1,8,8,8),source=0._dp)
    allocate(local_tree%nodes(1)%temperature(8,8,8),source=300._dp)
    allocate(input%fields(1)%values(1,8,8,8),halos%fields(1)%values(1,12,12,12))
    do c=1,8
    do b=1,8
    do a=1,8
      associate(g=>local_tree%nodes(1)%geometry)
        x=(a-.5_dp+g%cell_centroid_x(a,b,c))/8;y=(b-.5_dp)/8;z=(c-.5_dp)/8
        input%fields(1)%values(1,a,b,c)=2+.1_dp*x+.2_dp*y*y+.03_dp*z
        if(g%volume_fraction(a,b,c)<=0) input%fields(1)%values(:,a,b,c)=0
      end associate
    enddo
    enddo
    enddo
    do c=1,12
    do b=1,12
    do a=1,12
      da=max(1,min(8,a-2));db=max(1,min(8,b-2));dc=max(1,min(8,c-2))
      halos%fields(1)%values(:,a,b,c)=input%fields(1)%values(:,da,db,dc)
      if(b<=2.and.local_tree%nodes(1)%geometry%volume_fraction(da,db,dc)>0) then
        x=(da-.5_dp+local_tree%nodes(1)%geometry%cell_centroid_x(da,db,dc))/8;z=(dc-.5_dp)/8
        halos%fields(1)%values(1,a,b,c)=2+.1_dp*x+.03_dp*z
      endif
    enddo
    enddo
    enddo
    faces='outflow';faces(3)='inflow'
    before=sum(local_tree%nodes(1)%geometry%volume_fraction*input%fields(1)%values(1,:,:,:))
    do order=2,4,2
      call reactive_eb_weighted_state_redistribute_3d(local_tree%nodes(1)%geometry,input%fields(1)%values,ref, &
        good,max_order=order,boundary_state=halos%fields(1)%values,boundary_geometry=hg(1),boundary_faces=faces)
      call require(good,'boundary tree single-grid reference')
      call redistribute_eb_tree_state_3d(local_tree,input,output,good,boundary_faces=faces,max_order=order, &
        boundary_halos=halos,boundary_geometry=hg)
      call require(good,'boundary tree reconstruction')
      call require(maxval(abs(output%fields(1)%values-ref))<2.e-12_dp,'boundary tree single-grid agreement')
      after=sum(local_tree%nodes(1)%geometry%volume_fraction*output%fields(1)%values(1,:,:,:))
      call require(abs(after-before)<2.e-12_dp,'boundary tree conservative physical recipients')
#ifdef PELEF_TEST_MPI_REDIST
      do distribution=0,min(1,nranks-1)
        owners=distribution
        allocate(local_input%fields(1),local_halos%fields(1))
        if(rank==owners(1)) then
          local_input%fields(1)%values=input%fields(1)%values
          local_halos%fields(1)%values=halos%fields(1)%values
        endif
        call redistribute_mpi_eb_tree_state_3d(local_tree%topology,owners,[local_tree%nodes(1)%geometry], &
          local_input,1,mpi_output,MPI_COMM_WORLD,good,boundary_faces=faces,max_order=order, &
          boundary_halos=local_halos,boundary_geometry=hg)
        call require(good,'MPI physical boundary Qhat with idle rank')
        if(rank==owners(1)) call require(maxval(abs(mpi_output%fields(1)%values-output%fields(1)%values)) &
          <2.e-12_dp,'MPI serial boundary Qhat reconstruction agreement')
        if(nranks>1) then
          if(rank==0) then
            call redistribute_mpi_eb_tree_state_3d(local_tree%topology,owners,[local_tree%nodes(1)%geometry], &
              local_input,1,mpi_output,MPI_COMM_WORLD,good,boundary_faces=faces,max_order=order)
          else
            call redistribute_mpi_eb_tree_state_3d(local_tree%topology,owners,[local_tree%nodes(1)%geometry], &
              local_input,1,mpi_output,MPI_COMM_WORLD,good,boundary_faces=faces,max_order=order, &
              boundary_halos=local_halos,boundary_geometry=hg)
          endif
          call require(.not.good,'MPI rank-disagreeing boundary Qhat rejection')
          if(rank==owners(1)) call require(maxval(abs(mpi_output%fields(1)%values-output%fields(1)%values)) &
            <2.e-12_dp,'MPI boundary Qhat collective failure rollback')
        endif
        deallocate(local_input%fields,local_halos%fields)
      enddo
#endif

    enddo
    ref=output%fields(1)%values
    call redistribute_eb_tree_state_3d(local_tree,input,output,good,max_order=2,boundary_halos=halos)
    call require(.not.good.and.all(output%fields(1)%values==ref),'boundary tree incomplete metadata rollback')
    hg(1)%x_lower=hg(1)%x_lower+.001_dp;hg(1)%x_upper=hg(1)%x_upper+.001_dp
    call redistribute_eb_tree_state_3d(local_tree,input,output,good,max_order=2, &
      boundary_halos=halos,boundary_geometry=hg)
    call require(.not.good.and.all(output%fields(1)%values==ref),'boundary tree metric mismatch rollback')
  end subroutine check_boundary_reconstruction
  subroutine check_high_order_tree(input)
    type(amr_tree_state_3d), intent(in) :: input
    type(amr_tree_state_3d) :: affine,output,previous,delta,combined,full,rate,transfer
    real(dp) :: original(2),advanced(2),position(3),h,lower(2),upper(2)
    real(dp), allocatable :: grid_input(:,:,:,:),grid_reference(:,:,:,:)
    character(len=16) :: periodic_faces(6)
    integer :: q,a,b,c,d(3),order,depth
#ifdef PELEF_TEST_MPI_REDIST
    type(amr_tree_state_3d) :: local,local_state,mrate,mtransfer
    type(eb_geometry_3d) :: metrics(size(tree%nodes))
    integer :: owners(size(tree%nodes)),distribution
#endif
    logical :: good
    affine=input;h=.25_dp
    do q=1,size(tree%nodes)
      d=tree%topology%nodes(q)%shape
      do c=1,d(3)
      do b=1,d(2)
      do a=1,d(1)
        associate(g=>tree%nodes(q)%geometry)
          position=[g%x_lower,g%y_lower,g%z_lower]+[g%dx,g%dy,g%dz]* &
            (real([a,b,c],dp)-.5_dp+[g%cell_centroid_x(a,b,c),g%cell_centroid_y(a,b,c),g%cell_centroid_z(a,b,c)])
          affine%fields(q)%values(:,a,b,c)=[1+.1_dp*position(1)+.03_dp*position(2)-.02_dp*position(3), &
            2-.04_dp*position(1)+.01_dp*position(2)+.05_dp*position(3)]
          if(g%volume_fraction(a,b,c)<=0) affine%fields(q)%values(:,a,b,c)=0
        end associate
      enddo
      enddo
      enddo
    enddo
    do order=2,4,2
      call redistribute_eb_tree_state_3d(tree,affine,output,good,max_order=order)
      call require(good,'high-order affine composite redistribution')
#ifdef PELEF_TEST_MPI_REDIST
      call check_mpi(affine,output,order=order)
#endif
      do q=1,size(tree%nodes)
        d=tree%topology%nodes(q)%shape
        do c=1,d(3)
        do b=1,d(2)
        do a=1,d(1)
          if(fluid_volume(q,a,b,c)<=0) cycle
          call require(maxval(abs(output%fields(q)%values(:,a,b,c)-affine%fields(q)%values(:,a,b,c))) &
            <2.e-12_dp,'high-order composite affine preservation')
        enddo
        enddo
        enddo
      enddo
      previous=output
      call redistribute_eb_tree_state_3d(tree,affine,output,good,max_order=1)
      call require(.not.good,'bad tree reconstruction order rejected')
      do q=1,size(tree%nodes)
        call require(all(output%fields(q)%values==previous%fields(q)%values),'bad tree order rollback')
      enddo
      call redistribute_eb_tree_state_3d(tree,input,output,good,max_order=order)
      call require(good,'high-order nonlinear composite redistribution')
#ifdef PELEF_TEST_MPI_REDIST
      call check_mpi(input,output,order=order)
#endif
      if(layout==0.or.layout==2) then
        if(layout==0) then
          grid_input=input%fields(1)%values
          do c=1,4
          do b=1,4
          do a=1,4
            if(coarse%volume_fraction(a,b,c)<=0) grid_input(:,a,b,c)=0
          enddo
          enddo
          enddo
          allocate(grid_reference(2,4,4,4))
          call reactive_eb_weighted_state_redistribute_3d(coarse,grid_input,grid_reference,good,max_order=order)
          call require(good.and.maxval(abs(output%fields(1)%values-grid_reference))<2.e-12_dp, &
            'high-order single patch equals grid reconstruction')
        else
          grid_input=flat;allocate(grid_reference(2,8,8,8))
          call reactive_eb_weighted_state_redistribute_3d(finest,grid_input,grid_reference,good,max_order=order)
          call require(good,'high-order fine grid control')
          do q=2,3
            call require(maxval(abs(output%fields(q)%values-grid_reference(:,4*(q-2)+1:4*(q-1),:,:)))<2.e-12_dp, &
              'high-order peer split equals grid reconstruction')
          enddo
        endif
        deallocate(grid_input,grid_reference)
      endif
      original=0;advanced=0;lower=huge(1._dp);upper=-huge(1._dp)
      do q=1,size(tree%nodes)
        d=tree%topology%nodes(q)%shape
        do c=1,d(3)
        do b=1,d(2)
        do a=1,d(1)
          if(fluid_volume(q,a,b,c)<=0) cycle
          original=original+fluid_volume(q,a,b,c)*input%fields(q)%values(:,a,b,c)
          advanced=advanced+fluid_volume(q,a,b,c)*output%fields(q)%values(:,a,b,c)
          lower=min(lower,input%fields(q)%values(:,a,b,c));upper=max(upper,input%fields(q)%values(:,a,b,c))
        enddo
        enddo
        enddo
      enddo
      call require(maxval(abs(original-advanced))<3.e-14_dp,'high-order tree component inventory')
      do q=1,size(tree%nodes)
        d=tree%topology%nodes(q)%shape
        do c=1,d(3)
        do b=1,d(2)
        do a=1,d(1)
          if(fluid_volume(q,a,b,c)<=0) cycle
          call require(all(output%fields(q)%values(:,a,b,c)>=lower-2.e-12_dp).and. &
            all(output%fields(q)%values(:,a,b,c)<=upper+2.e-12_dp),'high-order tree component bounds')
        enddo
        enddo
        enddo
      enddo
      ! Nonlinear limited SDC corrections must partition one provisional state.
      do q=1,size(tree%nodes)
        tree%nodes(q)%state=affine%fields(q)%values
      enddo
      call redistribute_sdc_eb_tree_state_rate_3d(tree,input,h,full,good,max_order=order)
      call require(good,'high-order full SDC tree rate')
      combined=input
      do q=1,size(tree%nodes)
        combined%fields(q)%values=0
      enddo
      do depth=0,tree%topology%level_count()-1
        call redistribute_sdc_eb_tree_state_level_rate_3d(tree,input,depth,h,rate,transfer,good,max_order=order)
        call require(good,'high-order SDC tree level partition')
#ifdef PELEF_TEST_MPI_REDIST
        do distribution=1,2
          owners=[(mod(q-1,nranks),q=1,size(tree%nodes))]
          if(distribution==2) owners=0
          allocate(local%fields(size(tree%nodes)),local_state%fields(size(tree%nodes)))
          do q=1,size(tree%nodes)
            metrics(q)=tree%nodes(q)%geometry
            if(owners(q)/=rank) cycle
            local%fields(q)%values=input%fields(q)%values
            local_state%fields(q)%values=affine%fields(q)%values
          enddo
          call redistribute_sdc_mpi_eb_tree_state_level_rate_3d(tree%topology,owners,metrics,local_state,local, &
            2,depth,h,mrate,mtransfer,MPI_COMM_WORLD,good,max_order=order)
          call require(good,'high-order MPI SDC level partition')
          do q=1,size(tree%nodes)
            call require(allocated(mrate%fields(q)%values).eqv.(owners(q)==rank),'high-order MPI local level rate')
            if(owners(q)/=rank) cycle
            call require(maxval(abs(mrate%fields(q)%values-rate%fields(q)%values))<2.e-13_dp, &
              'high-order serial/MPI SDC level rate agreement')
            call require(maxval(abs(mtransfer%fields(q)%values-transfer%fields(q)%values))<2.e-13_dp, &
              'high-order serial/MPI SDC level transfer agreement')
          enddo
          call redistribute_sdc_mpi_eb_tree_state_rate_3d(tree%topology,owners,metrics,local_state,local, &
            2,h,mrate,MPI_COMM_WORLD,good,max_order=order)
          call require(good,'high-order MPI full SDC rate')
          do q=1,size(tree%nodes)
            if(owners(q)/=rank) cycle
            call require(maxval(abs(mrate%fields(q)%values-full%fields(q)%values))<2.e-13_dp, &
              'high-order serial/MPI full SDC agreement')
          enddo
          deallocate(local%fields,local_state%fields)
        enddo
#endif
        do q=1,size(tree%nodes)
          combined%fields(q)%values=combined%fields(q)%values+rate%fields(q)%values+transfer%fields(q)%values/h
        enddo
      enddo
      do q=1,size(tree%nodes)
        call require(maxval(abs(combined%fields(q)%values-full%fields(q)%values))<3.e-14_dp, &
          'high-order SDC level sum equals full limited operator')
      enddo
      periodic_faces='outflow';periodic_faces(3:6)='periodic'
      call redistribute_eb_tree_state_3d(tree,input,output,good,boundary_faces=periodic_faces,max_order=order)
      call require(good,'high-order periodic graph redistribution')
#ifdef PELEF_TEST_MPI_REDIST
      call check_mpi(input,output,labels=periodic_faces,order=order)
#endif
      advanced=0
      do q=1,size(tree%nodes)
        d=tree%topology%nodes(q)%shape
        do c=1,d(3)
        do b=1,d(2)
        do a=1,d(1)
          if(fluid_volume(q,a,b,c)<=0) cycle
          advanced=advanced+fluid_volume(q,a,b,c)*output%fields(q)%values(:,a,b,c)
          call require(all(output%fields(q)%values(:,a,b,c)>=lower-2.e-12_dp).and. &
            all(output%fields(q)%values(:,a,b,c)<=upper+2.e-12_dp),'high-order periodic graph bounds')
        enddo
        enddo
        enddo
      enddo
      call require(maxval(abs(original-advanced))<3.e-14_dp,'high-order periodic graph inventory')
    enddo
  end subroutine check_high_order_tree

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
  subroutine check_mpi(input,reference,labels,depth,order)
    type(amr_tree_state_3d), intent(in) :: input,reference
    character(len=*), optional, intent(in) :: labels(6)
    integer, optional, intent(in) :: depth,order
    type(amr_tree_state_3d) :: local,output,previous
    type(eb_geometry_3d) :: geometry(size(tree%nodes))
    integer :: owners(size(tree%nodes)),q,distribution,bad_order
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
        MPI_COMM_WORLD,good,boundary_faces=faces,correction_depth=depth,max_order=order)
      call require(good,'MPI composite redistribution')
      do q=1,size(tree%nodes)
        call require(allocated(output%fields(q)%values).eqv.(owners(q)==rank),'MPI owner-local result')
        if(owners(q)/=rank) cycle
        if(present(order)) then
          call require(maxval(abs(output%fields(q)%values-reference%fields(q)%values))<2.e-13_dp, &
            'high-order serial/MPI redistribution agreement')
        else
          call require(all(output%fields(q)%values==reference%fields(q)%values),'exact serial/MPI redistribution')
        endif
      end do
      previous=output
      bad_order=0
      if(rank==nranks-1) bad_order=2
      call redistribute_mpi_eb_tree_state_3d(tree%topology,owners,geometry,local,2,output, &
        MPI_COMM_WORLD,good,max_order=bad_order)
      call require(.not.good,'rank-dependent reconstruction order rejected')
      do q=1,size(tree%nodes)
        if(owners(q)/=rank) cycle
        call require(all(output%fields(q)%values==previous%fields(q)%values),'MPI reconstruction order rollback')
      enddo
      target=.5_dp
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
        MPI_COMM_WORLD,good,boundary_faces=faces,correction_depth=depth,max_order=order)
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
