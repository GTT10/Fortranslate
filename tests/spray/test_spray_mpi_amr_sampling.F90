program test_spray_mpi_amr_sampling
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use spray_amr_tree_coupling_mod, only: sample_amr_tree_spray_field
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,sample_mpi_amr_particle_halo, &
    deposit_mpi_amr_particle_halo
  implicit none
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: full,local,snapshot
  real(dp), allocatable :: positions(:,:),samples(:,:),reference(:,:),saved(:,:),good_positions(:,:)
  real(dp) :: points(3,4),zero(3)=0._dp,spacing(3)=1._dp,bounds(3)
  integer :: rank,nranks,ierr,node,n(3),i,j,k,np,owners(3),wrong(3),level
  logical :: ok,periodic(3)=[.true.,.true.,.true.]
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(full%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(full%fields(node)%values(2,n(1),n(2),n(3)))
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      full%fields(node)%values(:,i,j,k)=[real(10*node+i+3*j+7*k,dp),real(node*i*j+k,dp)]
    end do
    end do
    end do
  end do
  owners=[0,1,nranks-1];call distribute()
  points(:,1)=[1.9_dp,1.1_dp,1.1_dp]
  points(:,2)=[2.1_dp,1.1_dp,1.1_dp]
  points(:,3)=[.9_dp,.9_dp,.9_dp]
  points(:,4)=[3.9_dp,1.1_dp,1.1_dp]
  np=2
  if (rank==0) np=0
  if (rank==1) np=4
  positions=points(:,:np)
  snapshot=full;positions=.5_dp
  do node=1,3
    full%fields(node)%values=real(10*node,dp)
  end do
  call distribute()
  do level=0,1
    call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok,level)
    call require(ok.and.all(abs(reference-real(10*(level+1),dp))<1e-13_dp),'analytic serial level view')
    call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
      samples,ok,maximum_level=level)
    call require(ok.and.all(samples==reference),'exact MPI level view')
  end do
  saved=samples
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_level=merge(0,1,rank==0))
  call require(.not.ok.and.all(samples==saved),'level disagreement rollback')
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_level=3)
  call require(.not.ok.and.all(samples==saved),'invalid level rollback')
  if(rank==owners(3)) local%fields(3)%values=ieee_value(0._dp,ieee_quiet_nan)
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_level=0)
  call require(ok.and.all(abs(samples-10._dp)<1e-13_dp),'unrequested deeper payload ignored')
  full=snapshot;positions=points(:,:np);call distribute()
  call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok)
  call require(ok,'serial reference')
  call sample()
  call require(ok.and.all(samples==reference),'exact distributed composite sample')
  saved=samples
  ! Covered parent payload is not part of any composite interpolation request.
  if (rank==owners(1)) local%fields(1)%values(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call sample()
  call require(ok.and.all(samples==reference),'unrequested covered payload ignored')
  call distribute()
  if (rank==owners(2)) local%fields(2)%values=ieee_value(0._dp,ieee_quiet_nan)
  call sample()
  call require(.not.ok.and.all(samples==saved),'requested NaN collective rollback')
  call distribute()
  bounds=zero
  if (rank==0) bounds(1)=1._dp
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,bounds,spacing,periodic,samples,ok)
  call require(.not.ok.and.all(samples==saved),'physical bounds disagreement rollback')
  wrong=owners
  if (rank==0) wrong(1)=1
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,wrong,local,2,positions,zero,spacing,periodic,samples,ok)
  call require(.not.ok.and.all(samples==saved),'owner disagreement rollback')
  if (rank==owners(2)) deallocate(local%fields(2)%values)
  call sample()
  call require(.not.ok.and.all(samples==saved),'missing owned field rollback')
  call distribute()
  good_positions=positions
  if (rank==nranks-1) positions(1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call sample()
  call require(.not.ok.and.all(samples==saved),'bad local position collective rollback')
  positions=good_positions
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_requests=merge(0,1000000,rank==1))
  call require(.not.ok.and.all(samples==saved),'local request capacity rollback')
  owners=modulo(owners+1,nranks);call distribute()
  call sample()
  call require(ok.and.all(samples==reference),'ownership changes preserve samples')
  deallocate(positions);allocate(positions(3,0))
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic, &
    samples,ok,maximum_requests=0)
  call require(ok.and.size(samples,1)==2.and.size(samples,2)==0,'globally empty query')
  call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok,0)
  call require(ok.and.size(reference,2)==0,'empty serial bounded query')
  call sample_amr_tree_spray_field(tree,full,positions,zero,spacing,periodic,reference,ok,3)
  call require(.not.ok.and.size(reference,2)==0,'empty serial invalid level rollback')
  call check_particle_halo()
  if (rank==0) print *, 'PASS: owner-local AMR spray field sampling'
  call MPI_Finalize(ierr)
contains
  subroutine check_particle_halo()
    type(amr_tree_state_3d) :: halo,vf,packed,centroids,deposited,retained_source,cell_source
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),query_nodes(:)
    real(dp), allocatable :: query(:,:),values(:,:),retained(:,:),feedback(:,:)
    real(dp) :: inventory(2),paired(2),expected_pair(2),cell_volume
    real(dp) :: x(3),h(3),f,phi(0:16,0:16,0:16)
    integer :: p,a,b,c,n(3),count,q
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    call require(ok,'halo layout')
    allocate(halo%fields(3),vf%fields(3),packed%fields(3),centroids%fields(3))
    do p=1,3
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape;h=4._dp/real(domain(:,tree%nodes(p)%level+1),dp)
      allocate(halo%fields(p)%values(2,-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
      allocate(centroids%fields(p)%values(3,-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
      allocate(vf%fields(p)%values(1,-1:n(1)+2,-1:n(2)+2,-1:n(3)+2));vf%fields(p)%values=.5_dp
      do c=-1,n(3)+2
      do b=-1,n(2)+2
      do a=-1,n(1)+2
        x=(real(origin(:,p),dp)+real([a,b,c],dp)-.5_dp)*h
        centroids%fields(p)%values(:,a,b,c)=x
        f=2*x(1)+3*x(2)-x(3)+7*.25_dp
        halo%fields(p)%values(:,a,b,c)=[f,-2*f]
      end do
      end do
      end do
      allocate(packed%fields(p)%values(2,n(1)+4,n(2)+4,n(3)+4))
      packed%fields(p)%values=halo%fields(p)%values
    end do
    count=merge(2,0,rank==owners(2));allocate(query_nodes(count),query(3,count));query_nodes=2
    if(count>0) then
      query(:,1)=[2.125_dp,1.125_dp,1.125_dp]
      query(:,2)=[3.875_dp,1.125_dp,1.125_dp]
    end if
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,zero,spacing, &
      periodic,values,ok,vf)
    call require(ok,'complete ghost halo sampling')
    do q=1,count
      x=query(:,q)
      if(q==2) x(1)=x(1)-4
      f=2*x(1)+3*x(2)-x(3)+7*.25_dp
      call require(maxval(abs(values(:,q)-[f,-2*f]))<1e-13_dp,'affine ghost and periodic-image interpolation')
    end do
    retained=values
    allocate(feedback(2,count))
    do q=1,count
      feedback(:,q)=[2._dp,-3._dp]
    end do
    call deposit_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,feedback, &
      zero,spacing,periodic,deposited,ok,vf)
    call require(ok,'halo particle source deposition')
    inventory=0;paired=0;expected_pair=sum(feedback*values,dim=2)
    do p=1,3
      if(owners(p)/=rank) cycle
      cell_volume=product(4._dp/real(domain(:,tree%nodes(p)%level+1),dp))*.5_dp
      inventory=inventory+cell_volume*sum(sum(sum(deposited%fields(p)%values,dim=4),dim=3),dim=2)
      paired=paired+cell_volume*sum(sum(sum(deposited%fields(p)%values*halo%fields(p)%values,dim=4),dim=3),dim=2)
    end do
    call MPI_Allreduce(MPI_IN_PLACE,inventory,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,paired,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,expected_pair,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(maxval(abs(inventory-[4._dp,-6._dp]))<1e-13_dp,'halo source extensive inventory')
    call require(maxval(abs(paired-expected_pair))<1e-13_dp,'sampling and deposition use matching weights')
    call deposit_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,feedback, &
      zero,spacing,periodic,cell_source,ok,vf,containing_cell=.true.)
    call require(ok,'containing-cell deposition')
    if(size(query_nodes)>0) then
      call require(all(cell_source%fields(2)%values(:,5,3,3)==[32._dp,-48._dp]), &
        'off-center contribution goes entirely into containing fine ghost cell')
      call require(all(cell_source%fields(2)%values(:,0,3,3)==[32._dp,-48._dp]), &
        'periodic image contribution uses containing ghost cell')
      call require(sum(abs(cell_source%fields(2)%values))==160._dp,'no source on other donors')
      vf%fields(2)%values(1,5,3,3)=0
    end if
    retained_source=cell_source
    call deposit_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,feedback, &
      zero,spacing,periodic,cell_source,ok,vf,containing_cell=.true.)
    call require(.not.ok,'solid containing cell rejects despite neighboring fluid donors')
    do p=1,3
      if(owners(p)/=rank) cycle
      call require(all(cell_source%fields(p)%values==retained_source%fields(p)%values),'cell deposition rollback')
    end do
    if(size(query_nodes)>0) vf%fields(2)%values(1,5,3,3)=.5_dp
    retained_source=deposited
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,packed,2,2,query_nodes,query,zero,spacing, &
      periodic,values,ok)
    call require(ok.and.all(values==retained),'packed and physical halo bounds agree')
    if(count>0) vf%fields(2)%values=0
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,zero,spacing, &
      periodic,values,ok,vf)
    call require(.not.ok.and.all(values==retained),'zero-fluid support rollback')
    call deposit_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,feedback, &
      zero,spacing,periodic,deposited,ok,vf)
    call require(.not.ok,'zero-fluid deposition rejection')
    do p=1,3
      if(owners(p)/=rank) cycle
      call require(all(deposited%fields(p)%values==retained_source%fields(p)%values),'source halo rollback')
    end do

    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,merge(1,2,rank==0), &
      query_nodes,query,zero,spacing,periodic,values,ok)
    call require(.not.ok.and.all(values==retained),'halo width disagreement rollback')
    if(count>0) query(2,1)=10
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,zero,spacing, &
      [.true.,.false.,.true.],values,ok)
    call require(.not.ok.and.all(values==retained),'insufficient halo rollback')
    do a=0,16
      phi(a,:,:)=abs(.25_dp*a-.5_dp)-.1_dp
    end do
    do p=1,3
      if(owners(p)/=rank) cycle
      vf%fields(p)%values=1;n=tree%nodes(p)%shape
      do c=-1,n(3)+2
      do b=-1,n(2)+2
      do a=-1,n(1)+2
        f=merge(10._dp,100._dp,centroids%fields(p)%values(1,a,b,c)<.5_dp)
        halo%fields(p)%values(:,a,b,c)=[f,-2*f]
      end do
      end do
      end do
    end do
    if(count>0) then
      query(:,1)=[.375_dp,1.125_dp,1.125_dp]
      query(:,2)=query(:,1)
    end if
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,zero,spacing, &
      [.false.,.false.,.false.],values,ok,vf,phi,centroids)
    call require(ok,'occluded halo sampling')
    call deposit_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,feedback, &
      zero,spacing,[.false.,.false.,.false.],deposited,ok,vf,phi,centroids)
    call require(ok,'occluded halo deposition')
    inventory=0
    do p=1,3
      if(owners(p)/=rank) cycle
      cell_volume=product(4._dp/real(domain(:,tree%nodes(p)%level+1),dp))
      inventory=inventory+cell_volume*sum(sum(sum(deposited%fields(p)%values,dim=4),dim=3),dim=2)
      n=tree%nodes(p)%shape
      do c=-1,n(3)+2
      do b=-1,n(2)+2
      do a=-1,n(1)+2
        if(centroids%fields(p)%values(1,a,b,c)>.5_dp) &
          call require(all(deposited%fields(p)%values(:,a,b,c)==0),'blocked donors receive no source')
      end do
      end do
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,inventory,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(maxval(abs(inventory-[4._dp,-6._dp]))<1e-13_dp,'visible source inventory retained')

    if(count>0) call require(all(abs(values(1,:)-10)<1e-13_dp).and.all(abs(values(2,:)+20)<1e-13_dp), &
      'wall blocks otherwise fluid donor')
    retained=values
    if(count>0) query(1,1)=.5_dp
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,zero,spacing, &
      [.false.,.false.,.false.],values,ok,vf,phi,centroids)
    call require(.not.ok.and.all(values==retained),'solid particle rollback')
    if(count>0) query(1,1)=.375_dp
    if(rank==0) phi(16,16,16)=phi(16,16,16)+.01_dp
    call sample_mpi_amr_particle_halo(MPI_COMM_WORLD,tree,owners,halo,2,2,query_nodes,query,zero,spacing, &
      [.false.,.false.,.false.],values,ok,vf,phi,centroids)
    call require(.not.ok.and.all(values==retained),'nodal geometry disagreement rollback')
  end subroutine check_particle_halo
  subroutine sample()
    call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,local,2,positions,zero,spacing,periodic,samples,ok)
  end subroutine
  subroutine distribute()
    integer :: q
    if (allocated(local%fields)) deallocate(local%fields)
    allocate(local%fields(3))
    do q=1,3
      if (owners(q)==rank) local%fields(q)%values=full%fields(q)%values
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL rank ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
