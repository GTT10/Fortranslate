program test_les_mpi_amr_sdc
  use mpi_f08
  use les_mpi_amr_sdc_mod, only: advance_mpi_amr_tree_les_sdc_3d
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use les_amr_sdc_mod, only: advance_amr_patch_tree_les_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction) :: reactions(0)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patch(1)
  type(amr_tree_state_3d) :: u,t,saved,saved_t,local,local_t,retained,retained_t
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),q(12),sound,h(3),dt,bc(12),sc(12),correction,theta,diffusion
  real(dp) :: initial(12),final(12),scale(12),maximum
  integer :: i,j,k,geometry,mode,node,rank,nranks,ierr,ratio
  integer, allocatable :: owners(:)
  real(dp) :: mb(12),ms(12),mc,mt,md
  integer(int64), allocatable :: mpi_counts(:)
  character(len=16) :: faces(6),models(4)
  logical :: ok,flags(4)
  logical :: recursive_steps
  character(len=16) :: command
  integer(int64), allocatable :: counts(:)
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call get_command_argument(1,command)
  recursive_steps=command=='recursive'
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  models=[character(len=16)::'smagorinsky','wale','vreman','dynamic']
  do k=1,4
  do j=1,4
  do i=1,4
    q=0;q(1)=.3_dp;q(2)=2*j;q(3)=3*i;q(4)=k;q(5)=101325+10*j;q(12)=1
    call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
    call require(ok,'initial state')
  end do
  end do
  end do
  h=.001_dp;dt=1e-9_dp;faces='periodic';maximum=0
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok);call require(ok,'patch')
  do geometry=2,2
    if(geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patch(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
    end if
    call require(ok,'tree')
    if(allocated(counts)) deallocate(counts)
    allocate(counts(tree%level_count()),mpi_counts(tree%level_count()),owners(size(tree%nodes)))
    owners=mod([(node-1,node=1,size(tree%nodes))],nranks)
    allocate(local%fields(size(tree%nodes)),local_t%fields(size(tree%nodes)))
    call prolong_amr_patch_tree_3d(tree,root,saved,ok);call require(ok,'gas')
    call prolong_amr_patch_tree_3d(tree,rt,saved_t,ok);call require(ok,'temperature')
    call composite_integrals_amr_patch_tree_3d(tree,saved,h,initial,ok);call require(ok,'initial inventory')
    do mode=1,4
      flags=.true.
      if(mode==3) flags=.false.
      u=saved;t=saved_t
      call advance_amr_patch_tree_les_sdc_3d(tree,species,reactions,transport,u,t,h,dt,faces,2,.false.,flags, &
        [0._dp,0._dp,0._dp],1e-6_dp,1e-12_dp,.17_dp,.7_dp,.7_dp,bc,sc,correction,theta,diffusion,ok, &
        model=trim(models(mode)),subcycling=recursive_steps,level_steps=counts)
      call require(ok.and.diffusion>0,'coupled AMR LES SDC')
      call require(counts(1)==1,'root advances once')
      if(geometry==2) call require(counts(2)==merge(2,1,recursive_steps),'fine-level substep count')
      call composite_integrals_amr_patch_tree_3d(tree,u,h,final,ok);call require(ok,'final inventory')
      scale=max(abs(initial),abs(final),1e-20_dp);scale(2:4)=max(scale(2:4),sqrt(initial(1)*initial(5)))
      maximum=max(maximum,maxval(abs(final-initial-bc-sc)/scale))
      call require(all(abs(final-initial-bc-sc)<=3e-12_dp*scale),'refluxed composite LES inventory')
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        local%fields(node)=saved%fields(node);local_t%fields(node)=saved_t%fields(node)
      end do
      call advance_mpi_amr_tree_les_sdc_3d(tree,owners,species,reactions,transport,local,local_t,h,dt,faces, &
        2,.false.,flags,[0._dp,0._dp,0._dp],1e-6_dp,1e-12_dp,.17_dp,.7_dp,.7_dp,mb,ms,mc,mt,md,MPI_COMM_WORLD,ok, &
        model=trim(models(mode)),subcycling=recursive_steps,level_steps=mpi_counts)
      call require(ok,'MPI LES SDC')
      call require(all(mpi_counts==counts),'MPI substep counts')
      call require(all(mb==bc).and.all(ms==sc).and.mc==correction.and.mt==theta.and.md==diffusion,'MPI diagnostics exact')
      do node=1,size(tree%nodes)
        call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owner-local state')
        if(owners(node)/=rank) cycle
        call require(all(local%fields(node)%values==u%fields(node)%values),'serial MPI state exact')
        call require(all(local_t%fields(node)%values==t%fields(node)%values),'serial MPI temperature exact')
      end do
    end do
    retained=local;retained_t=local_t;mb=42;ms=43;md=44;mpi_counts=45
    ratio=2
    if(rank==nranks-1) ratio=-1
    call advance_mpi_amr_tree_les_sdc_3d(tree,owners,species,reactions,transport,local,local_t,h,dt,faces, &
      2,.false.,flags,[0._dp,0._dp,0._dp],1e-6_dp,1e-12_dp,.17_dp,.7_dp,.7_dp,mb,ms,mc,mt,md,MPI_COMM_WORLD,ok, &
      model='dynamic',test_ratio=ratio,subcycling=recursive_steps,level_steps=mpi_counts)
    call require(.not.ok.and.all(mb==42).and.all(ms==43).and.md==44,'idle-rank invalid filter rejects step')
    call require(all(mpi_counts==45),'failed MPI substep count rollback')
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      call require(all(local%fields(node)%values==retained%fields(node)%values).and. &
        all(local_t%fields(node)%values==retained_t%fields(node)%values),'MPI failure state rollback')
    end do
  end do
  if(rank==0) print *, 'PASS: MPI LES SDC four closures, serial parity, reflux; recursive, residual ',recursive_steps,maximum
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
