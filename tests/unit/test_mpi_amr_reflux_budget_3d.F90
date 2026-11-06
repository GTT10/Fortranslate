program test_mpi_amr_reflux_budget_3d
  use mpi_f08
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,get_face,reflux_amr_patch_tree_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  implicit none
  integer, parameter :: nv=8
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: original,limited,unlimited,saved
  type(amr_tree_state_3d) :: internal,internal_original,internal_saved
  type(amr_tree_flux_3d), allocatable :: flux(:),effective(:),sentinel(:)
  type(amr_tree_flux_3d), allocatable :: auxiliary_flux(:),auxiliary_effective(:),auxiliary_sentinel(:)
  integer :: rank,nranks,status,node,n(3),owners(3),i,j,k,axis,idx(3),upper(3),mode,first,coarse_i
  real(dp) :: gas_largest
  real(dp) :: budget(nv),largest,spacing(3),difference(nv)
  real(dp) :: local_sum(nv),global_sum(nv),fuel,rate(nv),scalar_delta,scalar_sum,scalar_global
  logical :: ok
  character(len=32) :: faces(6)
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call require(nranks>=2,'ranks')
  faces='slip_wall'
  call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(2,2,2,2,2,2,2,2,2,2,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  rate=[.01_dp,.0002_dp,-.0003_dp,.0001_dp,20._dp,.01_dp,0._dp,0._dp]
  do mode=1,6
    coarse_i=1;faces='slip_wall'
    if(mode==3) then
      call initialize_amr_patch_3d(4,4,4,1,1,2,2,2,2,2,patches(1),ok)
      call require(ok,'periodic child')
      call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
      call require(ok,'periodic tree');coarse_i=4;faces(1:2)='periodic'
    else if(mode==4) then
      call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,2,patches(1),ok)
      call require(ok,'first sibling')
      call initialize_amr_patch_3d(4,4,4,3,3,2,2,2,2,2,patches(2),ok)
      call require(ok,'second sibling')
      call initialize_amr_patch_tree_3d([4,4,4],[1,1],patches,tree,ok)
      call require(ok,'cross-parent tree')
    else if(mode==5) then
      call initialize_amr_patch_3d(2,2,2,2,2,2,2,2,2,2,patches(2),ok)
      call require(ok,'restore grandchild')
      call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
      call require(ok,'partial-budget tree')
    endif
    owners=[0,1,mod(2,nranks)]
    if(mode==2) owners=[nranks-1,0,0]
    allocate(original%fields(3),flux(3))
    do node=1,3
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(original%fields(node)%values(nv,n(1),n(2),n(3)), &
        flux(node)%x(nv,n(1)+1,n(2),n(3)),flux(node)%y(nv,n(1),n(2)+1,n(3)), &
        flux(node)%z(nv,n(1),n(2),n(3)+1))
      flux(node)%x=0;flux(node)%y=0;flux(node)%z=0
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            fuel=0
            if((node==2.and.i==2.and.mode/=4).or.node==3) fuel=.1_dp
            if(mode==5.and.node==1.and.all([i,j,k]==[1,2,2])) fuel=.005_dp
            original%fields(node)%values(:,i,j,k)=[1._dp,0._dp,0._dp,0._dp,1000._dp, &
              fuel,.4_dp,.6_dp-fuel]
          enddo
        enddo
      enddo
      if(node==2) then
        do k=1,n(3)
          do j=1,n(2)
            flux(node)%x(:,1,j,k)=rate
          enddo
        enddo
        if(mode==6) then
          do k=1,n(3)
            do i=1,n(1)
              flux(node)%y(:,i,1,k)=.7_dp*rate
            enddo
          enddo
        endif
      endif
    enddo
    unlimited=original;budget=0;largest=0
    call reflux_mpi_amr_tree_3d(tree,owners,unlimited,nv,flux,[1._dp,1._dp,1._dp],faces,budget,largest, &
      MPI_COMM_WORLD,ok)
    call require(ok,'ordinary reflux')
    if(rank==owners(1)) call require(unlimited%fields(1)%values(6,coarse_i,2,2)<0,'fixture exhausts coarse fuel')
    limited=original;budget=0;largest=0
    call reflux_mpi_amr_tree_3d(tree,owners,limited,nv,flux,[1._dp,1._dp,1._dp],faces,budget,largest, &
      MPI_COMM_WORLD,ok,synchronized_flux=effective,first_species=6,last_species=8)
    call require(ok,'mixed-level wide donor reflux')
    call require(all(budget==0),'no external boundary transfer')
    call check_serial(.false.,.false.)
    if(rank==owners(1)) then
      if(mode==5) then
        call require(maxval(abs(limited%fields(1)%values(:,coarse_i,2,2)- &
          original%fields(1)%values(:,coarse_i,2,2)+.5_dp*rate))<2.e-12_dp,'partial coarse vector retained')
      else
        call require(all(limited%fields(1)%values(:,coarse_i,2,2)== &
          original%fields(1)%values(:,coarse_i,2,2)),'zero-budget coarse cell preserved exactly')
      endif
    endif
    local_sum=0
    do node=1,3
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;spacing=1._dp/(2._dp**tree%nodes(node)%level)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(covered(tree,node,[i,j,k])) cycle
            call require(all(limited%fields(node)%values(6:8,i,j,k)>=0),'strict nonnegative species')
            call require(abs(sum(limited%fields(node)%values(6:8,i,j,k))- &
              limited%fields(node)%values(1,i,j,k))<1.e-15_dp,'density sum')
            idx=[i,j,k];difference=0
            do axis=1,3
              upper=idx;upper(axis)=upper(axis)+1
              difference=difference+(get_face(effective(node),axis,idx)-get_face(flux(node),axis,idx)- &
                get_face(effective(node),axis,upper)+get_face(flux(node),axis,upper))/spacing(axis)
            enddo
            call require(maxval(abs(limited%fields(node)%values(:,i,j,k)- &
              original%fields(node)%values(:,i,j,k)-difference))<2.e-12_dp,'independent face divergence audit')
            local_sum=local_sum+(limited%fields(node)%values(:,i,j,k)- &
              unlimited%fields(node)%values(:,i,j,k))*product(spacing)
          enddo
        enddo
      enddo
    enddo
    call MPI_Allreduce(local_sum,global_sum,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(maxval(abs(global_sum))<2.e-12_dp,'same composite inventory as unrestricted vector')
    if(rank==owners(2)) then
      call require(all(limited%fields(2)%values(6,1,:,:)==0),'empty boundary donors preserved')
      if(mode/=4) call require(limited%fields(2)%values(6,2,1,1)<.1_dp,'distant fine donor used')
    endif
    if(rank==owners(3)) call require(limited%fields(3)%values(6,1,1,1)<.1_dp,'finer descendant donor used')
    saved=limited
    gas_largest=largest
    allocate(internal%fields(3),auxiliary_flux(3))
    do node=1,3
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(internal%fields(node)%values(1,n(1),n(2),n(3)), &
        auxiliary_flux(node)%x(1,n(1)+1,n(2),n(3)),auxiliary_flux(node)%y(1,n(1),n(2)+1,n(3)), &
        auxiliary_flux(node)%z(1,n(1),n(2),n(3)+1))
      internal%fields(node)%values=400
      auxiliary_flux(node)%x=0;auxiliary_flux(node)%y=0;auxiliary_flux(node)%z=0
      if(node==2) auxiliary_flux(node)%x(1,1,:,:)=30
      if(node==2.and.mode==6) auxiliary_flux(node)%y(1,:,1,:)=21._dp
    enddo
    internal_original=internal;limited=original;budget=0;largest=0
    call reflux_mpi_amr_tree_3d(tree,owners,limited,nv,flux,[1._dp,1._dp,1._dp],faces,budget,largest, &
      MPI_COMM_WORLD,ok,synchronized_flux=effective,first_species=6,last_species=8, &
      auxiliary_state=internal,auxiliary_flux=auxiliary_flux,auxiliary_synchronized_flux=auxiliary_effective)
    call require(ok,'coupled independent internal energy')
    call require(largest==gas_largest,'held-energy flux does not change conserved maximum correction')
    call check_serial(.true.,.false.)
    scalar_sum=0
    do node=1,3
      if(owners(node)/=rank) cycle
      call require(all(limited%fields(node)%values==saved%fields(node)%values),'packing preserves accepted gas')
      n=tree%nodes(node)%shape;spacing=1._dp/(2._dp**tree%nodes(node)%level)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(covered(tree,node,[i,j,k])) cycle
            scalar_delta=internal%fields(node)%values(1,i,j,k)-400
            call require(abs(scalar_delta-3000*(limited%fields(node)%values(6,i,j,k)- &
              original%fields(node)%values(6,i,j,k)))<2.e-12_dp,'same complete-vector transfer coefficient')
            idx=[i,j,k];difference=0
            do axis=1,3
              upper=idx;upper(axis)=upper(axis)+1
              difference(1)=difference(1)+(sum(get_face(auxiliary_effective(node),axis,idx))- &
                sum(get_face(auxiliary_flux(node),axis,idx))-sum(get_face(auxiliary_effective(node),axis,upper))+ &
                sum(get_face(auxiliary_flux(node),axis,upper)))/spacing(axis)
            enddo
            call require(abs(scalar_delta-difference(1))<2.e-12_dp,'independent internal-energy face divergence')
            scalar_sum=scalar_sum+scalar_delta*product(spacing)
          enddo
        enddo
      enddo
    enddo
    call MPI_Allreduce(scalar_sum,scalar_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(abs(scalar_global+merge(51._dp,30._dp,mode==6))<2.e-12_dp,'independent internal energy composite inventory')
    limited=original
    do node=1,3
      if(owners(node)/=rank) cycle
      limited%fields(node)%values(8,:,:,:)=limited%fields(node)%values(8,:,:,:)+ &
        limited%fields(node)%values(6,:,:,:)
      limited%fields(node)%values(6,:,:,:)=0
    enddo
    saved=limited;sentinel=effective;budget=17;largest=19
    call reflux_mpi_amr_tree_3d(tree,owners,limited,nv,flux,[1._dp,1._dp,1._dp],faces,budget,largest, &
      MPI_COMM_WORLD,ok,synchronized_flux=effective,first_species=6,last_species=8)
    call require(.not.ok,'insufficient composite donor budget rejects')
    call check_serial(.false.,.true.)
    call unchanged()
    first=6
    if(rank==0) first=7
    call reflux_mpi_amr_tree_3d(tree,owners,limited,nv,flux,[1._dp,1._dp,1._dp],faces,budget,largest, &
      MPI_COMM_WORLD,ok,synchronized_flux=effective,first_species=first,last_species=8)
    call require(.not.ok,'rank-disagreeing species interval rejects')
    call unchanged()
    internal_saved=internal;auxiliary_sentinel=auxiliary_effective
    call reflux_mpi_amr_tree_3d(tree,owners,limited,nv,flux,[1._dp,1._dp,1._dp],faces,budget,largest, &
      MPI_COMM_WORLD,ok,synchronized_flux=effective,first_species=6,last_species=8, &
      auxiliary_state=internal,auxiliary_flux=auxiliary_flux,auxiliary_synchronized_flux=auxiliary_effective)
    call require(.not.ok,'coupled insufficient budget rejects atomically')
    call check_serial(.true.,.true.)
    call unchanged()
    do node=1,3
      if(owners(node)/=rank) cycle
      call require(all(internal%fields(node)%values==internal_saved%fields(node)%values),'internal energy rollback')
      call require(all(auxiliary_effective(node)%x==auxiliary_sentinel(node)%x).and. &
        all(auxiliary_effective(node)%y==auxiliary_sentinel(node)%y).and. &
        all(auxiliary_effective(node)%z==auxiliary_sentinel(node)%z),'internal energy flux rollback')
    enddo
    deallocate(original%fields,flux,effective,sentinel,internal%fields,auxiliary_flux,auxiliary_effective,auxiliary_sentinel)
  enddo
  if(rank==0) print *,'mixed-level conservative reflux budget PASS'
  call MPI_Finalize(status)
contains
  subroutine check_serial(coupled,failure)
    logical, intent(in) :: coupled,failure
    type(amr_tree_state_3d) :: ss,sa,ss_before,sa_before
    type(amr_tree_flux_3d), allocatable :: sf(:),saf(:),se(:),sae(:),se_before(:),sae_before(:)
    real(dp) :: sb(nv),sm
    logical :: success
    integer :: p
    if(failure) then
      call gather_state(saved,ss,nv)
    else
      call gather_state(original,ss,nv)
    endif
    call gather_flux(flux,sf,nv)
    call gather_flux(effective,se,nv)
    ss_before=ss;se_before=se;sb=17;sm=19
    if(coupled) then
      if(failure) then
        call gather_state(internal_saved,sa,1)
      else
        call gather_state(internal_original,sa,1)
      endif
      call gather_flux(auxiliary_flux,saf,1)
      call gather_flux(auxiliary_effective,sae,1)
      sa_before=sa;sae_before=sae
      call reflux_amr_patch_tree_3d(tree,ss,sf,[1._dp,1._dp,1._dp],faces,sb,sm,success, &
        first_species=6,last_species=8,synchronized_flux=se,auxiliary_state=sa,auxiliary_flux=saf, &
        auxiliary_synchronized_flux=sae)
    else
      call reflux_amr_patch_tree_3d(tree,ss,sf,[1._dp,1._dp,1._dp],faces,sb,sm,success, &
        first_species=6,last_species=8,synchronized_flux=se)
    endif
    call require(success.eqv.(.not.failure),'serial/MPI reflux acceptance parity')
    if(failure) then
      call require(all(sb==17).and.sm==19,'serial diagnostic rollback')
    else
      call require(all(sb==budget).and.sm==largest,'serial/MPI diagnostic parity')
    endif
    do p=1,3
      if(failure) then
        call require(all(ss%fields(p)%values==ss_before%fields(p)%values),'serial state rollback')
        call equal_flux(se(p),se_before(p),'serial effective flux rollback')
        if(coupled) then
          call require(all(sa%fields(p)%values==sa_before%fields(p)%values),'serial auxiliary state rollback')
          call equal_flux(sae(p),sae_before(p),'serial auxiliary flux rollback')
        endif
      else if(owners(p)==rank) then
        call require(maxval(abs(ss%fields(p)%values-limited%fields(p)%values))<2.e-12_dp, &
          'serial/MPI conservative repaired state')
        call equal_flux(se(p),effective(p),'serial/MPI effective flux')
        if(coupled) then
          call require(maxval(abs(sa%fields(p)%values-internal%fields(p)%values))<2.e-12_dp, &
            'serial/MPI held-energy transfer')
          call equal_flux(sae(p),auxiliary_effective(p),'serial/MPI auxiliary effective flux')
        endif
      endif
    enddo
    ! Invalid interval and mismatched auxiliary presence must be atomic as well.
    ss_before=ss;se_before=se;sb=17;sm=19
    call reflux_amr_patch_tree_3d(tree,ss,sf,[1._dp,1._dp,1._dp],faces,sb,sm,success, &
      first_species=6,synchronized_flux=se)
    call require(.not.success.and.all(sb==17).and.sm==19,'serial missing species bound refuses')
    do p=1,3
      call require(all(ss%fields(p)%values==ss_before%fields(p)%values),'invalid interval state rollback')
      call equal_flux(se(p),se_before(p),'invalid interval flux rollback')
    enddo
    if(coupled) then
      call reflux_amr_patch_tree_3d(tree,ss,sf,[1._dp,1._dp,1._dp],faces,sb,sm,success,auxiliary_state=sa)
      call require(.not.success.and.all(sb==17).and.sm==19,'serial incomplete auxiliary arguments refuse')
    endif
  end subroutine check_serial
  subroutine gather_state(local,full,nc)
    type(amr_tree_state_3d), intent(in) :: local
    type(amr_tree_state_3d), intent(out) :: full
    integer, intent(in) :: nc
    integer :: p,nn(3),ierr
    full=local
    do p=1,3
      nn=tree%nodes(p)%shape
      if(owners(p)/=rank) allocate(full%fields(p)%values(nc,nn(1),nn(2),nn(3)))
      call MPI_Bcast(full%fields(p)%values,size(full%fields(p)%values), &
        MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,ierr)
      call require(ierr==MPI_SUCCESS,'gather serial state')
    enddo
  end subroutine gather_state
  subroutine gather_flux(local,full,nc)
    type(amr_tree_flux_3d), intent(in) :: local(:)
    type(amr_tree_flux_3d), allocatable, intent(out) :: full(:)
    integer, intent(in) :: nc
    integer :: p,nn(3),ierr
    full=local
    do p=1,3
      nn=tree%nodes(p)%shape
      if(owners(p)/=rank) allocate(full(p)%x(nc,nn(1)+1,nn(2),nn(3)), &
        full(p)%y(nc,nn(1),nn(2)+1,nn(3)),full(p)%z(nc,nn(1),nn(2),nn(3)+1))
      call MPI_Bcast(full(p)%x,size(full(p)%x),MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,ierr)
      call require(ierr==MPI_SUCCESS,'gather x flux')
      call MPI_Bcast(full(p)%y,size(full(p)%y),MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,ierr)
      call require(ierr==MPI_SUCCESS,'gather y flux')
      call MPI_Bcast(full(p)%z,size(full(p)%z),MPI_DOUBLE_PRECISION,owners(p),MPI_COMM_WORLD,ierr)
      call require(ierr==MPI_SUCCESS,'gather z flux')
    enddo
  end subroutine gather_flux
  subroutine equal_flux(a,b,message)
    type(amr_tree_flux_3d), intent(in) :: a,b
    character(len=*), intent(in) :: message
    call require(all(a%x==b%x).and.all(a%y==b%y).and.all(a%z==b%z),message)
  end subroutine equal_flux
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: ierr
    if(condition) return
    write(*,*) 'FAIL rank',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine require
  subroutine unchanged()
    integer :: p
    call require(all(budget==17).and.largest==19,'diagnostic rollback')
    do p=1,3
      if(owners(p)/=rank) cycle
      call require(all(limited%fields(p)%values==saved%fields(p)%values),'state rollback')
      call require(all(effective(p)%x==sentinel(p)%x).and.all(effective(p)%y==sentinel(p)%y).and. &
        all(effective(p)%z==sentinel(p)%z),'synchronized flux rollback')
    enddo
  end subroutine unchanged
end program test_mpi_amr_reflux_budget_3d
