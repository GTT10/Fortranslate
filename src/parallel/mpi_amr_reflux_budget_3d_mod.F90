module mpi_amr_reflux_budget_3d_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,get_face,set_face
  use amr_tree_reflux_graph_3d_mod, only: build_amr_reflux_graph_3d
  use conservative_reflux_budget_mod, only: spread_reflux_return
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d
  implicit none
  private
  public :: balance_mpi_amr_reflux_species_3d
contains
  ! A reflux correction is limited as a complete conserved vector. The part
  ! returned to finer cells is funded by their connected composite inventory,
  ! with identical transfers recorded on every copy of every affected face.
  ! This runs on completed level states; fluid payload exchange is restricted
  ! to the fine component funding a failed coarse correction.
  subroutine balance_mpi_amr_reflux_species_3d(tree,owners,original,state,nvar,flux,root_spacing, &
      faces,first_species,last_species,minimum_depth,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nvar,first_species,last_species,minimum_depth
    type(amr_tree_state_3d), intent(in) :: original
    type(amr_tree_state_3d), intent(inout) :: state
    type(amr_tree_flux_3d), intent(inout) :: flux(:)
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    type(amr_tree_flux_3d), allocatable :: work(:)
    integer, allocatable :: targets(:,:),requests(:,:),queries(:,:),edges(:,:),edge_side(:),ratios(:)
    integer(int64), allocatable :: lower(:,:),upper(:,:),origin(:,:),domain(:,:)
    real(dp), allocatable :: values(:,:),accepted(:,:),integrals(:,:),volumes(:),seed(:)
    real(dp) :: before(nvar),current(nvar),coarse(nvar),delta(nvar),returned(nvar),theta,volume,spacing(3),lattice(3)
    integer :: rank,status,node,n(3),i,j,k,depth,nt,target(4),t,c,axis,side,attempt,e,u,v,peer,ios
    integer(int64) :: lo(3),hi(3),plane
    logical :: good,ready,negative,any_negative,funded,success
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    good=first_species>=2.and.last_species>=first_species.and.last_species<=nvar
    call agree(good,ready)
    if(.not.ready) return
    negative=.false.
    do node=1,size(tree%nodes)
      if(owners(node)/=rank.or.tree%nodes(node)%level<minimum_depth) cycle
      n=tree%nodes(node)%shape
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(covered(tree,node,[i,j,k])) cycle
            if(bad(state%fields(node)%values(:,i,j,k))) negative=.true.
          enddo
        enddo
      enddo
    enddo
    call MPI_Allreduce(negative,any_negative,1,MPI_LOGICAL,MPI_LOR,comm,status)
    if(status/=MPI_SUCCESS) return
    if(.not.any_negative) then
      ok=.true.;return
    endif
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    call agree(good,ready)
    if(.not.ready) return
    candidate=state;work=flux
    do depth=tree%level_count()-1,minimum_depth,-1
      spacing=root_spacing*real(domain(:,1),dp)/real(domain(:,depth+1),dp)
      volume=product(spacing)
      if(volume<=0.or..not.ieee_is_finite(volume)) return
      do node=1,size(tree%nodes)
        if(tree%nodes(node)%level/=depth) cycle
        n=tree%nodes(node)%shape;nt=0
        if(rank==owners(node)) then
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                if(covered(tree,node,[i,j,k])) cycle
                if(bad(candidate%fields(node)%values(:,i,j,k))) nt=nt+1
              enddo
            enddo
          enddo
        endif
        call MPI_Bcast(nt,1,MPI_INTEGER,owners(node),comm,status)
        if(status/=MPI_SUCCESS.or.nt>huge(1)/4) return
        if(nt==0) cycle
        allocate(targets(4,nt),stat=ios)
        call agree(ios==0,ready)
        if(.not.ready) return
        if(rank==owners(node)) then
          t=0
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                if(covered(tree,node,[i,j,k])) cycle
                if(.not.bad(candidate%fields(node)%values(:,i,j,k))) cycle
                t=t+1;targets(:,t)=[node,i,j,k]
              enddo
            enddo
          enddo
        endif
        call MPI_Bcast(targets,4*nt,MPI_INTEGER,owners(node),comm,status)
        if(status/=MPI_SUCCESS) return
        do t=1,nt
          target=targets(:,t);i=target(2);j=target(3);k=target(4)
          if(rank==owners(node)) then
            before=original%fields(node)%values(:,i,j,k);current=candidate%fields(node)%values(:,i,j,k)
          endif
          call MPI_Bcast(before,nvar,MPI_DOUBLE_PRECISION,owners(node),comm,status)
          if(status/=MPI_SUCCESS) return
          call MPI_Bcast(current,nvar,MPI_DOUBLE_PRECISION,owners(node),comm,status)
          if(status/=MPI_SUCCESS) return
          if(bad(before).or..not.all(ieee_is_finite(current))) return
          delta=current-before;theta=1
          do c=1,nvar
            if(c/=1.and.(c<first_species.or.c>last_species)) cycle
            if(delta(c)>=0) cycle
            theta=min(theta,before(c)/(-delta(c)))
          enddo
          if(theta<0.or.theta>1) return
          ! Recompute, rather than clip a component. Rounding at an exhausted
          ! budget is handled by reducing the common vector coefficient.
          do attempt=1,16
            coarse=before+theta*delta
            if(.not.bad(coarse)) exit
            if(theta==0) return
            theta=max(0._dp,nearest(theta,-1._dp))
          enddo
          if(bad(coarse)) return
          returned=(current-coarse)*volume;funded=.false.
          do axis=1,3
            do side=-1,1,2
              call build_amr_reflux_graph_3d(tree,root_spacing,faces,target,axis,side,requests,volumes, &
                edges,edge_side,seed,lower,upper,lattice,good)
              call agree(good,ready)
              if(.not.ready) cycle
              if(size(requests,2)>huge(1)/nvar.or.size(edges,2)>huge(1)/nvar) return
              if(rank==owners(node)) then
                allocate(queries,source=requests,stat=ios)
              else
                allocate(queries(4,0),stat=ios)
              endif
              call agree(ios==0,ready)
              if(.not.ready) return
              call fetch_mpi_amr_tree_cells_3d(tree,owners,candidate,nvar,queries,values,comm,good)
              deallocate(queries)
              if(.not.good) return
              allocate(accepted(nvar,size(requests,2)),integrals(nvar,size(edges,2)),stat=ios)
              call agree(ios==0,ready)
              if(.not.ready) return
              success=.false.
              if(rank==owners(node)) call spread_reflux_return(values,volumes,edges,seed,returned, &
                first_species,last_species,accepted,integrals,success)
              call MPI_Bcast(success,1,MPI_LOGICAL,owners(node),comm,status)
              if(status/=MPI_SUCCESS) return
              if(.not.success) then
                deallocate(accepted,integrals);cycle
              endif
              call MPI_Bcast(accepted,size(accepted),MPI_DOUBLE_PRECISION,owners(node),comm,status)
              if(status/=MPI_SUCCESS) return
              call MPI_Bcast(integrals,size(integrals),MPI_DOUBLE_PRECISION,owners(node),comm,status)
              if(status/=MPI_SUCCESS) return
              do c=1,size(requests,2)
                peer=requests(1,c)
                if(owners(peer)/=rank) cycle
                candidate%fields(peer)%values(:,requests(2,c),requests(3,c),requests(4,c))=accepted(:,c)
              enddo
              if(rank==owners(node)) candidate%fields(node)%values(:,i,j,k)=coarse
              ! External supply at each seed plus graph's internal transfers.
              seed=seed/maxval(seed);seed=seed/sum(seed)
              c=maxloc(seed,dim=1);seed(c)=0;seed(c)=1-sum(seed)
              do c=1,size(requests,2)
                if(seed(c)==0) cycle
                lo=lower(:,c);hi=upper(:,c);plane=lo(axis)
                if(side<0) plane=hi(axis)
                call transfer_face(lo,hi,axis,plane,real(side,dp)*seed(c)*returned,good)
                if(.not.good) return
              enddo
              do e=1,size(edges,2)
                u=edges(1,e);v=edges(2,e);peer=abs(edge_side(e))
                lo=max(lower(:,u),lower(:,v));hi=min(upper(:,u),upper(:,v))
                plane=lower(peer,u)
                if(edge_side(e)>0) plane=upper(peer,u)
                call transfer_face(lo,hi,peer,plane,sign(1._dp,real(edge_side(e),dp))*integrals(:,e),good)
                if(.not.good) return
              enddo
              deallocate(accepted,integrals);funded=.true.;exit
            enddo
            if(funded) exit
          enddo
          if(.not.funded) return
        enddo
        deallocate(targets)
      enddo
    enddo
    state=candidate;flux=work;ok=.true.
  contains
    logical function bad(value)
      real(dp), intent(in) :: value(:)
      bad=value(1)<=0.or.any(value(first_species:last_species)<0)
    end function bad
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
    subroutine transfer_face(face_lo,face_hi,normal,position,amount,ready)
      integer(int64), intent(in) :: face_lo(3),face_hi(3),position
      integer, intent(in) :: normal
      real(dp), intent(in) :: amount(:)
      logical, intent(out) :: ready
      integer(int64) :: base(3),last(3),w,pl,low_cell(3),high_cell(3),overlap(3),face_base(3)
      integer :: p,lev,b,d,ii,jj,kk,ix(3),wrap
      real(dp) :: region_area,cell_area,weight,next(nvar)
      logical :: valid
      b=modulo(normal,3)+1;d=modulo(normal+1,3)+1
      region_area=real(face_hi(b)-face_lo(b),dp)*lattice(b)*real(face_hi(d)-face_lo(d),dp)*lattice(d)
      call agree(region_area>0.and.ieee_is_finite(region_area),ready)
      if(.not.ready) return
      valid=.true.
      do wrap=0,1
        pl=position
        if(wrap==1) then
          if(faces(2*normal-1)/='periodic') cycle
          if(position==0) then
            pl=domain(normal,size(ratios))
          else if(position==domain(normal,size(ratios))) then
            pl=0
          else
            cycle
          endif
        endif
        do p=1,size(tree%nodes)
          if(owners(p)/=rank.or.tree%nodes(p)%level<minimum_depth) cycle
          lev=tree%nodes(p)%level+1;w=domain(1,size(ratios))/domain(1,lev)
          base=origin(:,p)*w;last=base+int(tree%nodes(p)%shape,int64)*w
          if(pl<base(normal).or.pl>last(normal)) cycle
          if(modulo(pl-base(normal),w)/=0) cycle
          if(face_hi(b)<=base(b).or.face_lo(b)>=last(b).or.face_hi(d)<=base(d).or.face_lo(d)>=last(d)) cycle
          low_cell=(max(face_lo,base)-base)/w+1;high_cell=(min(face_hi,last)-1-base)/w+1
          low_cell(normal)=(pl-base(normal))/w+1;high_cell(normal)=low_cell(normal)
          cell_area=real(w,dp)*lattice(b)*real(w,dp)*lattice(d)
          if(cell_area<=0.or..not.ieee_is_finite(cell_area)) then
            valid=.false.;cycle
          endif
          do kk=int(low_cell(3)),int(high_cell(3))
            do jj=int(low_cell(2)),int(high_cell(2))
              do ii=int(low_cell(1)),int(high_cell(1))
                ix=[ii,jj,kk];face_base=base+int(ix-1,int64)*w
                overlap=min(face_hi,face_base+w)-max(face_lo,face_base)
                weight=(real(overlap(b),dp)*lattice(b)*real(overlap(d),dp)*lattice(d)/region_area)/cell_area
                next=get_face(work(p),normal,ix)+weight*amount
                if(.not.all(ieee_is_finite(next))) valid=.false.
                call set_face(work(p),normal,ix,next)
              enddo
            enddo
          enddo
        enddo
      enddo
      call agree(valid,ready)
    end subroutine transfer_face
  end subroutine balance_mpi_amr_reflux_species_3d
end module mpi_amr_reflux_budget_3d_mod
