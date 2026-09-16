module amr_patch_tree_ghost_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d,valid_wall_velocity_3d
  implicit none
  private
  public :: fill_amr_patch_tree_ghosts_3d
contains
  subroutine fill_amr_patch_tree_ghosts_3d(tree,state,node,faces,ghost,ok,viscous_wall,inflow_state,wall_velocity, &
      coarse_future,coarse_fractions)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: state
    integer, intent(in) :: node
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: viscous_wall
    real(dp), optional, intent(in) :: inflow_state(:,:),wall_velocity(:,:)
    type(amr_tree_state_3d), optional, intent(in) :: coarse_future
    real(dp), optional, intent(in) :: coarse_fractions(:)
    real(dp), allocatable :: candidate(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: coordinate(3),mapped(3),point(3),extent(3)
    integer :: i,j,k,axis,face,inlet,level,source_level,source,parent,other,nvar,local(3),dims(3)
    real(dp) :: shift(3),old_momentum(3),rho,fraction
    logical :: reflect(3),viscous,found
    ok=.false.
    if (.not.state%is_valid(tree)) return
    if (node<1.or.node>size(tree%nodes)) return
    if (.not.valid_reactive_boundaries_3d(faces)) return
    nvar=size(state%fields(1)%values,1)
    if (nvar<5.or.any(tree%nodes(node)%shape>huge(1)-4)) return
    if (present(coarse_future).neqv.present(coarse_fractions)) return
    if (present(coarse_future)) then
      if (size(coarse_fractions)/=tree%nodes(node)%level) return
      if (.not.all(ieee_is_finite(coarse_fractions))) return
      if (any(coarse_fractions<0).or.any(coarse_fractions>1)) return
      if (.not.allocated(coarse_future%fields)) return
      if (size(coarse_future%fields)/=size(tree%nodes)) return
      ! Only coarser levels need a future endpoint; peers/finer fields may be absent.
      do other=1,size(tree%nodes)
        if (tree%nodes(other)%level>=tree%nodes(node)%level) cycle
        if (.not.allocated(coarse_future%fields(other)%values)) return
        if (any(shape(coarse_future%fields(other)%values)/=[nvar,tree%nodes(other)%shape])) return
        if (.not.all(ieee_is_finite(coarse_future%fields(other)%values))) return
      end do
    end if
    if (present(wall_velocity)) then
      if (.not.valid_wall_velocity_3d(faces,wall_velocity)) return
    end if
    if (any(faces=='inflow')) then
      if (.not.present(inflow_state)) return
      if (any(shape(inflow_state)/=[nvar,6])) return
      if (.not.all(ieee_is_finite(inflow_state))) return
    end if
    viscous=.false.
    if (present(viscous_wall)) viscous=viscous_wall
    allocate(origin(3,size(tree%nodes)),domain(3,tree%level_count()),ratios(tree%level_count()))
    origin=0;ratios=1;domain(:,1)=int(tree%nodes(1)%shape,int64)
    do other=2,size(tree%nodes)
      ratios(tree%nodes(other)%level+1)=tree%nodes(other)%patch%refinement_ratio
    end do
    do level=2,size(ratios)
      if (any(domain(:,level-1)>(huge(1_int64)-2)/int(ratios(level),int64))) return
      domain(:,level)=domain(:,level-1)*int(ratios(level),int64)
    end do
    do other=2,size(tree%nodes)
      parent=tree%nodes(other)%parent
      origin(:,other)=(origin(:,parent)+int([tree%nodes(other)%patch%coarse_i_lower-1, &
        tree%nodes(other)%patch%coarse_j_lower-1,tree%nodes(other)%patch%coarse_k_lower-1],int64))* &
        int(tree%nodes(other)%patch%refinement_ratio,int64)
    end do
    level=tree%nodes(node)%level+1;extent=domain(:,level);dims=tree%nodes(node)%shape+4
    allocate(candidate(nvar,dims(1),dims(2),dims(3)))
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          coordinate=origin(:,node)+int([i,j,k]-3,int64)
          mapped=coordinate;reflect=.false.;shift=0;inlet=0
          do axis=1,3
            if (coordinate(axis)>=0.and.coordinate(axis)<extent(axis)) cycle
            face=2*axis
            if (coordinate(axis)<0) face=face-1
            select case(trim(faces(face)))
            case('periodic')
              mapped(axis)=modulo(coordinate(axis),extent(axis))
            case('inflow')
              inlet=face
              mapped(axis)=max(0_int64,min(extent(axis)-1,coordinate(axis)))
            case('outflow')
              mapped(axis)=max(0_int64,min(extent(axis)-1,coordinate(axis)))
            case('slip_wall','no_slip_wall','isothermal_wall')
              mapped(axis)=-coordinate(axis)-1
              if (coordinate(axis)>=extent(axis)) &
                mapped(axis)=extent(axis)-1-(coordinate(axis)-extent(axis))
              if (faces(face)/='slip_wall'.and.viscous) then
                reflect=.not.reflect;shift=-shift
                if (present(wall_velocity)) shift=shift+2.0_dp*wall_velocity(:,face)
              else
                reflect(axis)=.not.reflect(axis);shift(axis)=-shift(axis)
              end if
            end select
          end do
          ! Search all peers at this resolution, including children of other parents.
          ! If no peer covers this cell, coarsen its global index until a donor exists.
          point=mapped;source_level=level;found=.false.
          do while(source_level>=1)
            do source=1,size(tree%nodes)
              if (tree%nodes(source)%level+1/=source_level) cycle
              if (any(point<origin(:,source))) cycle
              if (any(point-origin(:,source)>=int(tree%nodes(source)%shape,int64))) cycle
              local=int(point-origin(:,source))+1
              candidate(:,i,j,k)=state%fields(source)%values(:,local(1),local(2),local(3))
              if (present(coarse_future).and.source_level<level) then
                fraction=coarse_fractions(source_level)
                if (fraction==1.0_dp) then
                  candidate(:,i,j,k)=coarse_future%fields(source)%values(:,local(1),local(2),local(3))
                else if (fraction>0.0_dp) then
                  candidate(:,i,j,k)=(1.0_dp-fraction)*candidate(:,i,j,k)+ &
                    fraction*coarse_future%fields(source)%values(:,local(1),local(2),local(3))
                end if
              end if
              found=.true.;exit
            end do
            if (found) exit
            point=point/int(ratios(source_level),int64)
            source_level=source_level-1
          end do
          if (.not.found) return
          if (inlet>0) candidate(:,i,j,k)=inflow_state(:,inlet)
          old_momentum=candidate(2:4,i,j,k)
          do axis=1,3
            if (reflect(axis)) candidate(axis+1,i,j,k)=-candidate(axis+1,i,j,k)
          end do
          if (any(abs(shift)>0.0_dp)) then
            rho=candidate(1,i,j,k)
            if (rho<=0.0_dp) return
            candidate(2:4,i,j,k)=candidate(2:4,i,j,k)+rho*shift
            candidate(5,i,j,k)=candidate(5,i,j,k)+ &
              (sum(candidate(2:4,i,j,k)**2)-sum(old_momentum**2))/(2.0_dp*rho)
          end if
        end do
      end do
    end do
    if (.not.all(ieee_is_finite(candidate))) return
    call move_alloc(candidate,ghost);ok=.true.
  end subroutine fill_amr_patch_tree_ghosts_3d
end module amr_patch_tree_ghost_3d_mod
