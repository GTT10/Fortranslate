module amr_patch_tree_ghost_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use slope_limiter_mod, only: limited_slope
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d,valid_wall_velocity_3d, &
    fill_pressure_outlet_ghosts_3d
  implicit none
  private
  public :: fill_amr_patch_tree_ghosts_3d,sample_amr_patch_tree_ghosts_3d
  public :: apply_amr_tree_pressure_ghosts_3d
  public :: amr_tree_cell_source_3d
  type, abstract :: amr_tree_cell_source_3d
  contains
    procedure(tree_cell_reader), deferred :: read
  end type
  type, extends(amr_tree_cell_source_3d) :: serial_tree_source
    type(amr_tree_state_3d), pointer :: state=>null(),future=>null()
  contains
    procedure :: read=>read_serial_cell
  end type
  abstract interface
    subroutine tree_cell_reader(self,node,index,future,value,ok)
      import dp,amr_tree_cell_source_3d
      class(amr_tree_cell_source_3d), intent(inout) :: self
      integer, intent(in) :: node,index(3)
      logical, intent(in) :: future
      real(dp), intent(out) :: value(:)
      logical, intent(out) :: ok
    end subroutine
  end interface
contains
  ! Packed ghosts hold conserved components followed by one temperature guess.
  ! Call on a private buffer after ordinary tree ghost sampling.
  subroutine apply_amr_tree_pressure_ghosts_3d(tree,species,node,faces,pressure,ghost,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: node
    character(len=*), intent(in) :: faces(6)
    real(dp), optional, intent(in) :: pressure(6)
    real(dp), intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: n(3),nv,axis,level
    real(dp) :: local_pressure(6)
    ok=.true.
    if (.not.present(pressure)) return
    ok=.false.
    if (.not.all(ieee_is_finite(pressure))) return
    if (any(pressure<0).or.any(pressure>1e15_dp)) return
    if (any(pressure>0.and.faces/='outflow')) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) return
    ok=.false.
    if (node<1.or.node>size(tree%nodes)) return
    n=tree%nodes(node)%shape;nv=5+size(species)
    if (any(shape(ghost)/=[nv+1,n+4])) return
    level=tree%nodes(node)%level+1;local_pressure=0
    do axis=1,3
      if (origin(axis,node)==0) local_pressure(2*axis-1)=pressure(2*axis-1)
      if (origin(axis,node)+int(n(axis),int64)==domain(axis,level)) local_pressure(2*axis)=pressure(2*axis)
    end do
    call fill_pressure_outlet_ghosts_3d(species,n,faces,local_pressure,ghost(:nv,:,:,:),ghost(nv+1,:,:,:),ok)
  end subroutine

  subroutine fill_amr_patch_tree_ghosts_3d(tree,state,node,faces,ghost,ok,viscous_wall,inflow_state,wall_velocity, &
      coarse_future,coarse_fractions,coarse_limiter,species_count,ghost_width)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), target, intent(in) :: state
    integer, intent(in) :: node
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: viscous_wall
    real(dp), optional, intent(in) :: inflow_state(:,:),wall_velocity(:,:)
    type(amr_tree_state_3d), optional, target, intent(in) :: coarse_future
    real(dp), optional, intent(in) :: coarse_fractions(:)
    character(len=*), optional, intent(in) :: coarse_limiter
    integer, optional, intent(in) :: species_count
    integer, optional, intent(in) :: ghost_width
    type(serial_tree_source) :: reader
    integer :: nvar,other,nspecies
    ok=.false.
    if (.not.state%is_valid(tree)) return
    if (node<1.or.node>size(tree%nodes)) return
    if (.not.valid_reactive_boundaries_3d(faces)) return
    nvar=size(state%fields(1)%values,1)
    if (nvar<5.or.any(tree%nodes(node)%shape>huge(1)-4)) return
    nspecies=0
    if (present(species_count)) nspecies=species_count
    if (nspecies<0.or.nspecies>nvar-5) return
    if (present(coarse_limiter)) then
      if (coarse_limiter/='mc'.and.coarse_limiter/='minmod') return
    end if
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
    reader%state=>state
    if (present(coarse_future)) reader%future=>coarse_future
    call sample_amr_patch_tree_ghosts_3d(tree,nvar,node,faces,reader,ghost,ok,viscous_wall, &
      inflow_state,wall_velocity,coarse_fractions,coarse_limiter,species_count,ghost_width=ghost_width)
  end subroutine fill_amr_patch_tree_ghosts_3d

  subroutine read_serial_cell(self,node,index,future,value,ok)
    class(serial_tree_source), intent(inout) :: self
    integer, intent(in) :: node,index(3)
    logical, intent(in) :: future
    real(dp), intent(out) :: value(:)
    logical, intent(out) :: ok
    if (future) then
      value=self%future%fields(node)%values(:,index(1),index(2),index(3))
    else
      value=self%state%fields(node)%values(:,index(1),index(2),index(3))
    end if
    ok=.true.
  end subroutine

  subroutine sample_amr_patch_tree_ghosts_3d(tree,nvar,node,faces,reader,ghost,ok,viscous_wall,inflow_state,wall_velocity, &
      coarse_fractions,coarse_limiter,species_count,planning,validate_only,ghost_width)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nvar
    class(amr_tree_cell_source_3d), intent(inout) :: reader
    logical, optional, intent(in) :: planning,validate_only
    integer, optional, intent(in) :: ghost_width
    integer, intent(in) :: node
    character(len=*), intent(in) :: faces(6)
    real(dp), allocatable, intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: viscous_wall
    real(dp), optional, intent(in) :: inflow_state(:,:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: coarse_fractions(:)
    character(len=*), optional, intent(in) :: coarse_limiter
    integer, optional, intent(in) :: species_count
    real(dp), allocatable :: candidate(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: coordinate(3),mapped(3),point(3),extent(3)
    integer :: i,j,k,axis,face,inlet,level,source_level,source,local(3),dims(3),ios
    real(dp) :: endpoint(nvar)
    real(dp) :: shift(3),old_momentum(3),rho,fraction
    logical :: reflect(3),viscous,found,valid,plan
    integer :: nspecies,width
    ok=.false.
    plan=.false.
    if (present(planning)) plan=planning
    width=2
    if (present(ghost_width)) width=ghost_width
    if(width<1.or.width>128) return
    if (.not.tree%is_valid()) return
    if (node<1.or.node>size(tree%nodes)) return
    if (.not.valid_reactive_boundaries_3d(faces)) return
    if (nvar<5.or.any(tree%nodes(node)%shape>huge(1)-2*width)) return
    nspecies=0
    if (present(species_count)) nspecies=species_count
    if (nspecies<0.or.nspecies>nvar-5) return
    if (present(coarse_limiter)) then
      if (coarse_limiter/='mc'.and.coarse_limiter/='minmod') return
    end if
    if (present(coarse_fractions)) then
      if (size(coarse_fractions)/=tree%nodes(node)%level) return
      if (.not.all(ieee_is_finite(coarse_fractions))) return
      if (any(coarse_fractions<0).or.any(coarse_fractions>1)) return
    end if
    if (present(wall_velocity)) then
      if (.not.valid_wall_velocity_3d(faces,wall_velocity)) return
    end if
    if (any(faces=='inflow')) then
      if (.not.present(inflow_state)) return
      if (any(shape(inflow_state)/=[nvar,6])) return
      if (.not.all(ieee_is_finite(inflow_state))) return
    end if
    if (present(validate_only)) then
      if (validate_only) then
        ok=.true.;return
      end if
    end if
    viscous=.false.
    if (present(viscous_wall)) viscous=viscous_wall
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) return
    ok=.false.
    level=tree%nodes(node)%level+1;extent=domain(:,level);dims=tree%nodes(node)%shape+2*width
    ! Single wall reflection requires no more than one full domain of donor cells.
    do axis=1,3
      if(any(faces(2*axis-1:2*axis)=='slip_wall').or.any(faces(2*axis-1:2*axis)=='no_slip_wall').or. &
        any(faces(2*axis-1:2*axis)=='isothermal_wall')) then
        if(int(width,int64)>extent(axis)) return
      end if
    end do
    allocate(candidate(nvar,dims(1),dims(2),dims(3)),stat=ios)
    if (ios/=0) return
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          coordinate=origin(:,node)+int([i,j,k]-width-1,int64)
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
              call reader%read(source,local,.false.,candidate(:,i,j,k),valid)
              if (.not.valid) return
              if (present(coarse_fractions).and.source_level<level) then
                fraction=coarse_fractions(source_level)
                if (fraction==1.0_dp) then
                  call reader%read(source,local,.true.,candidate(:,i,j,k),valid)
                  if (.not.valid) return
                else if (fraction>0.0_dp) then
                  call reader%read(source,local,.true.,endpoint,valid)
                  if (.not.valid) return
                  candidate(:,i,j,k)=(1.0_dp-fraction)*candidate(:,i,j,k)+fraction*endpoint
                end if
              end if
              found=.true.;exit
            end do
            if (found) exit
            point=point/int(ratios(source_level),int64)
            source_level=source_level-1
          end do
          if (.not.found) return
          if (present(coarse_limiter).and.source_level<level.and.inlet==0) then
            call interpolate_donor(point,source_level,mapped,candidate(:,i,j,k),valid)
            if (.not.valid) return
          end if
          if (plan) cycle
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
  contains
    subroutine interpolate_donor(center_index,donor_level,fine_index,value,success)
      integer(int64), intent(in) :: center_index(3),fine_index(3)
      integer, intent(in) :: donor_level
      real(dp), intent(inout) :: value(:)
      logical, intent(out) :: success
      real(dp) :: center(nvar),neighbor(nvar),minus(nvar,3),plus(nvar,3),slopes(nvar,3)
      real(dp) :: minimum(nvar),maximum(nvar),offset(3),radius(nvar),theta,half_width(3)
      integer(int64) :: query(3),factor(3)
      integer :: di,dj,dk,direction,component,delta(3)
      logical :: good
      success=.false.;center=value;minimum=center;maximum=center
      factor=domain(:,level)/domain(:,donor_level)
      offset=(real(modulo(fine_index,factor),dp)+.5_dp)/real(factor,dp)-.5_dp
      half_width=.5_dp*(1.0_dp-1.0_dp/real(factor,dp))
      ! Corner bounds retain affine fields in multiple dimensions. Sampling
      ! only face extrema would unnecessarily flatten a diagonal linear field.
      do dk=-1,1
        do dj=-1,1
          do di=-1,1
            delta=[di,dj,dk]
            if (all(delta==0)) cycle
            query=center_index+int(delta,int64)
            call sample_neighbor(query,donor_level,neighbor,good)
            if (.not.good) return
            minimum=min(minimum,neighbor);maximum=max(maximum,neighbor)
            do direction=1,3
              if (sum(abs(delta))/=1.or.delta(direction)==0) cycle
              if (delta(direction)<0) minus(:,direction)=neighbor
              if (delta(direction)>0) plus(:,direction)=neighbor
            end do
          end do
        end do
      end do
      if (plan) then
        success=.true.;return
      end if
      do direction=1,3
        do component=1,nvar
          call limited_slope(center(component)-minus(component,direction), &
            plus(component,direction)-center(component),coarse_limiter,slopes(component,direction),good)
          if (.not.good) return
        end do
        ! Preserve sum(rhoY)=rho with one dependent species slope.
        if (nspecies>0) slopes(5+nspecies,direction)=slopes(1,direction)-sum(slopes(6:4+nspecies,direction))
      end do
      radius=matmul(abs(slopes),half_width);theta=1
      do component=1,nvar
        if (radius(component)<=0) cycle
        theta=min(theta,max(0.0_dp,min(maximum(component)-center(component), &
          center(component)-minimum(component)))/radius(component))
      end do
      ! One scale for all components and all children retains affine closure
      ! and a zero-mean correction over a complete donor-cell subdivision.
      value=center+theta*matmul(slopes,offset)
      success=all(ieee_is_finite(value))
    end subroutine
    subroutine sample_neighbor(query,requested_level,value,success)
      integer(int64), intent(in) :: query(3)
      integer, intent(in) :: requested_level
      real(dp), intent(out) :: value(:)
      logical, intent(out) :: success
      integer(int64) :: position(3),bounds(3)
      integer :: direction,boundary,reservoir,depth,peer,index(3)
      real(dp) :: future_value(nvar)
      real(dp) :: velocity_shift(3),momentum(3),blend,density
      logical :: flip(3)
      position=query;bounds=domain(:,requested_level);flip=.false.;velocity_shift=0;reservoir=0
      success=.false.
      do direction=1,3
        if (query(direction)>=0.and.query(direction)<bounds(direction)) cycle
        boundary=2*direction
        if (query(direction)<0) boundary=boundary-1
        select case(trim(faces(boundary)))
        case('periodic')
          position(direction)=modulo(query(direction),bounds(direction))
        case('inflow','outflow')
          position(direction)=max(0_int64,min(bounds(direction)-1,query(direction)))
          if (faces(boundary)=='inflow') reservoir=boundary
        case('slip_wall','no_slip_wall','isothermal_wall')
          position(direction)=-query(direction)-1
          if (query(direction)>=bounds(direction)) position(direction)=bounds(direction)-1-(query(direction)-bounds(direction))
          if (faces(boundary)/='slip_wall'.and.viscous) then
            flip=.not.flip;velocity_shift=-velocity_shift
            if (present(wall_velocity)) velocity_shift=velocity_shift+2*wall_velocity(:,boundary)
          else
            flip(direction)=.not.flip(direction);velocity_shift(direction)=-velocity_shift(direction)
          end if
        end select
      end do
      depth=requested_level
      do while(depth>=1)
        do peer=1,size(tree%nodes)
          if (tree%nodes(peer)%level+1/=depth) cycle
          if (any(position<origin(:,peer))) cycle
          if (any(position-origin(:,peer)>=int(tree%nodes(peer)%shape,int64))) cycle
          index=int(position-origin(:,peer))+1
          call reader%read(peer,index,.false.,value,success)
          if (.not.success) return
          if (present(coarse_fractions)) then
            blend=coarse_fractions(depth)
            if (blend==1) then
              call reader%read(peer,index,.true.,value,success)
              if (.not.success) return
            else if (blend>0) then
              call reader%read(peer,index,.true.,future_value,success)
              if (.not.success) return
              value=(1-blend)*value+blend*future_value
            end if
          end if
          success=.true.;exit
        end do
        if (success) exit
        position=position/int(ratios(depth),int64);depth=depth-1
      end do
      if (.not.success) return
      if (plan) return
      if (reservoir>0) value=inflow_state(:,reservoir)
      momentum=value(2:4)
      do direction=1,3
        if (flip(direction)) value(direction+1)=-value(direction+1)
      end do
      if (any(abs(velocity_shift)>0)) then
        density=value(1)
        if (density<=0) then
          success=.false.;return
        end if
        value(2:4)=value(2:4)+density*velocity_shift
        value(5)=value(5)+(sum(value(2:4)**2)-sum(momentum**2))/(2*density)
      end if
      success=all(ieee_is_finite(value))
    end subroutine
  end subroutine sample_amr_patch_tree_ghosts_3d
end module amr_patch_tree_ghost_3d_mod
