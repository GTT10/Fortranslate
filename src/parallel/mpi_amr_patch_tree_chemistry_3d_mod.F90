module mpi_amr_patch_tree_chemistry_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_chemistry_3d,advance_mpi_amr_tree_strang_3d,clone_mpi_amr_tree_inputs_3d
  public :: react_sdc_mpi_amr_tree_3d
  public :: initialize_mpi_amr_tree_reaction_source_3d
contains
  ! Owner-local counterpart of the serial fresh-start initializer. No public
  ! field is published until every level and rank accepts the complete trial.
  subroutine initialize_mpi_amr_tree_reaction_source_3d(tree,owners,species,reactions,state,temperature,dt,subcycling, &
      rtol,atol,reaction_source,comm,ok,evaluate_reaction,internal_energy,frame,domain_lower,root_spacing)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: subcycling
    type(amr_tree_state_3d), intent(inout) :: reaction_source
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    type(amr_tree_state_3d) :: candidate,temps,rates,sources
    type(amr_tree_state_3d), allocatable :: held
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,depth,nv,rank,status,mode,low,high
    real(dp) :: level_dt
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));mode=merge(1,0,subcycling)
    call MPI_Allreduce(mode,low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(mode,high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if(.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,valid)
    if(.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(comm,valid)
    if(.not.valid) return
    rates=state;sources=reaction_source
    if(present(internal_energy)) then
      allocate(held);held=internal_energy
    endif
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      rates%fields(node)%values=0;sources%fields(node)%values=0
    enddo
    level_dt=dt
    do depth=0,tree%level_count()-1
      call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,state,temperature,rates,level_dt,rtol,atol, &
        candidate,temps,sources,comm,valid,target_depth=depth,evaluate_reaction=evaluate_reaction, &
        internal_energy=held,old_internal_energy=internal_energy,frame=frame, &
        domain_lower=domain_lower,root_spacing=root_spacing)
      if(.not.valid) return
      if(depth==0) then
        call average_down_mpi_amr_tree_3d(tree,owners,sources,nv-4,comm,valid)
        if(.not.valid) return
      endif
      if(subcycling.and.depth<tree%level_count()-1) level_dt=level_dt/real(ratios(depth+2),dp)
    enddo
    call move_alloc(sources%fields,reaction_source%fields);ok=.true.
  end subroutine initialize_mpi_amr_tree_reaction_source_3d

  subroutine react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old_state,old_temperature,nonreacting_rates, &
      dt,rtol,atol,state,temperature,reaction_source,comm,ok,target_depth,evaluate_reaction,heat_release,internal_energy, &
      frame,domain_lower,root_spacing,old_internal_energy)
    type(amr_tree_state_3d), optional, intent(in) :: old_internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(amr_tree_state_3d), allocatable :: heat,auxiliary
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(in) :: old_state,old_temperature,nonreacting_rates
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature,reaction_source
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d) :: candidate,temps,sources
    real(dp) :: controls(3),low(3),high(3),work_stamp
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: frame_h(:,:),node_potential(:,:,:)
    real(dp) :: position(3)
    integer :: rank,status,ios,node,n(3),nv,depth,counts(7),min_counts(7),max_counts(7),i,j,k
    logical :: valid
    real(dp), allocatable :: node_heat(:,:,:),node_auxiliary(:,:,:),node_old_auxiliary(:,:,:)
    integer :: heat_supplied,heat_min,heat_max,heat_status
    ok=.false.;controls=[dt,rtol,atol];depth=-1;nv=reactive_nvar(size(species))
    heat_supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    valid=all(ieee_is_finite(controls)).and.min(dt,rtol,atol)>0
    valid=valid.and.(present(frame).eqv.present(domain_lower))
    if(present(frame)) then
      valid=valid.and.present(root_spacing)
      valid=valid.and.frame%is_valid()
      if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
      if(present(root_spacing)) then
        valid=valid.and.all(ieee_is_finite(root_spacing))
        valid=valid.and.all(root_spacing>0)
      endif
    endif
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0) valid=.false.
    end if
    call agree(comm,valid)
    if(.not.valid) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low/=high)) return
    counts=[size(species),size(reactions),depth,merge(1,0,present(evaluate_reaction)), &
      merge(1,0,present(internal_energy)),merge(1,0,present(frame)),merge(1,0,present(old_internal_energy))]
    call MPI_Allreduce(counts,min_counts,size(counts),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,max_counts,size(counts),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(min_counts/=max_counts)) return
    if(present(old_internal_energy)) then
      if(.not.present(internal_energy)) return
      call validate_mpi_amr_tree_local_3d(tree,owners,old_internal_energy,1,comm,valid)
      if(.not.valid) return
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(old_internal_energy%fields(node)%values))) valid=.false.
      enddo
      call agree(comm,valid)
      if(.not.valid) return
    endif
    if(present(frame)) then
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      call agree(comm,valid)
      if(.not.valid) return
      allocate(frame_h(3,size(tree%nodes)))
      frame_h(:,1)=root_spacing
      do node=2,size(tree%nodes)
        frame_h(:,node)=frame_h(:,tree%nodes(node)%parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      enddo
    endif
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),old_state,old_temperature, &
      candidate,temps,comm,rank,valid)
    if(.not.valid) return
    if(depth>=tree%level_count()) return
    call validate_mpi_amr_tree_local_3d(tree,owners,nonreacting_rates,nv,comm,valid)
    if(.not.valid) return
    if(depth>=0) then
      call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,valid)
      if(.not.valid) return
    end if
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(.not.all(ieee_is_finite(nonreacting_rates%fields(node)%values))) valid=.false.
      if(depth>=0) then
        if(.not.all(ieee_is_finite(reaction_source%fields(node)%values))) valid=.false.
      end if
    end do
    call agree(comm,valid)
    if(.not.valid) return
    allocate(sources%fields(size(tree%nodes)),stat=ios)
    if(ios==0) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(sources%fields(node)%values(nv-4,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) exit
        sources%fields(node)%values=0
        if(depth>=0) sources%fields(node)%values=reaction_source%fields(node)%values
      end do
    end if
    valid=ios==0;call agree(comm,valid)
    if(.not.valid) return
    ! All valid cells participate, including covered parents. Synchronization
    ! belongs to the SDC driver, after hydro/transport and child intervals.
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
      allocate(heat);heat=heat_release
    endif
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        valid=valid.and.all(ieee_is_finite(internal_energy%fields(node)%values))
      enddo
      call agree(comm,valid)
      if(.not.valid) return
      allocate(auxiliary);auxiliary=internal_energy
    endif
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      work_stamp=tree_work_start()
      if(present(heat_release)) allocate(node_heat(n(1),n(2),n(3)),source=0._dp)
      if(present(internal_energy)) allocate(node_auxiliary(n(1),n(2),n(3)),source=0._dp)
      if(present(old_internal_energy)) &
        allocate(node_old_auxiliary,source=old_internal_energy%fields(node)%values(1,:,:,:))
      if(present(frame)) then
        allocate(node_potential(n(1),n(2),n(3)))
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          position=domain_lower+(real(origin(:,node),dp)+real([i,j,k],dp)-.5_dp)*frame_h(:,node)
          call rotating_potential_3d(frame,position,node_potential(i,j,k),valid)
          if(.not.valid) exit
        enddo
        if(.not.valid) exit
        enddo
        if(.not.valid) exit
        enddo
        if(.not.valid) then
          call tree_work_stop(node,work_stamp)
          exit
        endif
      endif
      call react_sdc_state_3d(species,reactions,old_state%fields(node)%values, &
        old_temperature%fields(node)%values(1,:,:,:),nonreacting_rates%fields(node)%values,dt,rtol,atol, &
        candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:),sources%fields(node)%values,valid,evaluate_reaction, &
          heat_release=node_heat,internal_energy_density=node_auxiliary,rotating_potential=node_potential, &
          old_internal_energy_density=node_old_auxiliary)
      call tree_work_stop(node,work_stamp)
      if(.not.valid) exit
      if(allocated(node_potential)) deallocate(node_potential)
      if(allocated(node_old_auxiliary)) deallocate(node_old_auxiliary)
      if(present(internal_energy)) then
        auxiliary%fields(node)%values(1,:,:,:)=node_auxiliary
        deallocate(node_auxiliary)
      endif
      if(present(heat_release)) then
        heat%fields(node)%values(1,:,:,:)=node_heat
        deallocate(node_heat)
      endif
    end do
    call agree(comm,valid)
    if(.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    call move_alloc(sources%fields,reaction_source%fields)
    if(present(internal_energy)) internal_energy=auxiliary
    if(present(heat_release)) heat_release=heat
    ok=.true.
  end subroutine react_sdc_mpi_amr_tree_3d

  subroutine advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,state,temperature,interval,rtol,atol, &
      comm,ok,chemistry_integrator)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: interval,rtol,atol
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(amr_tree_state_3d) :: candidate,temps
    logical, allocatable :: active(:,:,:)
    real(dp) :: work_stamp
    real(dp) :: controls(3),low(3),high(3)
    character(len=64) :: method,reference_method
    integer :: rank,status,ios,counts(2),min_counts(2),max_counts(2)
    integer :: node,child,n(3),lo(3),hi(3)
    logical :: valid
    ok=.false.
    controls=[interval,rtol,atol]
    valid=all(ieee_is_finite(controls))
    if (valid) valid=interval>=0.and.min(rtol,atol)>0
    method='(default)'
    if (present(chemistry_integrator)) then
      if (len_trim(chemistry_integrator)>64) valid=.false.
      method=chemistry_integrator
    end if
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    counts=[size(species),size(reactions)]
    call MPI_Allreduce(counts,min_counts,2,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(counts,max_counts,2,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_counts/=max_counts)) return
    reference_method=method
    call MPI_Bcast(reference_method,64,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=method==reference_method;call agree(comm,valid)
    if (.not.valid) return
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(active(n(1),n(2),n(3)),stat=ios)
      if (ios/=0) then
        valid=.false.;exit
      end if
      active=.true.
      ! Only exposed cells own chemistry. Covered parents are supplied by
      ! restriction after every leaf and uncovered ancestor has succeeded.
      do child=node+1,size(tree%nodes)
        if (tree%nodes(child)%parent/=node) cycle
        lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
          tree%nodes(child)%patch%coarse_k_lower]
        hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
          tree%nodes(child)%patch%coarse_k_upper]
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      work_stamp=tree_work_start()
      call advance_reactive_chemistry_3d(species,reactions,candidate%fields(node)%values, &
        temps%fields(node)%values(1,:,:,:),n(1),n(2),n(3),interval,rtol,atol,valid,active,chemistry_integrator)
      call tree_work_stop(node,work_stamp)
      if (.not.valid) exit
      deallocate(active)
    end do
    call agree(comm,valid)
    if (.not.valid) return
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,reactive_nvar(size(species)),comm,valid)
    if (.not.valid) return
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,candidate,temps,comm,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    ok=.true.
  end subroutine advance_mpi_amr_tree_chemistry_3d

  subroutine advance_mpi_amr_tree_strang_3d(tree,owners,species,reactions,state,temperature,root_spacing,dt,faces, &
      reconstruction,limiter,solver,rtol,atol,boundary_change,maximum_correction,comm,ok, &
      chemistry_integrator,inflow_state,inflow_temperature,outlet_pressure,difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),outlet_pressure(6)
    real(dp), optional, intent(in) :: difmag
    type(amr_tree_state_3d) :: candidate,temps
    real(dp) :: budget(size(boundary_change)),largest
    integer :: rank
    logical :: valid
    ok=.false.
    valid=ieee_is_finite(dt).and.dt>0
    call agree(comm,valid)
    if (.not.valid) return
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if (.not.valid) return
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
      comm,valid,chemistry_integrator)
    if (.not.valid) return
    budget=0;largest=0
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,temps,root_spacing,dt,faces, &
      reconstruction,limiter,solver,budget,largest,comm,valid,inflow_state,inflow_temperature,outlet_pressure,difmag=difmag)
    if (.not.valid) return
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
      comm,valid,chemistry_integrator)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=budget;maximum_correction=largest;ok=.true.
  end subroutine advance_mpi_amr_tree_strang_3d
  subroutine clone_mpi_amr_tree_inputs_3d(tree,owners,nspecies,state,temperature,candidate,temps,comm,rank,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),nspecies
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), intent(out) :: candidate,temps
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: rank
    logical, intent(out) :: ok
    integer :: node,status,ios
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(nspecies),comm,ok)
    if (.not.ok) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ok)
    if (.not.ok) return
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    ok=.true.
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) ok=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) ok=.false.
      if (any(temperature%fields(node)%values<=0)) ok=.false.
    end do
    call agree(comm,ok)
    if (.not.ok) return
    allocate(candidate%fields(size(tree%nodes)),temps%fields(size(tree%nodes)),stat=ios)
    if (ios==0) then
      do node=1,size(tree%nodes)
        if (owners(node)/=rank) cycle
        allocate(candidate%fields(node)%values,source=state%fields(node)%values,stat=ios)
        if (ios/=0) exit
        allocate(temps%fields(node)%values,source=temperature%fields(node)%values,stat=ios)
        if (ios/=0) exit
      end do
    end if
    ok=ios==0;call agree(comm,ok)
  end subroutine
  subroutine agree(comm,condition)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(inout) :: condition
    logical :: result
    integer :: ierr
    call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    condition=.false.
    if (ierr==MPI_SUCCESS) condition=result
  end subroutine
end module mpi_amr_patch_tree_chemistry_3d_mod
