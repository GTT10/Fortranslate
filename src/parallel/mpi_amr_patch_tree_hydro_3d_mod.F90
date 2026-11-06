module mpi_amr_patch_tree_hydro_3d_mod
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: compute_reactive_face_fluxes_3d,compute_reactive_plm_face_fluxes_3d, &
    compute_reactive_ppm_halo_ctu_fluxes_3d,adjust_reactive_artificial_viscosity_3d
  use reactive_mol_3d_mod, only: compute_reactive_mol_halo_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_ghost_3d_mod, only: apply_amr_tree_pressure_ghosts_3d,amr_tree_inflow_faces_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_hydro_3d,compute_sdc_mpi_amr_tree_hydro_3d,compute_mol_mpi_amr_tree_hydro_3d
contains
  subroutine compute_mol_mpi_amr_tree_hydro_3d(tree,owners,species,state,temperature,root_spacing, &
      dt,faces,hydro_rate,integrated_flux,comm,ok,inflow_state,inflow_temperature,outlet_pressure, &
      inflow_mach,inflow_pressure, &
      target_depth,coarse_start,coarse_end,coarse_fractions,order,lax_friedrichs, &
      internal_energy,internal_energy_rate,internal_energy_integrated_flux,frame,domain_lower)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(amr_tree_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),outlet_pressure(6)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: coarse_fractions(:)
    integer, optional, intent(in) :: order
    logical, optional, intent(in) :: lax_friedrichs
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: internal_energy_integrated_flux(:)
    type(amr_tree_state_3d) :: local_state,local_temperature
    real(dp) :: budget(5+size(species)),maximum
    integer :: reconstruction
    reconstruction=2
    if(present(order)) reconstruction=order
    local_state=state;local_temperature=temperature
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local_state,local_temperature,root_spacing,dt,faces, &
      'characteristic_ppm','mc','pelec',budget,maximum,comm,ok,inflow_state,inflow_temperature, &
      outlet_pressure,inflow_mach,inflow_pressure,predictor_source=state,hydro_rate=hydro_rate, &
      integrated_flux=integrated_flux,target_depth=target_depth,coarse_start=coarse_start,coarse_end=coarse_end, &
      coarse_fractions=coarse_fractions,mol_order=reconstruction,mol_lax_friedrichs=lax_friedrichs, &
      internal_energy=internal_energy,internal_energy_rate=internal_energy_rate, &
      internal_energy_integrated_flux=internal_energy_integrated_flux,frame=frame,domain_lower=domain_lower)
  end subroutine compute_mol_mpi_amr_tree_hydro_3d

  subroutine compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,state,temperature,predictor_source,root_spacing, &
      dt,faces,hydro_rate,integrated_flux,comm,ok,inflow_state,inflow_temperature,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag, &
      target_depth,coarse_start,coarse_end,coarse_fractions,source_start,source_end,source_fractions, &
      internal_energy,internal_energy_rate,internal_energy_integrated_flux,internal_energy_pressure_work)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,predictor_source
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(amr_tree_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),outlet_pressure(6)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6),difmag
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end,source_start,source_end
    real(dp), optional, intent(in) :: coarse_fractions(:),source_fractions(:)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_pressure_work
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: internal_energy_integrated_flux(:)
    type(amr_tree_state_3d) :: local_state,local_temperature
    real(dp) :: budget(5+size(species)),maximum
    local_state=state;local_temperature=temperature
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local_state,local_temperature,root_spacing,dt,faces, &
      'characteristic_ppm','mc','pelec',budget,maximum,comm,ok,inflow_state,inflow_temperature, &
      outlet_pressure,inflow_mach,inflow_pressure,difmag,predictor_source,hydro_rate,integrated_flux, &
      target_depth,coarse_start,coarse_end,coarse_fractions,source_start,source_end,source_fractions, &
      internal_energy=internal_energy,internal_energy_rate=internal_energy_rate, &
      internal_energy_integrated_flux=internal_energy_integrated_flux,internal_energy_pressure_work=internal_energy_pressure_work)
  end subroutine compute_sdc_mpi_amr_tree_hydro_3d

  subroutine advance_mpi_amr_tree_hydro_3d(tree,owners,species,state,temperature,root_spacing,dt,faces, &
      reconstruction,limiter,solver,boundary_change,maximum_correction,comm,ok,inflow_state,inflow_temperature, &
      outlet_pressure,inflow_mach,inflow_pressure,difmag,predictor_source,hydro_rate,integrated_flux, &
      target_depth,coarse_start,coarse_end,coarse_fractions,source_start,source_end,source_fractions,mol_order,mol_lax_friedrichs, &
      internal_energy,internal_energy_rate,internal_energy_integrated_flux,frame,domain_lower,internal_energy_pressure_work)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    type(amr_tree_state_3d), optional, intent(in) :: predictor_source
    type(amr_tree_state_3d), optional, intent(inout) :: hydro_rate
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end,source_start,source_end
    real(dp), optional, intent(in) :: coarse_fractions(:),source_fractions(:)
    integer, optional, intent(in) :: mol_order
    logical, optional, intent(in) :: mol_lax_friedrichs
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: internal_energy_integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_pressure_work
    real(dp) :: pressure(6),mach(6),pin(6),tin(6),viscosity
    type(amr_tree_state_3d) :: initial,candidate,temps,stage_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)),first_max,second_max
    integer :: node,nvar,np,rank,status,ios,parent,width,depth,selection(15),selection_low(15),selection_high(15)
    integer :: reconstruction_order,face
    real(dp) :: controls(37),low(37),high(37),frame_lower(3)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(10),reference_options(10)
    character(len=:), allocatable :: donor_limiter
    logical :: valid,stage_only,mol,lf,stored
    ok=.false.
    stage_only=present(predictor_source)
    reconstruction_order=0;lf=.false.
    if(present(mol_order)) reconstruction_order=mol_order
    if(present(mol_lax_friedrichs)) lf=mol_lax_friedrichs
    mol=present(mol_order);stored=present(internal_energy)
    nvar=reactive_nvar(size(species));np=nvar+1
    if(stored) np=np+1
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    depth=-1
    if(present(target_depth)) depth=target_depth
    selection=[depth,merge(1,0,present(coarse_start)),merge(1,0,present(coarse_end)), &
      merge(1,0,present(coarse_fractions)),merge(1,0,present(source_start)), &
      merge(1,0,present(source_end)),merge(1,0,present(source_fractions)),reconstruction_order,merge(1,0,lf), &
      merge(1,0,stored),merge(1,0,present(internal_energy_rate)),merge(1,0,present(internal_energy_integrated_flux)), &
      merge(1,0,present(frame)),merge(1,0,present(domain_lower)),merge(1,0,present(internal_energy_pressure_work))]
    call MPI_Allreduce(selection,selection_low,size(selection),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selection,selection_high,size(selection),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(selection_low/=selection_high)) return
    valid=depth>=-1.and.depth<tree%level_count()
    valid=valid.and.all(selection(10:12)==selection(10))
    if(stored) then
      valid=valid.and.stage_only
      if(.not.mol) valid=valid.and.reconstruction=='characteristic_ppm'.and.solver=='pelec'
    endif
    if(present(internal_energy_pressure_work)) valid=valid.and.stored.and..not.mol
    valid=valid.and.reconstruction_order>=0.and.reconstruction_order<=2
    if(mol) valid=valid.and.stage_only.and.reconstruction_order>0
    if(present(frame)) then
      valid=valid.and.mol.and.present(domain_lower).and.frame%is_valid()
    endif
    if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
    if(present(target_depth)) valid=valid.and.depth>=0.and.stage_only
    valid=valid.and.all(selection(2:4)==selection(2)).and.all(selection(5:7)==selection(5))
    if(any(selection(2:7)>0)) valid=valid.and.depth>=0.and.stage_only
    call agree(valid)
    if(.not.valid) return
    if(stored) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        valid=valid.and.all(ieee_is_finite(internal_energy%fields(node)%values))
      enddo
      call agree(valid)
      if(.not.valid) return
    endif
    viscosity=.1_dp
    if(present(difmag)) viscosity=difmag
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    mach=-1;pin=0;tin=0
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_pressure)) pin=inflow_pressure
    valid=.true.
    if(present(inflow_temperature)) then
      valid=size(inflow_temperature)==6
      if(valid) then
        where(mach>=0) tin=inflow_temperature
      end if
    end if
    valid=valid.and.(.not.present(inflow_pressure).or.present(inflow_mach))
    valid=valid.and.all(ieee_is_finite(mach)).and.all(mach<=100)
    valid=valid.and..not.any(mach>=0.and.faces/='inflow')
    valid=valid.and.all(ieee_is_finite(pin)).and.all(pin>=0).and.all(pin<=1e15_dp)
    valid=valid.and..not.any(pin>0.and.mach<0)
    valid=valid.and.all(ieee_is_finite(tin)).and..not.any(tin<=0.and.mach>=0)
    call agree(valid)
    if(.not.valid) return
    controls=0
    controls(:29)=[root_spacing,dt,pressure,mach,tin,pin,viscosity]
    if(present(frame)) then
      controls(30)=real(frame%axis,dp)
      controls(31)=frame%angular_speed
      controls(32:34)=frame%axis_point
    endif
    if(present(domain_lower)) controls(35:37)=domain_lower
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0)
    valid=valid.and.viscosity>=0
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    valid=valid.and.size(boundary_change)==nvar.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.(reconstruction=='pcm'.or.reconstruction=='characteristic_plm'.or. &
      reconstruction=='characteristic_ppm')
    valid=valid.and.(limiter=='mc'.or.limiter=='minmod').and.len_trim(solver)<=64
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) valid=.false.
      if (any(temperature%fields(node)%values<=0)) valid=.false.
    end do
    call agree(valid)
    if (.not.valid) return
    call MPI_Allreduce(controls,low,37,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,37,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    valid=present(predictor_source).eqv.present(hydro_rate)
    valid=valid.and.(present(hydro_rate).eqv.present(integrated_flux))
    if(stage_only) valid=valid.and.reconstruction=='characteristic_ppm'.and.solver=='pelec'.and.limiter=='mc'
    call agree(valid)
    if(.not.valid) return
    options(:6)=faces;options(7)=reconstruction;options(8)=limiter;options(9)=solver
    options(10)='advance'
    if(stage_only) options(10)='sdc-stage'
    if(mol) options(10)='mol-stage'
    reference_options=options
    call MPI_Bcast(reference_options,64*10,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(options==reference_options)
    call agree(valid)
    if (.not.valid) return
    if(stage_only) then
      call validate_mpi_amr_tree_local_3d(tree,owners,predictor_source,nvar,comm,valid)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(predictor_source%fields(node)%values))) valid=.false.
      end do
      call agree(valid)
      if(.not.valid) return
    end if
    if (any(faces=='inflow'.and.mach<0)) then
      valid=present(inflow_state).and.present(inflow_temperature)
      if (valid) then
        valid=all(shape(inflow_state)==[nvar,6]).and.size(inflow_temperature)==6
        if (valid) valid=all(ieee_is_finite(inflow_state)).and.all(ieee_is_finite(inflow_temperature))
      end if
      call agree(valid)
      if (.not.valid) return
      allocate(inlet(np,6),reference_inlet(np,6),stat=ios)
      valid=ios==0;call agree(valid)
      if (.not.valid) return
      inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
      if(stored) then
        valid=.true.;inlet(np,:)=0
        do face=1,6
          if(inflow_state(1,face)>0) then
            inlet(np,face)=inflow_state(5,face)-.5_dp*sum(inflow_state(2:4,face)**2)/inflow_state(1,face)
          else if(faces(face)=='inflow'.and.mach(face)<0) then
            valid=.false.
          endif
        enddo
        call agree(valid)
        if(.not.valid) return
        valid=all(ieee_is_finite(inlet(np,:)));call agree(valid)
        if(.not.valid) return
      endif
      reference_inlet=inlet
      call MPI_Bcast(reference_inlet,size(reference_inlet),MPI_DOUBLE_PRECISION,0,comm,status)
      if (status/=MPI_SUCCESS) return
      valid=all(inlet==reference_inlet);call agree(valid)
      if (.not.valid) return
    end if
    if (reconstruction/='pcm') donor_limiter=limiter
    if (mol.or.stored) donor_limiter='linear'
    width=2
    if (reconstruction=='characteristic_ppm') width=4
    if(mol) width=3
    allocate(h(3,size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    h(:,1)=root_spacing
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if (any(h(:,node)<=0)) return
    end do
    if(present(frame)) then
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      call agree(valid)
      if(.not.valid) return
    endif
    call copy_local(state,initial,valid)
    if (.not.valid) return
    call copy_local(temperature,temps,valid)
    if (.not.valid) return
    if(.not.stage_only) then
      call average_down_mpi_amr_tree_3d(tree,owners,initial,nvar,comm,valid)
      if (.not.valid) return
      call recover_tree(initial,temps,valid)
      if (.not.valid) return
    end if
    call copy_local(initial,candidate,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,first_budget,first_max,valid)
    if (.not.valid) return
    if(stage_only) then
      ok=.true.
      return
    end if
    if (reconstruction=='characteristic_ppm') then
      call move_alloc(candidate%fields,state%fields)
      call move_alloc(temps%fields,temperature%fields)
      boundary_change=first_budget;maximum_correction=first_max;ok=.true.
      return
    end if
    call copy_local(temps,stage_temperature,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,second_budget,second_max,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      candidate%fields(node)%values=.5_dp*initial%fields(node)%values+.5_dp*candidate%fields(node)%values
    end do
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    if (.not.valid) return
    call move_alloc(stage_temperature%fields,temps%fields)
    call recover_tree(candidate,temps,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=.5_dp*(first_budget+second_budget)
    maximum_correction=max(first_max,second_max);ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success)
    end subroutine
    subroutine euler_stage(u,t,budget,largest,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      real(dp), intent(out) :: budget(:),largest
      real(dp) :: work_stamp
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated,ghosts,source_ghosts,predictor,auxiliary_rates,pressure_work
      type(amr_tree_flux_3d), allocatable :: auxiliary_flux(:)
      type(eb_integrated_flux_3d) :: auxiliary_faces
      type(amr_tree_flux_3d), allocatable :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),source_values(:,:,:,:)
      real(dp), allocatable :: ax(:,:,:),ay(:,:,:),az(:,:,:),pdiv(:,:,:)
      character(len=32) :: source_faces(6)
      real(dp) :: correction(nvar),rhs(nvar),theta,mol_boundary(nvar),flux_scale
      integer :: peer,n(3),dims(3),i,j,k,allocation_status
      success=.false.;budget=0;largest=0;correction=0
      allocate(packed%fields(size(tree%nodes)),flux(size(tree%nodes)),stat=allocation_status)
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(np,n(1),n(2),n(3)),stat=allocation_status)
        if (allocation_status/=0) exit
        packed%fields(peer)%values(:nvar,:,:,:)=u%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=t%fields(peer)%values(1,:,:,:)
        if(stored) packed%fields(peer)%values(np,:,:,:)=internal_energy%fields(peer)%values(1,:,:,:)
      end do
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      if(stored) then
        call copy_local(internal_energy,auxiliary_rates,success)
        if(.not.success) return
        if(present(internal_energy_pressure_work)) then
          call copy_local(internal_energy,pressure_work,success)
          if(.not.success) return
          do peer=1,size(tree%nodes)
            if(owners(peer)==rank) pressure_work%fields(peer)%values=0
          enddo
        endif
        allocate(auxiliary_flux(size(tree%nodes)),stat=allocation_status)
        success=allocation_status==0;call agree(success)
        if(.not.success) return
        do peer=1,size(tree%nodes)
          if(owners(peer)==rank) auxiliary_rates%fields(peer)%values=0
        enddo
      endif
      call copy_local(u,updated,success)
      if (.not.success) return
      if(stage_only) then
        do peer=1,size(tree%nodes)
          if(owners(peer)==rank) updated%fields(peer)%values=0
        end do
        predictor=predictor_source
        if(present(coarse_start)) then
          call set_ancestors(packed,coarse_start,np,success)
          if(.not.success) return
        end if
        if(present(source_start)) then
          call set_ancestors(predictor,source_start,nvar,success)
          if(.not.success) return
        end if
      end if
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,ghosts,comm,success, &
        viscous_wall=mol.and..not.present(frame), &
        inflow_state=inlet,coarse_limiter=donor_limiter,species_count=size(species),ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        depth=target_depth,coarse_future=coarse_end,coarse_fractions=coarse_fractions, &
        limit_internal_energy=mol.or.stored,stored_internal_energy=stored, &
        frame=frame,domain_lower=domain_lower,root_spacing=root_spacing)
      if (.not.success) return
      if(stage_only.and..not.mol) then
        source_faces=faces
        where(source_faces/='periodic') source_faces='outflow'
        call fill_mpi_amr_tree_ghosts_3d(tree,owners,predictor,nvar,source_faces,source_ghosts,comm,success, &
          coarse_limiter='mc',ghost_width=4,depth=target_depth,coarse_future=source_end,coarse_fractions=source_fractions)
        if(.not.success) return
      end if
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        n=tree%nodes(peer)%shape;dims=n+4
        frame_lower=0
        if(present(frame)) frame_lower=domain_lower+real(origin(:,peer),dp)*h(:,peer)
        if(stage_only.and..not.mol) source_values=source_ghosts%fields(peer)%values
        work_stamp=tree_work_start()
        associate(ghost=>ghosts%fields(peer)%values)
          call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,pressure,ghost,success,ghost_width=width, &
            stored_internal_energy=stored,frame=frame,domain_lower=domain_lower,root_spacing=root_spacing)
          if (.not.success) exit
          allocate(fx(nvar,dims(1),dims(2),dims(3)),fy(nvar,dims(1),dims(2),dims(3)), &
            fz(nvar,dims(1),dims(2),dims(3)),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          if(mol) then
            fx=0;fy=0;fz=0;mol_boundary=0
            if(stored) then
            call compute_reactive_mol_halo_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer),dt, &
              updated%fields(peer)%values,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
              fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2), &
              mol_boundary,success,reconstruction_order,lf, &
              internal_energy_rate=auxiliary_rates%fields(peer)%values(1,:,:,:), &
              internal_energy_integrated_flux=auxiliary_faces,frame=frame,domain_lower=frame_lower, &
              internal_energy_ghost=ghost(np,:,:,:))
              if(.not.success) exit
              allocate(auxiliary_flux(peer)%x(1,n(1)+1,n(2),n(3)), &
                auxiliary_flux(peer)%y(1,n(1),n(2)+1,n(3)),auxiliary_flux(peer)%z(1,n(1),n(2),n(3)+1), &
                stat=allocation_status)
              if(allocation_status/=0) then
                success=.false.;exit
              endif
              auxiliary_flux(peer)%x=auxiliary_faces%x
              auxiliary_flux(peer)%y=auxiliary_faces%y
              auxiliary_flux(peer)%z=auxiliary_faces%z
            else
            call compute_reactive_mol_halo_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer),dt, &
              updated%fields(peer)%values,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
              fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2), &
              mol_boundary,success,reconstruction_order,lf,frame=frame,domain_lower=frame_lower)
            endif
          else if (reconstruction=='characteristic_ppm') then
            block
              logical :: inflow_mask(6),wall_mask(6),reflecting_mask(6),outflow_mask(6)
              call amr_tree_inflow_faces_3d(tree,peer,faces,inflow_mask,success,wall_mask,reflecting_mask,outflow_mask)
              if (.not.success) exit
              fx=0;fy=0;fz=0
              if(stored) then
                allocate(ax(0:n(1),n(2),n(3)),ay(n(1),0:n(2),n(3)),az(n(1),n(2),0:n(3)), &
                  pdiv(n(1),n(2),n(3)),stat=allocation_status)
                if(allocation_status/=0) then
                  success=.false.;exit
                endif
                call compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
                  n,h(:,peer),dt,.false.,.true.,solver,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
                  fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),theta,success, &
                  inflow_faces=inflow_mask,no_slip_faces=wall_mask,difmag=viscosity, &
                  reflecting_faces=reflecting_mask,outflow_faces=outflow_mask,outlet_pressure=outlet_pressure, &
                  conserved_source=source_values,internal_energy_ghost=ghost(np,:,:,:), &
                  internal_energy_xflux=ax,internal_energy_yflux=ay,internal_energy_zflux=az, &
                  internal_energy_pressure_work=pdiv)
                if(.not.success) exit
                do k=1,n(3)
                do j=1,n(2)
                do i=1,n(1)
                  auxiliary_rates%fields(peer)%values(1,i,j,k)=-(ax(i,j,k)-ax(i-1,j,k))/h(1,peer)- &
                    (ay(i,j,k)-ay(i,j-1,k))/h(2,peer)-(az(i,j,k)-az(i,j,k-1))/h(3,peer)-pdiv(i,j,k)
                enddo
                enddo
                enddo
                allocate(auxiliary_flux(peer)%x(1,n(1)+1,n(2),n(3)), &
                  auxiliary_flux(peer)%y(1,n(1),n(2)+1,n(3)),auxiliary_flux(peer)%z(1,n(1),n(2),n(3)+1), &
                  stat=allocation_status)
                if(allocation_status/=0) then
                  success=.false.;exit
                endif
                auxiliary_flux(peer)%x(1,:,:,:)=dt*ax
                auxiliary_flux(peer)%y(1,:,:,:)=dt*ay
                auxiliary_flux(peer)%z(1,:,:,:)=dt*az
                if(present(internal_energy_pressure_work)) pressure_work%fields(peer)%values(1,:,:,:)=pdiv
                deallocate(ax,ay,az,pdiv)
              else
            call compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              n,h(:,peer),dt,.false.,.true.,solver,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
              fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),theta,success, &
              inflow_faces=inflow_mask,no_slip_faces=wall_mask,difmag=viscosity, &
              reflecting_faces=reflecting_mask,outflow_faces=outflow_mask,outlet_pressure=outlet_pressure, &
              conserved_source=source_values)
              endif
            end block
          else if (reconstruction=='pcm') then
            call compute_reactive_face_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              dims(1),dims(2),dims(3),solver,fx,fy,fz,success)
          else
            call compute_reactive_plm_face_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              dims(1),dims(2),dims(3),limiter,solver,fx,fy,fz,success)
          end if
          if (.not.success) exit
          if(viscosity>0.and.reconstruction/='characteristic_ppm') then
            block
              logical :: inflow_mask(6),wall_mask(6),reflecting_mask(6),outflow_mask(6)
              call amr_tree_inflow_faces_3d(tree,peer,faces,inflow_mask,success,wall_mask,reflecting_mask,outflow_mask)
              if (.not.success) exit
              call adjust_reactive_artificial_viscosity_3d(ghost(:nvar,2:n(1)+3,2:n(2)+3,2:n(3)+3), &
                n,h(:,peer),viscosity,wall_mask,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
                fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),success)
            end block
            if (.not.success) exit
          end if
          if(.not.mol) then
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                rhs=-(fx(:,i+2,j+2,k+2)-fx(:,i+1,j+2,k+2))/h(1,peer) &
                  -(fy(:,i+2,j+2,k+2)-fy(:,i+2,j+1,k+2))/h(2,peer) &
                  -(fz(:,i+2,j+2,k+2)-fz(:,i+2,j+2,k+1))/h(3,peer)
                if(stage_only) then
                  updated%fields(peer)%values(:,i,j,k)=rhs
                else
                  updated%fields(peer)%values(:,i,j,k)=u%fields(peer)%values(:,i,j,k)+dt*rhs
                end if
              end do
            end do
          end do
          endif
          flux_scale=dt
          if(mol) flux_scale=1
          allocate(flux(peer)%x,source=flux_scale*fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%y,source=flux_scale*fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%z,source=flux_scale*fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          deallocate(fx,fy,fz)
        end associate
        call tree_work_stop(peer,work_stamp)
      end do
      call agree(success)
      if (.not.success) return
      if(stage_only) then
        do peer=1,size(tree%nodes)
          if(owners(peer)/=rank) cycle
          if(.not.all(ieee_is_finite(updated%fields(peer)%values))) success=.false.
          if(stored) then
            if(.not.all(ieee_is_finite(auxiliary_rates%fields(peer)%values))) success=.false.
          endif
          if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
          if(stored) then
            if(.not.all(ieee_is_finite(auxiliary_flux(peer)%x)).or..not.all(ieee_is_finite(auxiliary_flux(peer)%y)).or. &
              .not.all(ieee_is_finite(auxiliary_flux(peer)%z))) success=.false.
          endif
          if(.not.all(ieee_is_finite(flux(peer)%x)).or..not.all(ieee_is_finite(flux(peer)%y)).or. &
            .not.all(ieee_is_finite(flux(peer)%z))) success=.false.
        end do
        call agree(success)
        if(.not.success) return
        if(stored) then
          call move_alloc(auxiliary_rates%fields,internal_energy_rate%fields)
          call move_alloc(auxiliary_flux,internal_energy_integrated_flux)
          if(present(internal_energy_pressure_work)) &
            call move_alloc(pressure_work%fields,internal_energy_pressure_work%fields)
        endif
        call move_alloc(updated%fields,hydro_rate%fields)
        call move_alloc(flux,integrated_flux)
        return
      end if
      if (rank==owners(1)) then
        n=tree%nodes(1)%shape
        do k=1,n(3)
          do j=1,n(2)
            budget=budget+(flux(1)%x(:,1,j,k)-flux(1)%x(:,n(1)+1,j,k))*h(2,1)*h(3,1)
          end do
        end do
        do k=1,n(3)
          do i=1,n(1)
            budget=budget+(flux(1)%y(:,i,1,k)-flux(1)%y(:,i,n(2)+1,k))*h(1,1)*h(3,1)
          end do
        end do
        do j=1,n(2)
          do i=1,n(1)
            budget=budget+(flux(1)%z(:,i,j,1)-flux(1)%z(:,i,j,n(3)+1))*h(1,1)*h(2,1)
          end do
        end do
      end if
      call MPI_Bcast(budget,nvar,MPI_DOUBLE_PRECISION,owners(1),comm,status)
      if (status/=MPI_SUCCESS) then
        success=.false.;return
      end if
      call reflux_mpi_amr_tree_3d(tree,owners,updated,nvar,flux,root_spacing,faces,correction,largest,comm,success)
      if (.not.success) return
      budget=budget+correction
      call recover_tree(updated,t,success)
      if (.not.success) return
      call move_alloc(updated%fields,u%fields)
    end subroutine
    subroutine set_ancestors(field,start,ncomponents,success)
      type(amr_tree_state_3d), intent(inout) :: field
      type(amr_tree_state_3d), intent(in) :: start
      integer, intent(in) :: ncomponents
      logical, intent(out) :: success
      integer :: ancestor
      call validate_mpi_amr_tree_local_3d(tree,owners,start,ncomponents,comm,success,maximum_depth=depth-1)
      if(.not.success) return
      do ancestor=1,size(tree%nodes)
        if(owners(ancestor)/=rank.or.tree%nodes(ancestor)%level>=depth) cycle
        if(.not.all(ieee_is_finite(start%fields(ancestor)%values))) success=.false.
      end do
      call agree(success)
      if(.not.success) return
      do ancestor=1,size(tree%nodes)
        if(owners(ancestor)/=rank.or.tree%nodes(ancestor)%level>=depth) cycle
        field%fields(ancestor)%values=start%fields(ancestor)%values
      end do
    end subroutine set_ancestors
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
    subroutine copy_local(source,destination,success)
      type(amr_tree_state_3d), intent(in) :: source
      type(amr_tree_state_3d), intent(inout) :: destination
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: copy
      integer :: peer,ierr
      allocate(copy%fields(size(tree%nodes)),stat=ierr)
      if (ierr==0) then
        do peer=1,size(tree%nodes)
          if (owners(peer)/=rank) cycle
          allocate(copy%fields(peer)%values,source=source%fields(peer)%values,stat=ierr)
          if (ierr/=0) exit
        end do
      end if
      success=ierr==0;call agree(success)
      if (success) call move_alloc(copy%fields,destination%fields)
    end subroutine
  end subroutine advance_mpi_amr_tree_hydro_3d
end module
