module amr_patch_tree_checkpoint_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: iso_c_binding, only: c_char,c_int,c_null_char
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_eos_policy_identity,valid_nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction,valid_elementary_reaction
  use gas_transport_mod, only: gas_transport_species,compatible_transport_database
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_difmag_3d
  use amr_reactive_3d_checkpoint_mod, only: write_selected_checkpoint_context,checkpoint_temperature_is_consistent
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    average_down_amr_patch_tree_3d
  use gas_transport_mod, only: write_transport_context
  implicit none
  private
  type, public :: amr_tree_checkpoint_3d
    type(amr_patch_tree_3d) :: tree
    type(amr_tree_state_3d) :: state,temperature,reaction_source,heat_release,internal_energy
    real(dp), allocatable :: initial(:),boundary(:),source(:)
    integer(int64), allocatable :: record_offsets(:)
    real(dp) :: time=0,maximum_reflux=0,minimum_theta=1,maximum_remap_error=0
    real(dp), allocatable :: energy_correction
    real(dp), allocatable :: pressure_work_change
    real(dp), allocatable :: previous_timestep
    logical :: post_regrid=.false.
    integer :: steps=0,regrids=0
  end type
  public :: checkpoint_tree_held_matches
  public :: valid_amr_tree_timestep_history_3d
  public :: checkpoint_record_offsets_3d
  public :: read_amr_tree_checkpoint_header_3d
  public :: build_amr_tree_checkpoint_context_3d,write_amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d
  ! Publish only fully closed checkpoints; retain the previous file on failure.
  interface
    function checkpoint_rename(old,new) bind(C,name='rename') result(status)
      import c_int,c_char
      character(kind=c_char), intent(in) :: old(*),new(*)
      integer(c_int) :: status
    end function
    function checkpoint_unlink(path) bind(C,name='unlink') result(status)
      import c_int,c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int) :: status
    end function
  end interface
contains
  subroutine build_amr_tree_checkpoint_context_3d(config,species,reactions,transport,mole_fractions,identity,integrator, &
      schedule_identity,context,ok)
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    character(len=*), intent(in) :: identity,integrator,schedule_identity
    character(len=:), allocatable, intent(inout) :: context
    logical, intent(out) :: ok
    type(reactive_3d_config) :: canonical
    character(len=65536) :: line
    character(len=:), allocatable :: buffer
    integer :: unit,status,bytes,cursor,length,i
    ok=.false.
    if (size(species)<1.or.size(mole_fractions)/=size(species)) return
    if(any(species%pinned_eos_extrapolation.neqv.config%pinned_eos_extrapolation)) return
    if (.not.all(ieee_is_finite(mole_fractions))) return
    if (any(mole_fractions<0).or.abs(sum(mole_fractions)-1)>5e-10_dp) return
    do i=1,size(species)
      if (.not.valid_nasa7_species(species(i))) return
    end do
    do i=1,size(reactions)
      if (.not.valid_elementary_reaction(reactions(i),size(species))) return
    end do
    if (.not.compatible_transport_database(species,transport)) return
    if (len_trim(identity)==0.or.(integrator/='explicit'.and.integrator/='implicit')) return
    if (.not.valid_context(trim(schedule_identity)//achar(10))) return
    ! Capture all configuration fields except continuation limits/output path.
    if(.not.ieee_is_finite(config%difmag)) return
    if(config%difmag<0.and.config%difmag/=-1._dp) return
    canonical=config;canonical%final_time=0;canonical%maximum_steps=0;canonical%output_file=''
    canonical%difmag=resolved_reactive_difmag_3d(config)
    if(canonical%reconstruction/='characteristic_ppm'.and.canonical%difmag==0) canonical%difmag=-1
    open(newunit=unit,status='scratch',access='stream',form='formatted',delim='quote',iostat=status)
    if (status/=0) return
    write(unit,*,iostat=status) &
      canonical%nx,canonical%ny,canonical%nz,canonical%maximum_steps, &
      canonical%sdc_iterations,canonical%sdc_reactor,canonical%time_integrator,canonical%mol_lax_friedrichs, &
      canonical%allow_small_energy,canonical%allow_negative_energy,canonical%dual_energy_update_E_from_e, &
      canonical%dual_energy_eta2,canonical%x_lower,canonical%x_upper,canonical%y_lower,canonical%y_upper, &
      canonical%z_lower,canonical%z_upper,canonical%final_time,canonical%cfl,canonical%problem, &
      canonical%reconstruction,canonical%limiter,canonical%thermo_model,canonical%riemann_solver, &
      canonical%boundary_condition,canonical%boundary_faces,canonical%wall_temperature,canonical%inflow_mach, &
      canonical%external_forcing,canonical%difmag,canonical%outlet_pressure,canonical%wall_velocity, &
      canonical%inflow_temperature,canonical%inflow_pressure,canonical%inflow_velocity, &
      canonical%inflow_composition_count,canonical%inflow_composition_species, &
      canonical%inflow_composition_mole_fractions,canonical%chemistry_enabled, &
      canonical%chemistry_relative_tolerance,canonical%chemistry_absolute_tolerance,canonical%transport_enabled, &
      canonical%viscosity_enabled,canonical%thermal_conduction_enabled,canonical%species_diffusion_enabled, &
      canonical%barodiffusion_enabled,canonical%transport_cfl,canonical%initial_temperature, &
      canonical%initial_pressure,canonical%initial_velocity_x,canonical%initial_velocity_y, &
      canonical%initial_velocity_z,canonical%density_wave_amplitude,canonical%wave_number_x, &
      canonical%wave_number_y,canonical%wave_number_z,canonical%hotspot_temperature_rise, &
      canonical%hotspot_center_x,canonical%hotspot_center_y,canonical%hotspot_center_z,canonical%hotspot_width, &
      canonical%x_h2,canonical%x_h,canonical%x_o,canonical%x_o2,canonical%x_oh,canonical%x_h2o,canonical%x_ho2, &
      canonical%x_h2o2,canonical%x_ar,canonical%x_n2,canonical%composition_count,canonical%composition_species, &
      canonical%composition_mole_fractions,canonical%output_file
    if(status==0.and.config%pinned_eos_extrapolation) &
      write(unit,'(a)',iostat=status) trim(nasa7_eos_policy_identity(species))
    if (status/=0) go to 900
    do i=1,size(species)
      write(unit,*,iostat=status) species(i)%name,species(i)%molecular_weight,species(i)%temperature_min, &
        species(i)%temperature_mid,species(i)%temperature_max,species(i)%low_coefficients,species(i)%high_coefficients
      if (status/=0) go to 900
      call write_transport_context(unit,transport(i),status)
      if (status/=0) go to 900
    end do
    call write_selected_checkpoint_context(unit,reactions,identity,integrator,mole_fractions, &
      config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,status)
    if (status/=0) go to 900
    write(unit,'(a)',iostat=status) schedule_identity
    if (status/=0) go to 900
    inquire(unit=unit,size=bytes,iostat=status)
    if (status/=0.or.bytes<1) go to 900
    allocate(character(len=bytes)::buffer);cursor=0
    rewind(unit,iostat=status)
    if (status/=0) go to 900
    do
      read(unit,'(a)',iostat=status) line
      if (status<0) exit
      if (status/=0.or.len_trim(line)==len(line)) go to 900
      length=len_trim(line)
      if (length+1>bytes-cursor) go to 900
      buffer(cursor+1:cursor+length+1)=trim(line)//achar(10);cursor=cursor+length+1
    end do
    context=buffer(:cursor);ok=.true.
900 close(unit)
  end subroutine

  subroutine write_amr_tree_checkpoint_3d(path,species,context,snapshot,ok,indexed,mol_energy_reset)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_checkpoint_3d), intent(in) :: snapshot
    logical, intent(out) :: ok
    logical, optional, intent(in) :: indexed,mol_energy_reset
    character(len=:), allocatable :: temporary
    integer(c_int) :: publish_status
    integer(int64) :: expected_position,actual_size
    logical :: with_index,valid
    integer(int64), allocatable :: offsets(:)
    integer :: unit,status,node,i,j,k,n(3),start,finish
    ok=.false.
    if (.not.valid_snapshot(snapshot,species,mol_energy_reset)) return
    if (.not.valid_context(context)) return
    with_index=.false.
    if (present(indexed)) with_index=indexed
    if (with_index) then
      call checkpoint_record_offsets_3d(snapshot%tree,size(snapshot%initial),offsets,valid, &
      allocated(snapshot%internal_energy%fields))
      if (.not.valid) return
    end if
    temporary=trim(path)//'.tmp'
    open(newunit=unit,file=temporary,status='new',action='write',access='stream', &
      form='formatted',iostat=status)
    if (status/=0) return
    if(allocated(snapshot%internal_energy%fields)) then
      if(allocated(snapshot%pressure_work_change)) then
        write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_INTERNAL_ENERGY_V2'
        if(status/=0) go to 900
        write(unit,'(2(es27.18e3,1x))',iostat=status) snapshot%energy_correction,snapshot%pressure_work_change
      else
      write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_INTERNAL_ENERGY_V1'
      if(status/=0) go to 900
      write(unit,'(es27.18e3)',iostat=status) snapshot%energy_correction
      endif
      if(status/=0) go to 900
    endif
    if(allocated(snapshot%previous_timestep)) then
      write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_TIMESTEP_HISTORY_V1'
      if(status/=0) go to 900
    endif
    if(allocated(snapshot%heat_release%fields)) then
      if(with_index) then
        write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V6'
      else
        write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V5'
      end if
    else if (allocated(snapshot%source)) then
      if (with_index) then
        write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V4'
      else
        write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V3'
      end if
    else if (with_index) then
      write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V2'
    else
      write(unit,'(a)',iostat=status) 'PELEF_AMR_TREE_3D_CHECKPOINT_V1'
    end if
    if (status/=0) go to 900
    write(unit,'(i0)',iostat=status) len(context)
    if (status/=0) go to 900
    start=1
    do while(start<=len(context))
      finish=start+index(context(start:),achar(10))-2
      write(unit,'(a)',iostat=status) context(start:finish)
      if (status/=0) go to 900
      start=finish+2
    end do
    write(unit,'(*(i0,1x))',iostat=status) size(snapshot%tree%nodes),size(snapshot%initial),snapshot%tree%nodes(1)%shape
    if (status/=0) go to 900
    write(unit,'(4(es27.18e3,1x),2(i0,1x))',iostat=status) snapshot%time,snapshot%maximum_reflux, &
      snapshot%minimum_theta,snapshot%maximum_remap_error,snapshot%steps,snapshot%regrids
    if (status/=0) go to 900
    if(allocated(snapshot%previous_timestep)) then
      write(unit,'(es27.18e3,1x,l1)',iostat=status) snapshot%previous_timestep,snapshot%post_regrid
      if(status/=0) go to 900
    endif
    write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%initial,snapshot%boundary
    if (status/=0) go to 900
    if (allocated(snapshot%source)) then
      write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%source
      if (status/=0) go to 900
      write(unit,'(i0)',iostat=status) merge(1,0,allocated(snapshot%reaction_source%fields))
      if (status/=0) go to 900
    end if
    do node=2,size(snapshot%tree%nodes)
      associate(p=>snapshot%tree%nodes(node)%patch)
        write(unit,'(*(i0,1x))',iostat=status) snapshot%tree%nodes(node)%parent,p%coarse_i_lower,p%coarse_i_upper, &
          p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper,p%refinement_ratio
      end associate
      if (status/=0) go to 900
    end do
    if (with_index) then
      do node=1,size(offsets)
        write(unit,'(i0)',iostat=status) offsets(node)
        if (status/=0) go to 900
      end do
    end if
    do node=1,size(snapshot%tree%nodes)
      n=snapshot%tree%nodes(node)%shape
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(allocated(snapshot%internal_energy%fields)) then
            write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%state%fields(node)%values(:,i,j,k), &
              snapshot%temperature%fields(node)%values(1,i,j,k),snapshot%internal_energy%fields(node)%values(1,i,j,k)
            else
            write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%state%fields(node)%values(:,i,j,k), &
              snapshot%temperature%fields(node)%values(1,i,j,k)
            endif
            if (status/=0) go to 900
          end do
        end do
      end do
    end do
    if (allocated(snapshot%reaction_source%fields)) then
      do node=1,size(snapshot%tree%nodes)
        n=snapshot%tree%nodes(node)%shape
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              write(unit,'(*(es27.18e3,1x))',iostat=status) snapshot%reaction_source%fields(node)%values(:,i,j,k)
              if (status/=0) go to 900
            end do
          end do
        end do
      end do
    end if
    if (allocated(snapshot%heat_release%fields)) then
      do node=1,size(snapshot%tree%nodes)
        n=snapshot%tree%nodes(node)%shape
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              write(unit,'(es27.18e3)',iostat=status) snapshot%heat_release%fields(node)%values(:,i,j,k)
              if (status/=0) go to 900
            end do
          end do
        end do
      end do
    end if
    write(unit,'(a)',iostat=status) 'END_AMR_TREE_CHECKPOINT'
    if(status==0) inquire(unit=unit,pos=expected_position,iostat=status)
900 close(unit,iostat=node)
    ok=status==0.and.node==0
    if(ok) then
      ! GNU formatted buffering can lose a failed write without setting IOSTAT.
      ! Check the closed payload length before replacing the public checkpoint.
      inquire(file=temporary,size=actual_size,iostat=status)
      ok=status==0
      if(ok) ok=actual_size==expected_position-1_int64
    endif
    if(ok) then
      publish_status=checkpoint_rename(temporary//c_null_char,trim(path)//c_null_char)
      ok=publish_status==0
    endif
    if(.not.ok) publish_status=checkpoint_unlink(temporary//c_null_char)
  end subroutine

  subroutine read_amr_tree_checkpoint_3d(path,species,context,root_shape,snapshot,ok,allow_internal_energy,mol_energy_reset)
    character(len=*), intent(in) :: path,context
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: root_shape(3)
    type(amr_tree_checkpoint_3d), intent(inout) :: snapshot
    logical, intent(out) :: ok
    logical, optional, intent(in) :: allow_internal_energy,mol_energy_reset
    type(amr_tree_checkpoint_3d) :: candidate
    character(len=65536) :: line
    integer :: unit,status,count,nvar,n(3),node,i,j,k
    integer(int64) :: payload_start,position
    logical :: valid
    ok=.false.
    if (.not.valid_context(context)) return
    open(newunit=unit,file=path,status='old',action='read',access='stream',form='formatted',iostat=status)
    if (status/=0) return
    call read_amr_tree_checkpoint_header_3d(unit,species,context,root_shape,candidate,valid, &
      allow_sdc=.true.,allow_internal_energy=allow_internal_energy)
    if (.not.valid) go to 900
    count=size(candidate%tree%nodes);nvar=size(candidate%initial)
    inquire(unit=unit,pos=payload_start,iostat=status)
    if (status/=0) go to 900
    allocate(candidate%state%fields(count),candidate%temperature%fields(count))
    do node=1,count
      n=candidate%tree%nodes(node)%shape
      allocate(candidate%state%fields(node)%values(nvar,n(1),n(2),n(3)), &
        candidate%temperature%fields(node)%values(1,n(1),n(2),n(3)),stat=status)
      if (status/=0) go to 900
      if(allocated(candidate%internal_energy%fields)) then
        allocate(candidate%internal_energy%fields(node)%values(1,n(1),n(2),n(3)),stat=status)
        if(status/=0) go to 900
      endif
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(allocated(candidate%internal_energy%fields)) then
            read(unit,*,iostat=status) candidate%state%fields(node)%values(:,i,j,k), &
              candidate%temperature%fields(node)%values(1,i,j,k),candidate%internal_energy%fields(node)%values(1,i,j,k)
            else
            read(unit,*,iostat=status) candidate%state%fields(node)%values(:,i,j,k), &
              candidate%temperature%fields(node)%values(1,i,j,k)
            endif
            if (status/=0) go to 900
          end do
        end do
      end do
      if (allocated(candidate%record_offsets)) then
        inquire(unit=unit,pos=position,iostat=status)
        if (status/=0) go to 900
        if (position-payload_start/=candidate%record_offsets(node+1)) go to 900
      end if
    end do
    if (allocated(candidate%record_offsets)) then
      inquire(unit=unit,pos=position,iostat=status)
      if (status/=0) go to 900
      if (position-payload_start/=candidate%record_offsets(count+1)) go to 900
    end if
    if (allocated(candidate%reaction_source%fields)) then
      do node=1,count
        n=candidate%tree%nodes(node)%shape
        allocate(candidate%reaction_source%fields(node)%values(nvar-4,n(1),n(2),n(3)),stat=status)
        if (status/=0) go to 900
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              read(unit,*,iostat=status) candidate%reaction_source%fields(node)%values(:,i,j,k)
              if (status/=0) go to 900
            end do
          end do
        end do
      end do
    end if
    if (allocated(candidate%heat_release%fields)) then
      do node=1,count
        n=candidate%tree%nodes(node)%shape
        allocate(candidate%heat_release%fields(node)%values(1,n(1),n(2),n(3)),stat=status)
        if (status/=0) go to 900
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              read(unit,*,iostat=status) candidate%heat_release%fields(node)%values(:,i,j,k)
              if (status/=0) go to 900
            end do
          end do
        end do
      end do
    end if
    read(unit,'(a)',iostat=status) line
    if (status/=0) go to 900
    if (trim(line)/='END_AMR_TREE_CHECKPOINT') go to 900
    read(unit,'(a)',iostat=status) line
    if (status>=0) go to 900
    if (.not.valid_snapshot(candidate,species,mol_energy_reset)) go to 900
    snapshot=candidate;ok=.true.
900 close(unit)
  end subroutine

  ! Read metadata; SDC V3/V4 requires explicit caller support.
  ! Read V1/V2 metadata from an already-open formatted unit.
  subroutine read_amr_tree_checkpoint_header_3d(unit,species,context,root_shape,snapshot,ok,allow_sdc,allow_internal_energy)
    logical, optional, intent(in) :: allow_sdc,allow_internal_energy
    integer, intent(in) :: unit,root_shape(3)
    type(nasa7_species), intent(in) :: species(:)
    character(len=*), intent(in) :: context
    type(amr_tree_checkpoint_3d), intent(inout) :: snapshot
    logical, intent(out) :: ok
    type(amr_tree_checkpoint_3d) :: candidate
    type(amr_patch_3d), allocatable :: patches(:)
    integer, allocatable :: parents(:),shapes(:,:)
    character(len=65536) :: line
    integer :: status,count,nvar,n(3),node,start,finish,context_length,lo(3),hi(3),ratio,parent,has_reaction
    integer(int64), allocatable :: expected_offsets(:)
    logical :: valid,indexed,sdc,with_heat,with_history,with_internal,with_work
    ok=.false.
    if (.not.valid_context(context)) return
    read(unit,'(a)',iostat=status) line
    if (status/=0) return
    with_work=trim(line)=='PELEF_AMR_TREE_INTERNAL_ENERGY_V2'
    with_internal=trim(line)=='PELEF_AMR_TREE_INTERNAL_ENERGY_V1'.or.with_work
    if(with_internal) then
      if(.not.present(allow_internal_energy)) return
      if(.not.allow_internal_energy) return
      allocate(candidate%energy_correction)
      if(with_work) then
        allocate(candidate%pressure_work_change)
        read(unit,*,iostat=status) candidate%energy_correction,candidate%pressure_work_change
        if(status/=0) return
        if(.not.ieee_is_finite(candidate%pressure_work_change)) return
      else
      read(unit,*,iostat=status) candidate%energy_correction
      endif
      if(status/=0) return
      if(.not.ieee_is_finite(candidate%energy_correction)) return
      read(unit,'(a)',iostat=status) line
      if(status/=0) return
    endif
    with_history=trim(line)=='PELEF_AMR_TREE_TIMESTEP_HISTORY_V1'
    if(with_history) then
      read(unit,'(a)',iostat=status) line
      if(status/=0) return
    endif
    with_heat=trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V5'.or.trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V6'
    sdc=trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V3'.or.trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V4'.or.with_heat
    if (sdc) then
      if (.not.present(allow_sdc)) return
      if (.not.allow_sdc) return
    end if
    indexed=trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V2'.or.trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V4'.or. &
      trim(line)=='PELEF_AMR_TREE_3D_CHECKPOINT_V6'
    if (.not.sdc.and..not.indexed.and.trim(line)/='PELEF_AMR_TREE_3D_CHECKPOINT_V1') return
    read(unit,*,iostat=status) context_length
    if (status/=0) return
    if (context_length/=len(context)) return
    start=1
    do while(start<=len(context))
      finish=start+index(context(start:),achar(10))-2
      read(unit,'(a)',iostat=status) line
      if (status/=0) return
      if (trim(line)/=context(start:finish)) return
      start=finish+2
    end do
    read(unit,*,iostat=status) count,nvar,n
    if (status/=0) return
    if (count<1.or.nvar/=reactive_nvar(size(species)).or.any(n/=root_shape)) return
    allocate(patches(count-1),parents(count-1),shapes(3,count),stat=status)
    if (status/=0) return
    if(with_heat) allocate(candidate%heat_release%fields(count))
    if(with_internal) allocate(candidate%internal_energy%fields(count))
    shapes(:,1)=n
    read(unit,*,iostat=status) candidate%time,candidate%maximum_reflux,candidate%minimum_theta, &
      candidate%maximum_remap_error,candidate%steps,candidate%regrids
    if (status/=0) return
    if(with_history) then
      allocate(candidate%previous_timestep)
      read(unit,*,iostat=status) candidate%previous_timestep,candidate%post_regrid
      if(status/=0) return
    endif
    if(.not.valid_amr_tree_timestep_history_3d(candidate)) return
    allocate(candidate%initial(nvar),candidate%boundary(nvar))
    read(unit,*,iostat=status) candidate%initial,candidate%boundary
    if (status/=0) return
    if (sdc) then
      allocate(candidate%source(nvar))
      read(unit,*,iostat=status) candidate%source
      if (status/=0) return
      read(unit,*,iostat=status) has_reaction
      if (status/=0) return
      if (has_reaction<0.or.has_reaction>1) return
      if (has_reaction==1) allocate(candidate%reaction_source%fields(count))
    end if
    do node=2,count
      read(unit,*,iostat=status) parent,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),ratio
      if (status/=0) return
      if (parent<1.or.parent>=node) return
      n=shapes(:,parent);parents(node-1)=parent
      call initialize_amr_patch_3d(n(1),n(2),n(3),lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),ratio,patches(node-1),valid)
      if (.not.valid) return
      if (any(int(hi-lo+1,int64)*int(ratio,int64)>int(huge(1),int64))) return
      shapes(:,node)=int(int(hi-lo+1,int64)*int(ratio,int64))
    end do
    call initialize_amr_patch_tree_3d(root_shape,parents,patches,candidate%tree,valid)
    if (.not.valid) return
    if (indexed) then
      call checkpoint_record_offsets_3d(candidate%tree,nvar,expected_offsets,valid,with_internal)
      if (.not.valid) return
      allocate(candidate%record_offsets(size(expected_offsets)),stat=status)
      if (status/=0) return
      do node=1,size(expected_offsets)
        read(unit,*,iostat=status) candidate%record_offsets(node)
        if (status/=0) return
      end do
      if (any(candidate%record_offsets/=expected_offsets)) return
    end if
    snapshot=candidate;ok=.true.
  end subroutine

  logical function valid_amr_tree_timestep_history_3d(s) result(valid)
    type(amr_tree_checkpoint_3d),intent(in) :: s
    valid=.false.
    if(.not.allocated(s%previous_timestep)) then
      valid=.not.s%post_regrid
      return
    endif
    if(.not.all(ieee_is_finite([s%time,s%previous_timestep]))) return
    if(s%previous_timestep<0.or.s%previous_timestep>s%time) return
    if((s%steps==0).neqv.(s%previous_timestep==0)) return
    valid=.true.
  end function valid_amr_tree_timestep_history_3d

  subroutine checkpoint_record_offsets_3d(tree,nvar,offsets,ok,with_internal_energy)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nvar
    integer(int64), allocatable, intent(out) :: offsets(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: with_internal_energy
    integer(int64) :: cells,width
    integer :: node,axis,status
    ok=.false.
    if (.not.tree%is_valid().or.nvar<1) return
    width=28_int64*(int(nvar,int64)+1_int64)
    if(present(with_internal_energy)) then
      if(with_internal_energy) width=width+28_int64
    endif
    allocate(offsets(size(tree%nodes)+1),stat=status)
    if (status/=0) return
    offsets(1)=0
    do node=1,size(tree%nodes)
      cells=1
      do axis=1,3
        if (cells>huge(cells)/int(tree%nodes(node)%shape(axis),int64)) return
        cells=cells*int(tree%nodes(node)%shape(axis),int64)
      end do
      if (cells>(huge(cells)-offsets(node))/width) return
      offsets(node+1)=offsets(node)+cells*width
    end do
    ok=.true.
  end subroutine

  logical function valid_context(context) result(valid)
    character(len=*), intent(in) :: context
    integer :: start,finish
    valid=.false.
    if (len(context)<1) return
    if (context(len(context):)/=achar(10)) return
    start=1
    do while(start<=len(context))
      finish=start+index(context(start:),achar(10))-2
      if (finish-start+1>=65536) return
      start=finish+2
    end do
    valid=.true.
  end function

  ! SDC/MOL reset covered-cell internal energy from restricted total energy.
  ! Restriction and subtraction of kinetic energy do not commute. The explicit
  ! nonrotating SDC/MOL caller permits that identity, without modifying the payload
  ! or weakening the restriction check for conserved state or the EOS check.
  logical function checkpoint_tree_held_matches(state,held,restricted,mol_energy_reset) result(valid)
    real(dp), intent(in) :: state(:,:,:,:),held(:,:,:),restricted(:,:,:)
    logical, optional, intent(in) :: mol_energy_reset
    real(dp) :: rho,velocity(3),kinetic,thermal,tolerance,limit
    integer :: i,j,k
    logical :: reset
    valid=.false.;reset=.false.
    if(present(mol_energy_reset)) reset=mol_energy_reset
    if(size(state,1)<5) return
    if(any(shape(state(1,:,:,:))/=shape(held)).or.any(shape(restricted)/=shape(held))) return
    if(.not.all(ieee_is_finite(held)).or..not.all(ieee_is_finite(restricted))) return
    do k=1,size(held,3)
    do j=1,size(held,2)
    do i=1,size(held,1)
      tolerance=64*epsilon(1._dp)*max(1._dp,abs(held(i,j,k)))
      if(abs(restricted(i,j,k)-held(i,j,k))<=tolerance) cycle
      if(.not.reset) return
      if(.not.all(ieee_is_finite(state(1:5,i,j,k)))) return
      rho=state(1,i,j,k)
      if(rho<=0) return
      limit=sqrt(huge(1._dp)/16)
      if(rho<1) then
        if(any(abs(state(2:4,i,j,k))>limit*rho)) return
      endif
      velocity=state(2:4,i,j,k)/rho
      if(any(abs(velocity)>limit)) return
      kinetic=.5_dp*sum(velocity**2)
      if(rho>1) then
        if(kinetic>huge(1._dp)/rho/4) return
      endif
      kinetic=rho*kinetic
      if(abs(state(5,i,j,k))>huge(1._dp)/2) return
      thermal=state(5,i,j,k)-kinetic
      if(abs(thermal-held(i,j,k))>tolerance) return
    enddo
    enddo
    enddo
    valid=.true.
  end function checkpoint_tree_held_matches

  logical function valid_snapshot(s,species,mol_energy_reset) result(valid)
    type(amr_tree_checkpoint_3d), intent(in) :: s
    type(nasa7_species), intent(in) :: species(:)
    logical, optional, intent(in) :: mol_energy_reset
    type(amr_tree_state_3d) :: synchronized
    integer :: node,nvar
    logical :: consistent
    valid=.false.;nvar=reactive_nvar(size(species))
    if(.not.valid_amr_tree_timestep_history_3d(s)) return
    if (.not.s%state%is_valid(s%tree).or..not.s%temperature%is_valid(s%tree)) return
    if (size(s%state%fields(1)%values,1)/=nvar.or.size(s%temperature%fields(1)%values,1)/=1) return
    if (.not.allocated(s%initial).or..not.allocated(s%boundary)) return
    if (size(s%initial)/=nvar.or.size(s%boundary)/=nvar) return
    if (.not.all(ieee_is_finite(s%initial)).or..not.all(ieee_is_finite(s%boundary))) return
    if (.not.all(ieee_is_finite([s%time,s%maximum_reflux,s%minimum_theta,s%maximum_remap_error]))) return
    if (min(s%time,s%maximum_reflux,s%minimum_theta,s%maximum_remap_error)<0.or.s%minimum_theta>1) return
    if (min(s%steps,s%regrids)<0.or.s%regrids>s%steps) return
    if ((s%steps==0.and.s%time/=0).or.(s%steps>0.and.s%time<=0)) return
    if (allocated(s%source)) then
      if (size(s%source)/=nvar) return
      if (.not.all(ieee_is_finite(s%source))) return
    end if
    if (allocated(s%reaction_source%fields)) then
      if (.not.allocated(s%source)) return
      if (.not.s%reaction_source%is_valid(s%tree)) return
      do node=1,size(s%tree%nodes)
        if (size(s%reaction_source%fields(node)%values,1)/=nvar-4) return
        if (.not.all(ieee_is_finite(s%reaction_source%fields(node)%values))) return
      end do
    end if
    if (allocated(s%heat_release%fields)) then
      if (.not.allocated(s%source)) return
      if (.not.s%heat_release%is_valid(s%tree)) return
      do node=1,size(s%tree%nodes)
        if (size(s%heat_release%fields(node)%values,1)/=1) return
        if (.not.all(ieee_is_finite(s%heat_release%fields(node)%values))) return
      end do
    end if
    if(allocated(s%internal_energy%fields).neqv.allocated(s%energy_correction)) return
    if(allocated(s%pressure_work_change)) then
      if(.not.allocated(s%internal_energy%fields)) return
      if(.not.ieee_is_finite(s%pressure_work_change)) return
    endif
    if(allocated(s%internal_energy%fields)) then
      if(.not.ieee_is_finite(s%energy_correction)) return
      if(.not.s%internal_energy%is_valid(s%tree)) return
      if(size(s%internal_energy%fields(1)%values,1)/=1) return
      synchronized=s%internal_energy
      call average_down_amr_patch_tree_3d(s%tree,synchronized,consistent)
      if(.not.consistent) return
      do node=1,size(s%tree%nodes)
        if(.not.checkpoint_tree_held_matches(s%state%fields(node)%values, &
          s%internal_energy%fields(node)%values(1,:,:,:),synchronized%fields(node)%values(1,:,:,:), &
          mol_energy_reset)) return
      enddo
    endif
    synchronized=s%state
    call average_down_amr_patch_tree_3d(s%tree,synchronized,consistent)
    if (.not.consistent) return
    do node=1,size(s%tree%nodes)
      if (any(abs(synchronized%fields(node)%values-s%state%fields(node)%values)> &
          64*epsilon(1.0_dp)*max(1.0_dp,abs(s%state%fields(node)%values)))) return
      if(allocated(s%internal_energy%fields)) then
      call checkpoint_temperature_is_consistent(species,s%state%fields(node)%values, &
        s%temperature%fields(node)%values(1,:,:,:),consistent, &
          s%internal_energy%fields(node)%values(1,:,:,:))
      else
      call checkpoint_temperature_is_consistent(species,s%state%fields(node)%values, &
        s%temperature%fields(node)%values(1,:,:,:),consistent)
      endif
      if (.not.consistent) return
    end do
    valid=.true.
  end function
end module amr_patch_tree_checkpoint_3d_mod
