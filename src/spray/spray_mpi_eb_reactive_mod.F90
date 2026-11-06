module spray_mpi_eb_reactive_mod
  use mpi_f08
  use precision_mod, only: dp
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use eb_reactive_patch_tree_3d_mod, only: pack_eb_held_auxiliary_state_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use spray_mpi_amr_balance_mod, only: sum_mpi_spray_rotating_potential
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,valid_spray_options,valid_spray_parcel
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling, &
    advance_mpi_eb_recursive_population_stage,spray_sdc_population_stage
  use spray_mpi_amr_tree_mod, only: spray_level_population,initialize_mpi_amr_sdc_populations, &
    collect_mpi_amr_level_particles,redistribute_mpi_amr_subtree_particles
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_particle_predictor_halo_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d,advance_mpi_eb_tree_sdc_3d, &
    advance_mpi_eb_tree_sdc_subcycle_3d,mpi_eb_sdc_diffusion_source,mpi_eb_sdc_spatial_diffusion_source
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options
  use spray_mpi_amr_film_mod, only: spray_level_film_population,initialize_mpi_eb_film_levels,migrate_mpi_amr_wall_films
  use spray_mpi_amr_film_step_mod, only: advance_mpi_amr_film_coupling,spray_sdc_film_stage, &
    advance_mpi_eb_recursive_film_stage,finish_mpi_eb_film_level,append_mpi_eb_film_births
  use spray_mpi_mod, only: split_mpi_spray_parcels
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d, &
    integrate_mpi_eb_tree_3d,average_down_mpi_eb_tree_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use mpi_amr_patch_tree_transport_3d_mod, only: agree_mpi_constant_transport_coefficients
  implicit none
  private
  public :: advance_mpi_eb_spray_reactive,advance_mpi_eb_spray_sdc,advance_mpi_eb_spray_sdc_subcycle
contains
  ! Shared common/recursive level gas, particle and stationary film coordinator;
  ! outer publication is atomic.
  subroutine advance_mpi_eb_spray_sdc_subcycle(comm,tree,owners,geometry,extended_geometry,phi, &
      species,reactions,transport,liquids,options,state,temperature,parcels,component_mass,lower,spacing, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary_change,wall_change, &
      chemical_change,source_change,minimum_theta,particle_steps,last_id,ok,reaction_source,difmag, &
      redistribution,target_volume_fraction,gas_particle_change,level_steps,gas_acceleration, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,particle_acceleration,split_created, &
      films,film_mass,film_settings,film_steps,film_impacts,model_rejected,evaluate_diffusion,evaluate_reaction, &
      embedded_wall_temperature,embedded_noslip,mol_mode,mol_lax_friedrichs,frame, &
      state_redist_max_order,evaluate_spatial_diffusion,spatial_diffusion_halo_width,common_step, &
      constant_coefficients,heat_release,frame_energy_change, &
      internal_energy,energy_policy,energy_correction,internal_energy_pressure_work_change)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction,internal_energy_pressure_work_change
    real(dp), optional, intent(in) :: constant_coefficients(:)
    logical, optional, intent(in) :: common_step
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(inout) :: frame_energy_change
    real(dp) :: potential_before,potential_after,frame_source,frame_boundary,frame_change
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluate_diffusion
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: evaluate_spatial_diffusion
    integer, optional, intent(in) :: spatial_diffusion_halo_width
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,force(3),rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, intent(out) :: particle_steps
    integer(int64), intent(inout) :: last_id
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release
    real(dp), optional, intent(in) :: difmag,target_volume_fraction,gas_acceleration(3)
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(inout) :: gas_particle_change(:)
    integer(int64), optional, intent(inout) :: level_steps(:)
    integer(int64), optional, intent(out) :: split_created
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: embedded_noslip,mol_mode,mol_lax_friedrichs
    real(dp), optional, intent(in) :: wall_temperature(6),particle_acceleration(3)
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    logical, optional, intent(in) :: film_impacts
    logical, optional, intent(out) :: model_rejected
    logical :: impacts,rejected
    type(spray_level_film_population), allocatable :: film_physical(:),film_virtuals(:),film_ghosts(:)
    type(spray_sdc_film_stage), allocatable :: film_work(:)
    type(spray_wall_film), allocatable :: final_films(:)
    real(dp), allocatable :: final_film_mass(:,:)
    integer :: film_flag,film_low,film_high,film_count,cursor
    type(spray_level_population), allocatable :: physical(:),virtuals(:),ghosts(:)
    type(spray_sdc_population_stage), allocatable :: work(:)
    type(amr_tree_state_3d) :: gas,temp,metrics,centroids,particle_first,particle_last,source_correction
    integer :: particle_iterations,geometry_padding,predictor_width
    logical :: mol,common
    integer :: common_flag,common_low,common_high
    type(amr_tree_state_3d), allocatable :: ir,heat,stored
    type(amr_tree_state_3d) :: zero_held
    type(eb_integrated_flux_3d), allocatable :: held_flux(:)
    type(dual_energy_policy) :: held_policy
    real(dp), allocatable :: correction_prefix
    real(dp) :: pressure_work_change
    integer :: held_flag,held_low,held_high,naux
    real(dp), allocatable :: exchange(:)
    type(spray_parcel), allocatable :: final(:)
    real(dp), allocatable :: masses(:,:)
    integer(int64), allocatable :: counts(:)
    integer(int64) :: watermark,low_id,high_id,maximum_id,global_id,total_steps,before_split
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),particle_wall(size(source_change))
    real(dp) :: low_multiplicity,high_multiplicity
    real(dp) :: exterior(size(source_change)),particle_force_total(size(source_change)),theta,h(3),origin(3)
    integer :: nv,levels,rank,status,node,n(3),i,j,k,p,ios,frame_flag,frame_min,frame_max
    integer :: redist_order,redist_low,redist_high
    logical :: valid,agreed,gas_accepted,periodic(3)
    ok=.false.;particle_steps=0
    held_flag=merge(1,0,present(internal_energy))+2*merge(1,0,present(energy_policy))+ &
      4*merge(1,0,present(energy_correction))+8*merge(1,0,present(internal_energy_pressure_work_change))
    call MPI_Allreduce(held_flag,held_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_flag,held_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_low/=held_high) return
    valid=(.not.present(energy_policy).or.present(energy_correction)).and. &
      (.not.present(energy_correction).or.present(internal_energy))
    call agree();if(.not.agreed) return
    held_policy=dual_energy_policy();pressure_work_change=0
    if(present(energy_policy)) held_policy=energy_policy
    if(present(energy_correction)) then
      allocate(correction_prefix,stat=ios)
      valid=ios==0;call agree();if(.not.agreed) return
      correction_prefix=energy_correction
    endif
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,agreed)
    if(.not.agreed) return
    common=.false.
    if(present(common_step)) common=common_step
    common_flag=merge(1,0,common)
    call MPI_Allreduce(common_flag,common_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(common_flag,common_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.common_low/=common_high) return
    mol=.false.
    if(present(mol_mode)) mol=mol_mode
    call MPI_Allreduce(merge(1,0,mol),common_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(merge(1,0,mol),common_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.common_low/=common_high) return
    frame_flag=merge(1,0,present(frame))+2*merge(1,0,present(frame_energy_change))
    call MPI_Allreduce(frame_flag,frame_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame_flag,frame_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.frame_min/=frame_max) then
      ok=.false.;particle_steps=0
      return
    end if
    if(present(frame).and..not.mol) then
      ok=.false.;particle_steps=0
      return
    end if
    potential_before=0;potential_after=0;frame_source=0;frame_boundary=0;frame_change=0
    particle_iterations=iterations
    if(mol) then
      particle_iterations=1
    end if
    ok=.false.;particle_steps=0;nv=reactive_nvar(size(species))
    redist_order=0
    if(present(state_redist_max_order)) redist_order=state_redist_max_order
    valid=redist_order==0.or.redist_order==2.or.redist_order==4
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(redist_order,redist_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(redist_order,redist_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.redist_low/=redist_high) return
    impacts=.false.
    if(present(film_impacts)) impacts=film_impacts
    if(present(model_rejected)) model_rejected=.false.
    if(present(split_created)) split_created=0
    if(present(film_steps)) film_steps=0
    film_flag=merge(1,0,present(films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))+ &
      8*merge(1,0,impacts)
    call MPI_Allreduce(film_flag,film_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,film_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.film_low/=film_high) return
    if(film_flag/=0.and.film_flag/=7.and.film_flag/=15) return
    valid=allocated(parcels).and.allocated(component_mass).and.tree%is_valid().and.valid_spray_options(options)
    if(present(films)) valid=valid.and.allocated(films).and.allocated(film_mass)
    valid=valid.and.all(boundary_faces=='periodic'.or.boundary_faces=='outflow'.or. &
      boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
      boundary_faces=='inflow'.or.boundary_faces=='isothermal_wall')
    valid=valid.and.all([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]==nv)
    if(present(gas_particle_change)) valid=valid.and.size(gas_particle_change)==nv
    periodic=boundary_faces([1,3,5])=='periodic'
    valid=valid.and.all(periodic.eqv.(boundary_faces([2,4,6])=='periodic'))
    call agree();if(.not.agreed) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    levels=tree%level_count()
    if(present(level_steps)) valid=size(level_steps)==levels
    call agree();if(.not.agreed) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,low_id,1,MPI_INTEGER8,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,high_id,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low_id/=high_id) return
    maximum_id=max(0_int64,maxval(parcels%id))
    call MPI_Allreduce(maximum_id,global_id,1,MPI_INTEGER8,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.last_id<global_id) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,agreed)
    if(.not.agreed) return
    valid=all(lower==[geometry(1)%x_lower,geometry(1)%y_lower,geometry(1)%z_lower]).and. &
      all(spacing==[geometry(1)%dx,geometry(1)%dy,geometry(1)%dz])
    call agree();if(.not.agreed) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,agreed)
    if(.not.agreed) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,agreed)
    if(.not.agreed) return
    naux=nv+merge(1,0,present(internal_energy))
    allocate(exchange(naux),stat=ios)
    valid=ios==0;call agree();if(.not.agreed) return
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,agreed)
      if(.not.agreed) return
      allocate(stored,stat=ios)
      valid=ios==0;call agree();if(.not.agreed) return
      stored=internal_energy;zero_held=internal_energy
      do node=1,size(owners)
        if(owners(node)==rank) zero_held%fields(node)%values=0
      enddo
    endif
    valid=size(extended_geometry)==size(owners)
    call agree();if(.not.agreed) return
    geometry_padding=(extended_geometry(1)%nx-tree%nodes(1)%shape(1))/2
    valid=geometry_padding==6.or.(.not.mol.and.(geometry_padding==9.or.geometry_padding==10))
    call agree();if(.not.agreed) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,agreed,padding=geometry_padding)
    if(.not.agreed) return
    predictor_width=geometry_padding-3
    allocate(work(levels),counts(levels),metrics%fields(size(owners)), &
        centroids%fields(size(owners)),stat=ios)
    valid=ios==0;call agree();if(.not.agreed) return
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;h=[geometry(node)%dx,geometry(node)%dy,geometry(node)%dz]
      origin=[geometry(node)%x_lower,geometry(node)%y_lower,geometry(node)%z_lower]
      allocate(metrics%fields(node)%values(1,1-predictor_width:n(1)+predictor_width, &
          1-predictor_width:n(2)+predictor_width,1-predictor_width:n(3)+predictor_width), &
        centroids%fields(node)%values(3,1-predictor_width:n(1)+predictor_width, &
          1-predictor_width:n(2)+predictor_width,1-predictor_width:n(3)+predictor_width),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction( &
        4:n(1)+2*geometry_padding-3,4:n(2)+2*geometry_padding-3,4:n(3)+2*geometry_padding-3)
      do k=1-predictor_width,n(3)+predictor_width
      do j=1-predictor_width,n(2)+predictor_width
      do i=1-predictor_width,n(1)+predictor_width
        centroids%fields(node)%values(:,i,j,k)=origin+(real([i,j,k],dp)-.5_dp+ &
          [extended_geometry(node)%cell_centroid_x(i+geometry_padding,j+geometry_padding,k+geometry_padding), &
           extended_geometry(node)%cell_centroid_y(i+geometry_padding,j+geometry_padding,k+geometry_padding), &
           extended_geometry(node)%cell_centroid_z(i+geometry_padding,j+geometry_padding,k+geometry_padding)])*h
      end do
      end do
      end do
    end do
    call agree();if(.not.agreed) return
    if(present(frame)) then
      call sum_mpi_spray_rotating_potential(comm,frame,parcels,potential_before,agreed,films)
      if(.not.agreed) return
    endif
    call initialize_mpi_amr_sdc_populations(comm,tree,owners,parcels,component_mass,lower,spacing, &
      periodic,2,options%maximum_parcels,physical,virtuals,ghosts,agreed,phi)
    if(.not.agreed) return
    if(present(films)) then
      call initialize_mpi_eb_film_levels(comm,tree,owners,films,film_mass,lower,spacing, &
        film_settings%maximum_films,phi,film_physical,agreed,periodic)
      if(.not.agreed) return
      allocate(film_virtuals(levels),film_ghosts(levels),film_work(levels))
      do p=1,levels
        allocate(film_virtuals(p)%films(0),film_virtuals(p)%component_mass(size(liquids),0), &
          film_virtuals(p)%nodes(0),film_virtuals(p)%probes(3,0))
        film_ghosts(p)=film_virtuals(p)
      end do
    end if
    gas=state;temp=temperature;particle_first=state
    do node=1,size(owners)
      if(owners(node)==rank) particle_first%fields(node)%values=0
    end do
    source_correction=particle_first
    particle_force_total=0
    particle_last=particle_first;particle_wall=0;exterior=0;watermark=last_id
    if(options%etab_parcel_factor>=0) before_split=watermark
    if(present(reaction_source)) then
      allocate(ir);ir=reaction_source
    end if
    if(present(heat_release)) then
      allocate(heat);heat=heat_release
    end if
    if(allocated(stored)) then
    if(common) then
      counts=1_int64
      if(mol) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,gas,temp,extended_geometry, &
        dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary,wall,chemical,external,theta, &
        comm,gas_accepted,ir,difmag,redistribution=redistribution,target_volume_fraction=target_volume_fraction, &
        auxiliary_source=common_particle_source,auxiliary_change=exchange, &
        gas_acceleration=gas_acceleration,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,wall_temperature=wall_temperature, &
        evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
        evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat, &
        internal_energy=stored,internal_energy_flux=held_flux,energy_policy=energy_policy, &
        energy_correction=correction_prefix,internal_energy_pressure_work_change=pressure_work_change, &
        auxiliary_has_internal_energy=.true.)
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,gas,temp,extended_geometry, &
        dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary,wall,chemical,external,theta, &
        comm,gas_accepted,ir,difmag,redistribution=redistribution,target_volume_fraction=target_volume_fraction, &
        auxiliary_source=common_particle_source,auxiliary_change=exchange,predictor_builder=common_particle_predictor, &
        gas_acceleration=gas_acceleration,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,wall_temperature=wall_temperature, &
        evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
        evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat, &
        internal_energy=stored,internal_energy_flux=held_flux,energy_policy=energy_policy, &
        energy_correction=correction_prefix,internal_energy_pressure_work_change=pressure_work_change, &
        auxiliary_has_internal_energy=.true.)
      end if
    else if(mol) then
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,gas,temp, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary,wall,chemical,external,theta,comm,gas_accepted,ir,difmag,level_steps=counts,redistribution=redistribution, &
      target_volume_fraction=target_volume_fraction,halo_source=particle_source,auxiliary_change=exchange, &
      gas_acceleration=gas_acceleration,inflow_state=inflow_state, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,post_level=particle_post_level,evaluate_diffusion=evaluate_diffusion, &
      evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
      embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
      constant_coefficients=constant_coefficients,heat_release=heat, &
        internal_energy=stored,internal_energy_flux=held_flux,energy_policy=energy_policy, &
        energy_correction=correction_prefix,internal_energy_pressure_work_change=pressure_work_change, &
        auxiliary_has_internal_energy=.true.)
    else
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,gas,temp, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary,wall,chemical,external,theta,comm,gas_accepted,ir,difmag,level_steps=counts,redistribution=redistribution, &
      target_volume_fraction=target_volume_fraction,halo_source=particle_source,auxiliary_change=exchange, &
      gas_acceleration=gas_acceleration,predictor_builder=particle_predictor,inflow_state=inflow_state, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,post_level=particle_post_level,evaluate_diffusion=evaluate_diffusion, &
      evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
      embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
      constant_coefficients=constant_coefficients,heat_release=heat, &
        internal_energy=stored,internal_energy_flux=held_flux,energy_policy=energy_policy, &
        energy_correction=correction_prefix,internal_energy_pressure_work_change=pressure_work_change, &
        auxiliary_has_internal_energy=.true.)
    end if
    else
    if(common) then
      counts=1_int64
      if(mol) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,gas,temp,extended_geometry, &
        dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary,wall,chemical,external,theta, &
        comm,gas_accepted,ir,difmag,redistribution=redistribution,target_volume_fraction=target_volume_fraction, &
        auxiliary_source=common_particle_source,auxiliary_change=exchange, &
        gas_acceleration=gas_acceleration,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,wall_temperature=wall_temperature, &
        evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
        evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat)
      else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,owners,geometry,gas,temp,extended_geometry, &
        dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary,wall,chemical,external,theta, &
        comm,gas_accepted,ir,difmag,redistribution=redistribution,target_volume_fraction=target_volume_fraction, &
        auxiliary_source=common_particle_source,auxiliary_change=exchange,predictor_builder=common_particle_predictor, &
        gas_acceleration=gas_acceleration,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,wall_temperature=wall_temperature, &
        evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame,state_redist_max_order=redist_order, &
        evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat)
      end if
    else if(mol) then
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,gas,temp, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary,wall,chemical,external,theta,comm,gas_accepted,ir,difmag,level_steps=counts,redistribution=redistribution, &
      target_volume_fraction=target_volume_fraction,halo_source=particle_source,auxiliary_change=exchange, &
      gas_acceleration=gas_acceleration,inflow_state=inflow_state, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,post_level=particle_post_level,evaluate_diffusion=evaluate_diffusion, &
      evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
      embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
      constant_coefficients=constant_coefficients,heat_release=heat)
    else
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,owners,geometry,gas,temp, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,extended_geometry,boundary_faces, &
      boundary,wall,chemical,external,theta,comm,gas_accepted,ir,difmag,level_steps=counts,redistribution=redistribution, &
      target_volume_fraction=target_volume_fraction,halo_source=particle_source,auxiliary_change=exchange, &
      gas_acceleration=gas_acceleration,predictor_builder=particle_predictor,inflow_state=inflow_state, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      wall_temperature=wall_temperature,post_level=particle_post_level,evaluate_diffusion=evaluate_diffusion, &
      evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
      embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
      constant_coefficients=constant_coefficients,heat_release=heat)
    end if
    endif
    ! Callbacks update agreed; keep the gas result separate so an early gas
    ! rejection cannot be overwritten by the last successful particle callback.

    if(.not.gas_accepted) return
    call synchronize_particle_sources()

    if(.not.agreed) return
    total_steps=0
    do p=1,levels
      valid=work(p)%next_stage==0.and.abs(work(p)%time+work(p)%interval-dt)<=32*epsilon(dt)*dt
      valid=valid.and.counts(p)<=int(huge(particle_steps)/2,int64)-total_steps
      call agree();if(.not.agreed) return
      total_steps=total_steps+counts(p)
    end do
    call collect_mpi_amr_level_particles(comm,tree,owners,physical,lower,spacing,periodic, &
      options%maximum_parcels,final,masses,agreed,phi)

    if(.not.agreed) return
    if(present(films)) then
      ! Only physical populations are returned; virtual/ghost payloads never enter
      ! the accepted film inventory. Reuse migration to validate the flat output.
      film_count=0
      do p=1,levels
        valid=size(film_physical(p)%films)<=film_settings%maximum_films-film_count
        call agree();if(.not.agreed) return
        film_count=film_count+size(film_physical(p)%films)
      end do
      allocate(final_films(film_count),final_film_mass(size(liquids),film_count));cursor=0
      do p=1,levels
        film_count=size(film_physical(p)%films)
        final_films(cursor+1:cursor+film_count)=film_physical(p)%films
        final_film_mass(:,cursor+1:cursor+film_count)=film_physical(p)%component_mass
        cursor=cursor+film_count
      end do
      call migrate_mpi_amr_wall_films(comm,tree,owners,final_films,final_film_mass,lower,spacing, &
        film_settings%maximum_films,agreed,phi,periodic=periodic)
      if(.not.agreed) return
    end if
    if(options%etab_parcel_factor<0) before_split=watermark
    if(options%maximum_multiplicity>0.and.options%etab_parcel_factor<0) then
      call split_mpi_spray_parcels(comm,final,options%maximum_multiplicity,options%maximum_parcels,watermark,agreed,masses)
      if(.not.agreed) return
    end if
    if(present(frame)) then
      call sum_mpi_spray_rotating_potential(comm,frame,final,potential_after,agreed,final_films)
      if(.not.agreed) return
      frame_change=frame_source+potential_before-potential_after-frame_boundary
      valid=ieee_is_finite(frame_change);call agree();if(.not.agreed) return
    endif
    valid=safe_sum(boundary,exterior).and.safe_sum(wall,particle_wall).and.safe_sum(external,particle_force_total)
    call agree();if(.not.agreed) return
    boundary=boundary+exterior;wall=wall+particle_wall;external=external+particle_force_total
    if(present(internal_energy)) call move_alloc(stored%fields,internal_energy%fields)
    if(present(energy_correction)) energy_correction=correction_prefix
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=pressure_work_change
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(final,parcels);call move_alloc(masses,component_mass)
    if(present(films)) then
      call move_alloc(final_films,films);call move_alloc(final_film_mass,film_mass)
      if(present(film_steps)) film_steps=int(2*total_steps)
    end if
    if(present(heat_release)) heat_release=heat
    if(present(reaction_source)) reaction_source=ir
    if(present(gas_particle_change)) gas_particle_change=exchange(:nv)
    if(present(frame_energy_change)) frame_energy_change=frame_change
    if(present(level_steps)) level_steps=counts
    boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
    if(present(split_created)) split_created=watermark-before_split
    minimum_theta=theta;last_id=watermark;particle_steps=int(2*total_steps);ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
      if(status/=MPI_SUCCESS) agreed=.false.
    end subroutine agree
    logical function safe_sum(a,b)
      real(dp), intent(in) :: a(:),b(:)
      safe_sum=all(ieee_is_finite(a)).and.all(ieee_is_finite(b))
      if(safe_sum) safe_sum=.not.any(a>huge(1._dp)-max(b,0._dp)).and..not.any(a< -huge(1._dp)-min(b,0._dp))
    end function safe_sum
    subroutine add_source_rate(destination,rate,factor,good)
      real(dp), intent(inout) :: destination(:,:,:,:)
      real(dp), intent(in) :: rate(:,:,:,:),factor
      logical, intent(out) :: good
      real(dp) :: term(size(rate,1))
      integer :: x,y,z
      good=ieee_is_finite(factor).and.all(ieee_is_finite(rate))
      if(.not.good) return
      if(abs(factor)>1) then
        good=all(abs(rate)<=huge(1._dp)/abs(factor))
        if(.not.good) return
      end if
      do z=1,size(rate,4)
      do y=1,size(rate,3)
      do x=1,size(rate,2)
        term=factor*rate(:,x,y,z)
        good=safe_sum(destination(:,x,y,z),term)
        if(.not.good) return
        destination(:,x,y,z)=destination(:,x,y,z)+term
      end do
      end do
      end do
    end subroutine add_source_rate
    ! Conservative source synchronization: replace provisional copy-supported
    ! leaf exchange by the actual physical half-kick exchange at its receiver.
    ! This is a conservation extension to the pinned same-level source transfer.
    subroutine synchronize_particle_sources()
      real(dp) :: correction(nv),primitive(reactive_nprim(size(species))),recovered,sound
      real(dp) :: local_reset,global_reset,delta,potential,cell_volume
      integer :: p,x,y,z,d(3)
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,source_correction,correction,comm,agreed)
      if(.not.agreed) return
      valid=safe_sum(exchange(:nv),correction)
      call agree();if(.not.agreed) return
      local_reset=0
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        d=tree%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(covered(tree,p,[x,y,z]).or.geometry(p)%volume_fraction(x,y,z)<=0) cycle
          if(all(source_correction%fields(p)%values(:,x,y,z)==0)) cycle
          valid=safe_sum(gas%fields(p)%values(:,x,y,z),source_correction%fields(p)%values(:,x,y,z))
          if(.not.valid) exit
          gas%fields(p)%values(:,x,y,z)=gas%fields(p)%values(:,x,y,z)+source_correction%fields(p)%values(:,x,y,z)
          if(allocated(correction_prefix)) then
            potential=0
            if(present(frame)) then
              call rotating_potential_3d(frame, &
                [geometry(p)%x_lower+(x-.5_dp)*geometry(p)%dx, &
                 geometry(p)%y_lower+(y-.5_dp)*geometry(p)%dy,geometry(p)%z_lower+(z-.5_dp)*geometry(p)%dz], &
                potential,valid)
              if(.not.valid) exit
            endif
            delta=0
            call reset_nasa7_dual_energy_cell(species,gas%fields(p)%values(:,x,y,z), &
              stored%fields(p)%values(1,x,y,z),temp%fields(p)%values(1,x,y,z),held_policy,delta,valid,potential)
            if(.not.valid) exit
            cell_volume=geometry(p)%volume_fraction(x,y,z)*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
            if(cell_volume>1) then
              valid=abs(delta)<=huge(1._dp)/cell_volume
              if(.not.valid) exit
            endif
            delta=cell_volume*delta
            valid=safe_sum([local_reset],[delta]);if(.not.valid) exit
            local_reset=local_reset+delta
          endif
          if(allocated(stored)) then
            call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],gas%fields(p)%values(:,x,y,z), &
              temp%fields(p)%values(1,x,y,z),primitive,recovered,sound,valid,frame,stored%fields(p)%values(1,x,y,z))
          else
          call recover_reactive_eb_cell_3d(species,geometry(p),[x,y,z],gas%fields(p)%values(:,x,y,z), &
            temp%fields(p)%values(1,x,y,z),primitive,recovered,sound,valid,frame)
          endif
          if(.not.valid) exit
          temp%fields(p)%values(1,x,y,z)=recovered
        end do
        if(.not.valid) exit
        end do
        if(.not.valid) exit
        end do
        if(.not.valid) exit
      end do
      call agree();if(.not.agreed) return
      if(allocated(correction_prefix)) then
        call MPI_Allreduce(local_reset,global_reset,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
        if(status/=MPI_SUCCESS) then
          agreed=.false.;return
        endif
        valid=safe_sum([correction_prefix],[global_reset]);call agree();if(.not.agreed) return
        correction_prefix=correction_prefix+global_reset
      endif
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,gas,temp,comm,agreed,frame,internal_energy=stored)
      if(.not.agreed) return
      exchange(:nv)=exchange(:nv)+correction
    end subroutine synchronize_particle_sources
    subroutine particle_post_level(depth,time,interval,iteration,ncycle,post_comm,good)
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(MPI_Comm), intent(in) :: post_comm
      logical, intent(out) :: good
      logical :: synchronized
      integer :: p,status
      synchronized=.true.
      do p=depth+1,levels
        synchronized=synchronized.and.work(p)%next_stage==0.and. &
          abs(work(p)%time+work(p)%interval-time-interval)<=32*epsilon(dt)*max(dt,abs(time+interval))
      end do
      call MPI_Allreduce(synchronized,good,1,MPI_LOGICAL,MPI_LAND,post_comm,status)
      if(status/=MPI_SUCCESS) good=.false.
      if(.not.good) return
      if(depth==0.or.(iteration<ncycle.and.depth<levels-1)) then
        call redistribute_mpi_amr_subtree_particles(post_comm,tree,owners,physical,lower,spacing, &
          periodic,depth,iteration,options%maximum_parcels,good,phi)
        if(.not.good) return
      end if
      if(present(films)) then
        call finish_mpi_eb_film_level(post_comm,depth,time,interval,iteration,ncycle, &
          film_virtuals,film_ghosts,film_work,good)
        if(.not.good) return
      end if
      virtuals(depth+1)%parcels=virtuals(depth+1)%parcels(:0)
      virtuals(depth+1)%component_mass=virtuals(depth+1)%component_mass(:,:0)
      virtuals(depth+1)%nodes=virtuals(depth+1)%nodes(:0)
      if(iteration==ncycle) then
        ghosts(depth+1)%parcels=ghosts(depth+1)%parcels(:0)
        ghosts(depth+1)%component_mass=ghosts(depth+1)%component_mass(:,:0)
        ghosts(depth+1)%nodes=ghosts(depth+1)%nodes(:0)
      end if
    end subroutine particle_post_level
    subroutine unpack_particle_stage(packed,gas,held,good)
      type(amr_tree_state_3d), intent(in) :: packed
      type(amr_tree_state_3d), intent(out) :: gas
      type(amr_tree_state_3d), allocatable, intent(out) :: held
      logical, intent(out) :: good
      integer :: p,d(3),allocation_status
      if(.not.allocated(stored)) then
        gas=packed;good=.true.;return
      endif
      local_block: block
        logical :: local_ok
        local_ok=allocated(packed%fields)
        if(local_ok) local_ok=size(packed%fields)==size(owners).and.lbound(packed%fields,1)==1
        valid=local_ok;call agree();good=agreed;if(.not.good) return
        allocate(gas%fields(size(owners)),held,stat=allocation_status)
        valid=allocation_status==0;call agree();good=agreed;if(.not.good) return
        allocate(held%fields(size(owners)),stat=allocation_status)
        valid=allocation_status==0;call agree();good=agreed;if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          if(.not.allocated(packed%fields(p)%values)) then
            local_ok=.false.;exit
          endif
          d=tree%nodes(p)%shape
          local_ok=all(shape(packed%fields(p)%values)==[nv+1,d]).and.all(lbound(packed%fields(p)%values)==1)
          if(.not.local_ok) exit
          allocate(gas%fields(p)%values(nv,d(1),d(2),d(3)),held%fields(p)%values(1,d(1),d(2),d(3)), &
            stat=allocation_status)
          if(allocation_status/=0) then
            local_ok=.false.;exit
          endif
          gas%fields(p)%values=packed%fields(p)%values(:nv,:,:,:)
          held%fields(p)%values(1,:,:,:)=packed%fields(p)%values(nv+1,:,:,:)
        enddo
        valid=local_ok;call agree();good=agreed
      end block local_block
    end subroutine unpack_particle_stage
    subroutine particle_source(stage_gas,stage_temp,interval,stage,depth,time,rate,source_comm,good, &
        ancestor_start,ancestor_end,temperature_start,temperature_end,fractions)
      type(amr_tree_state_3d), intent(in) :: stage_gas,stage_temp,ancestor_start,ancestor_end,temperature_start,temperature_end
      real(dp), intent(in) :: interval,time,fractions(:)
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: remainder,physical_rate,film_rate,film_remainder,film_actual
      real(dp) :: wall_delta(5),film_wall_delta(5),outside(nv),wall_full(nv),force_full(nv)
      type(spray_wall_film), allocatable :: impact_films(:)
      real(dp), allocatable :: impact_mass(:,:)
      type(amr_tree_state_3d) :: sampled_gas,sampled_start,sampled_end,level_rate
      type(amr_tree_state_3d), allocatable :: sampled_held,held_start,held_end
      integer :: p,old_films
      call unpack_particle_stage(stage_gas,sampled_gas,sampled_held,good)
      if(.not.good) return
      call unpack_particle_stage(ancestor_start,sampled_start,held_start,good)
      if(.not.good) return
      call unpack_particle_stage(ancestor_end,sampled_end,held_end,good)
      if(.not.good) return

      if(impacts) then
        ! The particle snapshot is for impact thickness and appended deposits;
        ! it must never replace the independently evolved film hierarchy.
        impact_films=film_physical(depth+1)%films;impact_mass=film_physical(depth+1)%component_mass
        old_films=size(impact_films)
      call advance_mpi_eb_recursive_population_stage(source_comm,species,transport,liquids,options,tree,owners,geometry, &
        sampled_gas,stage_temp,predictor_width,depth,2,options%maximum_parcels,physical,virtuals,ghosts,work,time,interval,stage, &
        particle_iterations, &
        boundary_faces,level_rate,remainder,wall_delta,outside,good,phi,centroids,metrics, &
        sampled_start,sampled_end,temperature_start,temperature_end,fractions, &
        inflow_state,inflow_mach,inflow_temperature,inflow_pressure,physical_source=physical_rate, &
        acceleration=particle_acceleration,last_id=watermark,wall_films=impact_films,film_mass=impact_mass, &
        film_settings=film_settings,model_rejected=rejected,frame=frame,common_step=common, &
        internal_energy=sampled_held,ancestor_internal_energy_start=held_start, &
        ancestor_internal_energy_end=held_end)
        if(rejected.and.present(model_rejected)) model_rejected=.true.
      else
      call advance_mpi_eb_recursive_population_stage(source_comm,species,transport,liquids,options,tree,owners,geometry, &
        sampled_gas,stage_temp,predictor_width,depth,2,options%maximum_parcels,physical,virtuals,ghosts,work,time,interval,stage, &
        particle_iterations, &
        boundary_faces,level_rate,remainder,wall_delta,outside,good,phi,centroids,metrics, &
        sampled_start,sampled_end,temperature_start,temperature_end,fractions, &
        inflow_state,inflow_mach,inflow_temperature,inflow_pressure,physical_source=physical_rate, &
        acceleration=particle_acceleration,last_id=watermark,frame=frame,common_step=common, &
        internal_energy=sampled_held,ancestor_internal_energy_start=held_start, &
        ancestor_internal_energy_end=held_end)
      end if

      if(.not.good) return
      if(present(films)) then
        call advance_mpi_eb_recursive_film_stage(source_comm,species,transport,liquids,options,film_settings, &
          tree,owners,geometry,sampled_gas,stage_temp,predictor_width,depth,2,film_settings%maximum_films, &
          film_physical,film_virtuals,film_ghosts,film_work,time,interval,stage,particle_iterations,boundary_faces, &
          film_rate,film_remainder,film_wall_delta,good,phi,centroids,metrics, &
          sampled_start,sampled_end,temperature_start,temperature_end,fractions, &
          inflow_state,inflow_mach,inflow_temperature,inflow_pressure,physical_source=film_actual,frame=frame,common_step=common, &
        internal_energy=sampled_held,ancestor_internal_energy_start=held_start, &
        ancestor_internal_energy_end=held_end)

        if(.not.good) return
        if(impacts.and.stage==0) then
          ! The impact kernel preserves the old prefix and appends births.
          ! They miss the old film kick and join only the final half kick.
          call append_mpi_eb_film_births(source_comm,tree,owners,lower,spacing,2,film_settings%maximum_films,phi, &
            work(depth+1)%midpoint_films(old_films+1:),work(depth+1)%midpoint_film_mass(:,old_films+1:), &
            film_work(depth+1),good,periodic)
          if(.not.good) return
        end if
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          call add_source_rate(level_rate%fields(p)%values,film_rate%fields(p)%values,1._dp,valid)
          if(.not.valid) exit
        end do
        if(valid) valid=safe_sum(wall_delta,film_wall_delta)
        call agree();good=agreed
        if(.not.good) return
        wall_delta=wall_delta+film_wall_delta
      end if
      do p=1,size(owners)
        if(owners(p)/=rank.or.tree%nodes(p)%level/=depth) cycle
        if(stage==0) particle_first%fields(p)%values=level_rate%fields(p)%values
        particle_last%fields(p)%values=level_rate%fields(p)%values
      end do
      if(allocated(stored)) then
        call pack_eb_held_auxiliary_state_3d(level_rate,zero_held,rate,valid)
        call agree();good=agreed;if(.not.good) return
      else
        rate=level_rate
      endif

      if(stage/=particle_iterations) return
      ! All fields are rates. Count each accepted physical half kick once;
      ! remove only this logical level's provisional valid-cell source.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        associate(delta=>source_correction%fields(p)%values,actual=>physical_rate%fields(p)%values, &
          old_actual=>work(depth+1)%old_physical_source%fields(p)%values)
          call add_source_rate(delta,old_actual,.5_dp*interval,valid)
          if(valid) call add_source_rate(delta,actual,.5_dp*interval,valid)
          if(valid.and.tree%nodes(p)%level==depth) then
            call add_source_rate(delta,particle_first%fields(p)%values,-.5_dp*interval,valid)
            ! Later reacting MOL passes omit the new auxiliary family.
            ! Keep the already completed physical MK exchange in the accepted
            ! particle/gas inventory using the existing source correction.
            if(valid.and..not.(mol.and.chemistry.and.iterations>1)) &
              call add_source_rate(delta,level_rate%fields(p)%values,-.5_dp*interval,valid)
          end if
          if(valid.and.present(films)) then
            call add_source_rate(delta,film_work(depth+1)%old_physical_source%fields(p)%values,.5_dp*interval,valid)
            if(valid) call add_source_rate(delta,film_actual%fields(p)%values,.5_dp*interval,valid)
          end if
        end associate
        if(.not.valid) exit
      end do
      call agree();good=agreed
      if(.not.good) return
      wall_full=0;wall_full(:5)=wall_delta
      force_full=0;force_full(:5)=work(depth+1)%force_change
      valid=safe_sum(particle_wall,wall_full).and.safe_sum(exterior,outside).and. &
        safe_sum(particle_force_total,force_full)
      call agree();good=agreed
      if(.not.good) return
      particle_wall=particle_wall+wall_full;exterior=exterior+outside
      particle_force_total=particle_force_total+force_full
      frame_source=frame_source+work(depth+1)%frame_source_change
      frame_boundary=frame_boundary+work(depth+1)%frame_boundary_change
      if(present(films)) frame_source=frame_source+film_work(depth+1)%frame_source_change
    end subroutine particle_source
    ! All logical levels share stage time and interval; old stages precede
    ! every predictor and new stage. The shared level kernels retain separate
    ! physical, virtual and ghost populations and deposited halo remainders.
    subroutine common_particle_source(stage_gas,stage_temp,interval,stage,depth,rate,source_comm,good)
      type(amr_tree_state_3d), intent(in) :: stage_gas,stage_temp
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: level_rate
      real(dp) :: fractions(levels-1)
      integer :: level,p
      good=.false.
      if(depth/=-1) return
      if(allocated(stored)) then
        call pack_eb_held_auxiliary_state_3d(particle_first,zero_held,rate,valid)
        call agree();good=agreed;if(.not.good) return
      else
      rate=particle_first
      endif
      do p=1,size(owners)
        if(owners(p)==rank) rate%fields(p)%values=0
      end do
      fractions=0
      do level=0,levels-1
        call particle_source(stage_gas,stage_temp,interval,stage,level,0._dp,level_rate,source_comm,good, &
          stage_gas,stage_gas,stage_temp,stage_temp,fractions)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=level) cycle
          rate%fields(p)%values=level_rate%fields(p)%values
        end do
      end do
      good=.true.
    end subroutine common_particle_source
    subroutine common_particle_predictor(rate,interval,iteration,depth,halo,source_comm,good)
      type(amr_tree_state_3d), intent(in) :: rate
      real(dp), intent(in) :: interval
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: candidate,level_halo,combined_old,combined_remainder
      real(dp) :: fractions(levels-1)
      integer :: level,p,ios
      good=.false.
      if(depth/=-1) return
      allocate(candidate%fields(size(owners)),stat=ios)
      valid=ios==0;call agree();if(.not.agreed) return
      fractions=0
      do level=1,levels
        valid=work(level)%time==0.and.work(level)%interval==interval.and.work(level)%next_stage==iteration
        call agree();if(.not.agreed) return
        combined_old=work(level)%old_source;combined_remainder=work(level)%old_remainder
        if(present(films)) then
          valid=film_work(level)%time==0.and.film_work(level)%interval==interval.and. &
            film_work(level)%next_stage==iteration
          call agree();if(.not.agreed) return
          do p=1,size(owners)
            if(owners(p)/=rank) cycle
            call add_source_rate(combined_old%fields(p)%values,film_work(level)%old_source%fields(p)%values,1._dp,valid)
            if(valid) call add_source_rate(combined_remainder%fields(p)%values, &
              film_work(level)%old_remainder%fields(p)%values,1._dp,valid)
            if(.not.valid) exit
          end do
          call agree();if(.not.agreed) return
        end if
        call fill_sdc_mpi_eb_particle_predictor_halo_3d(tree,owners,geometry,rate,combined_old, &
          combined_remainder,nv,predictor_width,level_halo,source_comm,good,boundary_faces,metrics, &
          rate,rate,particle_first,particle_first,fractions)
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank.or.tree%nodes(p)%level/=level-1) cycle
          call move_alloc(level_halo%fields(p)%values,candidate%fields(p)%values)
        end do
      end do
      call move_alloc(candidate%fields,halo%fields);good=.true.
    end subroutine common_particle_predictor
    subroutine particle_predictor(rate,interval,iteration,depth,time,halo,source_comm,good, &
        source_start,source_end,fractions)
      type(amr_tree_state_3d), intent(in) :: rate,source_start,source_end
      real(dp), intent(in) :: interval,time,fractions(:)
      integer, intent(in) :: iteration,depth
      type(amr_tree_state_3d), intent(inout) :: halo
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: combined_old,combined_remainder
      integer :: p
      valid=work(depth+1)%time==time.and.work(depth+1)%interval==interval.and.work(depth+1)%next_stage==iteration
      call agree();good=agreed
      if(.not.good) return
      combined_old=work(depth+1)%old_source;combined_remainder=work(depth+1)%old_remainder
      if(present(films)) then
        valid=film_work(depth+1)%time==time.and.film_work(depth+1)%interval==interval.and. &
          film_work(depth+1)%next_stage==iteration
        call agree();good=agreed
        if(.not.good) return
        do p=1,size(owners)
          if(owners(p)/=rank) cycle
          call add_source_rate(combined_old%fields(p)%values,film_work(depth+1)%old_source%fields(p)%values,1._dp,valid)
          if(valid) call add_source_rate(combined_remainder%fields(p)%values, &
            film_work(depth+1)%old_remainder%fields(p)%values,1._dp,valid)
          if(.not.valid) exit
        end do
        call agree();good=agreed
        if(.not.good) return
      end if
      call fill_sdc_mpi_eb_particle_predictor_halo_3d(tree,owners,geometry,rate,combined_old, &
        combined_remainder,nv,predictor_width,halo,source_comm,good,boundary_faces,metrics, &
        source_start,source_end,particle_first,particle_last,fractions)
    end subroutine particle_predictor
  end subroutine advance_mpi_eb_spray_sdc_subcycle

  ! S/2 - F/2 - EB reaction/transport/hydro - F/2 - S/2; atomic publication.
  ! Supplying last_id with film state enables EB impacts in both particle halves.
  ! Replicated tables/options and a consistent geometry/nodal field are caller contracts.
  ! Ledgers use change-in-system signs, including the particle wall impulse.
  subroutine advance_mpi_eb_spray_reactive(comm,tree,owners,geometry,phi,species,reactions,transport, &
      liquids,options,state,temperature,parcels,component_mass,lower,spacing,dt,solver, &
      chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      particle_steps,minimum_theta,boundary_change,wall_change,chemical_change,ok,integrator,subcycling, &
      films,film_mass,film_settings,film_steps,last_id,model_rejected,split_created,particle_outflow, &
      acceleration,source_change,inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry, &
      boundary_faces,inflow_state,wall_temperature)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    character(len=*), intent(in) :: solver
    logical, intent(in) :: chemistry,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    integer, intent(out) :: particle_steps
    real(dp), intent(inout) :: minimum_theta,boundary_change(:),wall_change(:),chemical_change(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient,low_coefficient,high_coefficient
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),wall_temperature(6)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    integer, optional, intent(out) :: film_steps
    integer(int64), optional, intent(inout) :: last_id
    logical, optional, intent(out) :: model_rejected
    integer(int64), optional, intent(out) :: split_created
    real(dp), optional, intent(in) :: acceleration(3)
    real(dp), optional, intent(inout) :: source_change(:)
    real(dp) :: accel(3),low_accel(3),high_accel(3),force_ledger(reactive_nvar(size(species)))
    logical, optional, intent(in) :: particle_outflow
    character(len=32), allocatable :: particle_faces(:)
    logical :: outflow
    real(dp), allocatable :: half_exterior(:)
    real(dp) :: exterior(reactive_nvar(size(species)))
    integer(int64) :: before_split
    real(dp) :: low_multiplicity,high_multiplicity
    integer(int64) :: watermark
    type(spray_wall_film), allocatable :: patches(:)
    real(dp), allocatable :: patch_mass(:,:)
    real(dp) :: first_film(reactive_nvar(size(species))),second_film(reactive_nvar(size(species)))
    integer :: flags,lowflags,highflags,first_film_steps,second_film_steps
    character(len=32) :: film_faces(6)
    type(amr_tree_state_3d) :: gas,temp
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: masses(:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    real(dp) :: chemical(reactive_nvar(size(species))),first(reactive_nvar(size(species)))
    real(dp) :: second(reactive_nvar(size(species))),theta
    integer :: first_steps,second_steps,ierr,nv
    logical :: valid,agreed
    ok=.false.;particle_steps=0;nv=size(boundary)
    if(present(film_steps)) film_steps=0
    if(present(split_created)) split_created=0
    if(present(model_rejected)) model_rejected=.false.
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(coefficient).and.coefficient>=0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(coefficient,low_coefficient,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,high_coefficient,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_coefficient/=high_coefficient) return
    accel=0;force_ledger=0
    if(present(acceleration)) accel=acceleration
    valid=all(ieee_is_finite(accel))
    valid=valid.and.all(abs(accel)<=1e6_dp)
    if(any(accel/=0)) valid=valid.and.present(source_change)
    if(present(source_change)) valid=valid.and.size(source_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(accel,low_accel,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(accel,high_accel,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(low_accel/=high_accel)) return
    outflow=.false.;exterior=0
    if(present(particle_outflow)) outflow=particle_outflow
    valid=.true.
    if(present(boundary_faces)) then
      valid=all(boundary_faces=='inflow'.or.boundary_faces=='outflow'.or.boundary_faces=='periodic'.or. &
        boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or.boundary_faces=='isothermal_wall')
      if(any(boundary_faces=='periodic'.or.boundary_faces=='slip_wall'.or. &
        boundary_faces=='no_slip_wall'.or.boundary_faces=='isothermal_wall')) then
        particle_faces=boundary_faces
        where(particle_faces=='inflow') particle_faces='outflow'
        where(particle_faces=='isothermal_wall') particle_faces='no_slip_wall'
      end if
      valid=valid.and.outflow
    else
      valid=.not.present(inflow_state)
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    watermark=0
    if(present(last_id)) watermark=last_id
    if(options%etab_parcel_factor>=0) before_split=watermark
    flags=merge(1,0,present(films))+2*merge(1,0,present(film_mass))+4*merge(1,0,present(film_settings))+ &
      8*merge(1,0,present(last_id))+16*merge(1,0,outflow)
    call MPI_Allreduce(flags,lowflags,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,highflags,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.lowflags/=highflags) return
    if(outflow) then
      if(flags/=16.and.flags/=24.and.flags/=31) return
      allocate(half_exterior(nv));half_exterior=0
      flags=flags-16
    end if
    if(flags/=0.and.flags/=7.and.flags/=8.and.flags/=15) return
    valid=valid_spray_options(options)
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.(options%maximum_multiplicity==0.or.present(last_id))
    if(present(films)) valid=valid.and.allocated(films).and.allocated(film_mass)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(options%maximum_multiplicity,low_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(options%maximum_multiplicity,high_multiplicity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.low_multiplicity/=high_multiplicity) return
    if(present(films)) then
      patches=films;patch_mass=film_mass
    end if
    first_film=0;second_film=0;first_film_steps=0;second_film_steps=0;film_faces='slip_wall'
    if(present(boundary_faces)) then
      where(boundary_faces=='periodic') film_faces='periodic'
    end if
    valid=allocated(parcels).and.allocated(component_mass).and.ieee_is_finite(dt)
    if(valid) valid=dt>0.and.dt<=1e6_dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.size(chemical_change)==nv
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    boundary=0;wall=0;chemical=0;first=0;second=0;theta=1
    if(any(accel/=0)) then
      call force_kick()
      if(.not.valid) return
    end if
    call particle_half(first_steps,first)
    if(.not.valid) return
    if(present(films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,film_faces,first_film_steps,first_film,valid,geometry,phi)
      if(.not.valid) return
    end if
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,tree,owners,geometry,solver,dt, &
      chemistry,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      gas,temp,boundary,wall,chemical,theta,comm,valid,integrator,subcycling, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      difmag=coefficient,extended_geometry=extended_geometry,boundary_faces=boundary_faces, &
      inflow_state=inflow_state,wall_temperature=wall_temperature)
    if(.not.valid) return
    if(present(films)) then
      call advance_mpi_amr_film_coupling(comm,species,transport,liquids,options,film_settings,gas,temp, &
        tree,owners,patches,patch_mass,lower,spacing,.5_dp*dt,film_faces,second_film_steps,second_film,valid,geometry,phi)
      if(.not.valid) return
    end if
    call particle_half(second_steps,second)
    if(.not.valid) return
    if(options%etab_parcel_factor<0) before_split=watermark
    if(options%maximum_multiplicity>0.and.options%etab_parcel_factor<0) then
      call split_mpi_spray_parcels(comm,drops,options%maximum_multiplicity,options%maximum_parcels,watermark,valid,masses)
      if(.not.valid) return
    end if
    if(any(accel/=0)) then
      call force_kick()
      if(.not.valid) return
    end if
    state=gas;temperature=temp
    call move_alloc(drops,parcels)
    call move_alloc(masses,component_mass)
    if(present(films)) then
      call move_alloc(patches,films);call move_alloc(patch_mass,film_mass)
    end if
    if(present(last_id)) last_id=watermark
    if(present(split_created)) split_created=watermark-before_split
    if(present(film_steps)) film_steps=first_film_steps+second_film_steps
    boundary_change=boundary+exterior;wall_change=wall+first+second+first_film+second_film;chemical_change=chemical
    if(present(source_change)) source_change=force_ledger
    minimum_theta=theta;particle_steps=first_steps+second_steps;ok=.true.
  contains
    subroutine force_kick()
      real(dp) :: local(nv),global(nv)
      call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,gas,nv,comm,valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temp,1,comm,valid)
      if(.not.valid) return
      call local_force_kick(local)
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      valid=ierr==MPI_SUCCESS.and.agreed
      if(.not.valid) return
      call MPI_Allreduce(local,global,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      valid=ierr==MPI_SUCCESS
      if(valid) force_ledger=force_ledger+global
    end subroutine
    subroutine local_force_kick(local)
      real(dp), intent(out) :: local(nv)
      real(dp) :: dv(3),mass,momentum(3),energy,volume,velocity(3),prim(nv),checked_t,sound
      integer :: rank,status,node,i,j,k,p,n(3)
      valid=.false.;local=0;dv=.5_dp*dt*accel
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          mass=gas%fields(node)%values(1,i,j,k)
          if(.not.ieee_is_finite(mass).or.mass<=0) return
          velocity=gas%fields(node)%values(2:4,i,j,k)/mass
          if(any(abs(velocity+dv)>1e7_dp)) return
          momentum=mass*dv
          energy=dot_product(gas%fields(node)%values(2:4,i,j,k)+.5_dp*momentum,dv)
          gas%fields(node)%values(2:4,i,j,k)=gas%fields(node)%values(2:4,i,j,k)+momentum
          gas%fields(node)%values(5,i,j,k)=gas%fields(node)%values(5,i,j,k)+energy
          call reactive_conserved_to_primitive(species,gas%fields(node)%values(:,i,j,k), &
            temp%fields(node)%values(1,i,j,k),prim,checked_t,sound,valid)
          if(.not.valid) return
          valid=.false.
          if(abs(checked_t-temp%fields(node)%values(1,i,j,k))> &
            1e-9_dp*max(1._dp,abs(temp%fields(node)%values(1,i,j,k)))) return
          if(covered(tree,node,[i,j,k])) cycle
          volume=geometry(node)%volume_fraction(i,j,k)*geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
          local(2:4)=local(2:4)+volume*momentum;local(5)=local(5)+volume*energy
        end do
        end do
        end do
      end do
      do p=1,size(drops)
        if(drops(p)%mass<=0) cycle
        mass=drops(p)%mass*drops(p)%multiplicity;velocity=drops(p)%velocity
        drops(p)%velocity=velocity+dv
        if(.not.valid_spray_parcel(drops(p))) return
        local(2:4)=local(2:4)+mass*dv
        local(5)=local(5)+mass*dot_product(velocity+.5_dp*dv,dv)
      end do
      valid=all(ieee_is_finite(local))
    end subroutine
    subroutine particle_half(steps,exchange)
      integer, intent(out) :: steps
      real(dp), intent(inout) :: exchange(:)
      logical :: rejected
      if(present(last_id).and.present(films)) then
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,last_id=watermark,boundary_change=exchange, &
          wall_films=patches,film_mass=patch_mass,film_settings=film_settings, &
          eb_geometry=geometry,nodal_phi=phi,model_rejected=rejected,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
        if(present(model_rejected)) model_rejected=rejected
      else if(present(last_id)) then
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,last_id=watermark,boundary_change=exchange, &
          eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
      else
        call advance_mpi_amr_mixture_coupling(comm,species,transport,liquids,options,gas,temp,tree,owners, &
          drops,masses,lower,spacing,.5_dp*dt,steps,valid,boundary_change=exchange, &
          eb_geometry=geometry,nodal_phi=phi,eb_outflow_change=half_exterior,boundary_faces=particle_faces)
      end if
      if(valid.and.outflow) exterior=exterior+half_exterior
    end subroutine
  end subroutine
  ! Common-step gas/particle SDC. New kicks are reevaluated from one saved
  ! kick-drift endpoint, never cumulatively advanced across SDC iterations.
  ! Existing stationary films evolve thermally; particle contacts remain specular.
  subroutine advance_mpi_eb_spray_sdc(comm,tree,owners,geometry,extended_geometry,phi,species,reactions,transport, &
      liquids,options,state,temperature,parcels,component_mass,lower,spacing,dt,iterations,chemistry,transport_flags, &
      force,rtol,atol,boundary_faces,boundary_change,wall_change,chemical_change,source_change,minimum_theta, &
      particle_steps,last_id,ok,reaction_source,difmag,redistribution,target_volume_fraction,gas_particle_change, &
      inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,split_created,particle_acceleration, &
      films,film_mass,film_settings,film_steps,film_impacts,model_rejected,gas_acceleration,evaluate_diffusion,evaluate_reaction, &
      embedded_wall_temperature,embedded_noslip,mol_mode,mol_lax_friedrichs,frame, &
      state_redist_max_order,evaluate_spatial_diffusion,spatial_diffusion_halo_width, &
      constant_coefficients,heat_release,frame_energy_change, &
      internal_energy,energy_policy,energy_correction,internal_energy_pressure_work_change)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction,internal_energy_pressure_work_change
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(inout) :: frame_energy_change
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(mpi_eb_sdc_diffusion_source), optional :: evaluate_diffusion
    procedure(mpi_eb_sdc_spatial_diffusion_source), optional :: evaluate_spatial_diffusion
    integer, optional, intent(in) :: spatial_diffusion_halo_width
    type(spray_wall_film), allocatable, optional, intent(inout) :: films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    logical, optional, intent(in) :: film_impacts
    logical, optional, intent(out) :: model_rejected
    integer, optional, intent(out) :: film_steps
    real(dp), optional, intent(in) :: gas_acceleration(3)
    type(MPI_Comm), intent(in) :: comm
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),iterations
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    real(dp), intent(in) :: phi(0:,0:,0:),lower(3),spacing(3),dt,force(3),rtol,atol
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquids(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    logical, intent(in) :: chemistry,transport_flags(4)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, intent(out) :: particle_steps
    integer(int64), intent(inout) :: last_id
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release
    real(dp), optional, intent(in) :: difmag,target_volume_fraction
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(inout) :: gas_particle_change(:)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: embedded_noslip,mol_mode,mol_lax_friedrichs
    real(dp), optional, intent(in) :: wall_temperature(6)
    real(dp), optional, intent(in) :: particle_acceleration(3)
    integer(int64), optional, intent(out) :: split_created
    call advance_mpi_eb_spray_sdc_subcycle(comm,tree,owners,geometry,extended_geometry,phi, &
      species,reactions,transport,liquids,options,state,temperature,parcels,component_mass,lower,spacing, &
      dt,iterations,chemistry,transport_flags,force,rtol,atol,boundary_faces,boundary_change,wall_change, &
      chemical_change,source_change,minimum_theta,particle_steps,last_id,ok,reaction_source,difmag, &
      redistribution,target_volume_fraction,gas_particle_change,gas_acceleration=gas_acceleration, &
      inflow_state=inflow_state,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature, &
      inflow_pressure=inflow_pressure,wall_temperature=wall_temperature,particle_acceleration=particle_acceleration, &
      split_created=split_created,films=films,film_mass=film_mass,film_settings=film_settings,film_steps=film_steps, &
      film_impacts=film_impacts,model_rejected=model_rejected,evaluate_diffusion=evaluate_diffusion, &
      evaluate_reaction=evaluate_reaction,embedded_wall_temperature=embedded_wall_temperature, &
      embedded_noslip=embedded_noslip,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,frame=frame, &
      state_redist_max_order=state_redist_max_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
      spatial_diffusion_halo_width=spatial_diffusion_halo_width,common_step=.true., &
      constant_coefficients=constant_coefficients,heat_release=heat_release,frame_energy_change=frame_energy_change, &
        internal_energy=internal_energy,energy_policy=energy_policy,energy_correction=energy_correction, &
        internal_energy_pressure_work_change=internal_energy_pressure_work_change)
  end subroutine advance_mpi_eb_spray_sdc
end module
