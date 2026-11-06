module mpi_eb_tree_transport_3d_mod
  use mpi_f08
  use mpi_eb_wall_heat_3d_mod, only: compute_mpi_eb_wall_transport_3d
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d, &
    compute_sdc_eb_transport_3d,compute_sdc_eb_transport_halo_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d, &
    validate_mpi_eb_tree_geometry_3d,integrate_mpi_eb_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_transport_source_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d
  use mpi_eb_tree_redistribution_3d_mod, only: compute_state_redist_mpi_level_stage_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: agree_mpi_constant_transport_coefficients
  implicit none
  private
  public :: advance_mpi_eb_tree_transport_3d,compute_sdc_mpi_eb_tree_transport_3d
  public :: compute_sdc_mpi_eb_tree_transport_halos_3d
contains
  ! Common molecular-transport interval with closed physical diffusion faces.
  subroutine advance_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,dt,viscosity,conduction, &
      diffusion,barodiffusion,state,temperature,boundary_change,minimum_theta,comm,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
      constant_coefficients)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    logical :: constant_agreed
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    real(dp) :: thermal(7),thermal_low(7),thermal_high(7)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,ghost,gt
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: work_stamp
    real(dp) :: boundary(size(boundary_change)),theta,local_theta,stage_theta,low,high
    integer :: nv,rank,status,node,n(3),c,ios,depth,flags(4),minflags(4),maxflags(4)
    logical :: faces(6)
    character(len=16) :: patch_faces(6)
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed,species=species)
    if(.not.constant_agreed) return
    valid=ieee_is_finite(dt)
    if (valid) valid=dt>0._dp
    valid=valid.and.size(boundary_change)==nv
    thermal=0
    if(present(wall_temperature)) then
      thermal(1)=1;thermal(2:)=wall_temperature
      valid=valid.and.present(boundary_faces)
      valid=valid.and.all(ieee_is_finite(wall_temperature)).and.all(wall_temperature>0)
    end if
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) valid=valid.and.conduction.and.present(wall_temperature)
    end if
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(thermal,thermal_low,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(thermal,thermal_high,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(thermal_low/=thermal_high)) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high) return
    flags=merge(1,0,[viscosity,conduction,diffusion,barodiffusion])
    call MPI_Allreduce(flags,minflags,4,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,maxflags,4,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(minflags/=maxflags)) return
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    candidate=snapshot;ct=st
    allocate(flux(size(owners)),stat=ios)
    valid=ios==0
    call agree(valid,ready)
    if (.not.ready) return
    local_theta=1._dp
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;depth=tree%nodes(node)%level+1
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      work_stamp=tree_work_start()
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
        patch_faces='outflow'
        if(present(boundary_faces)) then
          where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
            boundary_faces=='isothermal_wall')) &
            patch_faces=boundary_faces
          faces=.true.
        end if
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,geometry(node),dt, &
        viscosity,conduction,diffusion,barodiffusion,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        flux(node),stage_theta,valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),faces, &
          physical_faces=patch_faces,wall_temperature=wall_temperature, &
        constant_coefficients=constant_coefficients)
      call tree_work_stop(node,work_stamp)
      if (.not.valid) exit
      local_theta=min(local_theta,stage_theta)
    end do
    call agree(valid.and.ieee_is_finite(local_theta),ready)
    if (.not.ready) return
    call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,flux,comm,ready, &
      boundary_faces=boundary_faces)
    if (.not.ready) return
    boundary=0._dp
    if (rank==owners(1)) then
      n=tree%nodes(1)%shape
      associate(g=>geometry(1),f=>flux(1))
      do c=1,nv
        boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
          g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
          g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
          g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
          g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
          g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
      end do
      end associate
    end if
    call MPI_Bcast(boundary,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(local_theta,theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.ieee_is_finite(theta),ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;minimum_theta=theta;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
  ! Owner-local SDC rates and dt-integrated flux densities; no state update,
  ! restriction, redistribution or reflux. Publish only after collective success.
  subroutine compute_sdc_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,dt,viscosity,conduction, &
      diffusion,barodiffusion,state,temperature,transport_rate,integrated_flux,minimum_theta,comm,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      embedded_wall_temperature,extended_geometry,embedded_wall_energy,embedded_frequency, &
      embedded_noslip,embedded_wall_impulse,frame, &
      constant_coefficients,internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), allocatable :: held_ghost
    real(dp), allocatable :: patch_held(:,:,:),ghost_held(:,:,:)
    logical :: constant_agreed
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(amr_tree_state_3d), intent(inout) :: transport_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    real(dp), intent(inout) :: minimum_theta
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    integer :: selected_depth,min_depth,max_depth
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(inout) :: embedded_wall_energy,embedded_frequency
    logical, optional, intent(in) :: embedded_noslip
    real(dp), optional, intent(inout) :: embedded_wall_impulse(3)
    logical :: noslip,wall_active
    real(dp) :: wall_values(reactive_nvar(size(species))),wall_impulse(3),heat_frequency,heat_energy
    real(dp) :: thermal(10),thermal_low(10),thermal_high(10)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: rates,ghost,gt,heat
    real(dp), allocatable :: ordinary(:,:,:,:),ordinary_ghost(:,:,:,:)
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: work_stamp
    real(dp) :: theta,local_theta,stage_theta,low,high
    integer :: nv,rank,status,node,n(3),ios,depth,flags(4),minflags(4),maxflags(4)
    logical :: faces(6)
    character(len=16) :: patch_faces(6)
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed,species=species)
    if(.not.constant_agreed) return
    valid=ieee_is_finite(dt)
    if (valid) valid=dt>0._dp
    thermal=0;noslip=.false.
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=present(embedded_wall_temperature).or.noslip
    thermal(10)=merge(1._dp,0._dp,noslip)
    valid=valid.and.(noslip.eqv.present(embedded_wall_impulse))
    if(noslip) valid=valid.and.viscosity
    if(present(wall_temperature)) then
      thermal(1)=1;thermal(2:7)=wall_temperature
      valid=valid.and.present(boundary_faces)
      valid=valid.and.all(ieee_is_finite(wall_temperature)).and.all(wall_temperature>0)
    end if
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) valid=valid.and.conduction.and.present(wall_temperature)
    end if
    valid=valid.and.(wall_active.eqv.present(extended_geometry))
    valid=valid.and.(wall_active.eqv.present(embedded_wall_energy))
    valid=valid.and.(wall_active.eqv.present(embedded_frequency))
    if(present(embedded_wall_temperature)) then
      valid=valid.and.conduction.and.ieee_is_finite(embedded_wall_temperature).and.embedded_wall_temperature>0
      thermal(8:9)=[1._dp,embedded_wall_temperature]
    end if
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(thermal,thermal_low,10,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(thermal,thermal_high,10,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(thermal_low/=thermal_high)) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high) return
    flags=merge(1,0,[viscosity,conduction,diffusion,barodiffusion])
    call MPI_Allreduce(flags,minflags,4,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,maxflags,4,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(minflags/=maxflags)) return
    selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    valid=selected_depth>=-1
    if(present(target_depth)) valid=selected_depth>=0
    if(present(ancestor_start).or.present(ancestor_end).or.present(temperature_start).or. &
      present(temperature_end).or.present(ancestor_fraction)) valid=valid.and.selected_depth>=0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(selected_depth,min_depth,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selected_depth,max_depth,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.min_depth/=max_depth) return
    if(present(internal_energy)) allocate(held_ghost)
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,restrict_snapshot=.false., &
      ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
      temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,frame=frame, &
      internal_energy=internal_energy,ghost_internal_energy=held_ghost, &
      ancestor_internal_energy_start=ancestor_internal_energy_start, &
      ancestor_internal_energy_end=ancestor_internal_energy_end)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call agree(selected_depth<tree%level_count(),ready)
    if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    allocate(flux(size(owners)),rates%fields(size(owners)),stat=ios)
    valid=ios==0
    call agree(valid,ready)
    if (.not.ready) return
    local_theta=1._dp
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;depth=tree%nodes(node)%level+1
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      if(selected_depth>=0.and.tree%nodes(node)%level/=selected_depth) cycle
      work_stamp=tree_work_start()
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      patch_faces='outflow'
      if(present(boundary_faces)) then
        where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
          boundary_faces=='isothermal_wall')) &
          patch_faces=boundary_faces
        faces=.true.
      end if
      if(present(internal_energy)) then
        patch_held=internal_energy%fields(node)%values(1,:,:,:)
        ghost_held=held_ghost%fields(node)%values(1,:,:,:)
      endif
      ordinary=state%fields(node)%values
      ordinary_ghost=ghost%fields(node)%values
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,geometry(node),[1,1,1], &
        geometry(node)%volume_fraction>0,ordinary,valid)
      endif
      if(.not.valid) exit
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,geometry(node),[0,0,0], &
        ordinary_ghost(1,:,:,:)>0,ordinary_ghost,valid)
      endif
      if(.not.valid) exit
      call compute_sdc_eb_transport_3d(species,transport,ordinary, &
        temperature%fields(node)%values(1,:,:,:),geometry(node),dt, &
        [viscosity,conduction,diffusion,barodiffusion],rates%fields(node)%values, &
        flux(node),stage_theta,valid,ordinary_ghost,gt%fields(node)%values(1,:,:,:),faces, &
        physical_faces=patch_faces,wall_temperature=wall_temperature, &
        constant_coefficients=constant_coefficients,internal_energy_density=patch_held, &
        ghost_internal_energy_density=ghost_held)
      call tree_work_stop(node,work_stamp)
      if (.not.valid) exit
      local_theta=min(local_theta,stage_theta)
    end do
    call agree(valid.and.ieee_is_finite(local_theta),ready)
    if (.not.ready) return
    call MPI_Allreduce(local_theta,theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS.and.ieee_is_finite(theta),ready)
    if (.not.ready) return
    if(wall_active) then
      call compute_mpi_eb_wall_transport_3d(species,transport,tree,owners,geometry,extended_geometry,state,temperature, &
        present(embedded_wall_temperature),noslip,thermal(9),heat,wall_values,heat_frequency,comm,ready, &
        boundary_faces=boundary_faces, &
        inflow_state=inflow_state,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature, &
        inflow_pressure=inflow_pressure,target_depth=target_depth,ancestor_start=ancestor_start, &
        ancestor_end=ancestor_end,temperature_start=temperature_start,temperature_end=temperature_end, &
        ancestor_fraction=ancestor_fraction,frame=frame, &
        constant_coefficients=constant_coefficients,internal_energy=internal_energy, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      if(.not.ready) return
      heat_energy=dt*wall_values(5);wall_impulse=dt*wall_values(2:4)
      valid=all(ieee_is_finite([heat_energy,wall_impulse]))
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        rates%fields(node)%values=rates%fields(node)%values+heat%fields(node)%values
        valid=valid.and.all(ieee_is_finite(rates%fields(node)%values))
      end do
      call agree(valid,ready)
      if(.not.ready) return
    end if
    call move_alloc(rates%fields,transport_rate%fields)
    call move_alloc(flux,integrated_flux)
    minimum_theta=theta;ok=.true.
    if(wall_active) then
      embedded_wall_energy=heat_energy;embedded_frequency=heat_frequency
      if(noslip) embedded_wall_impulse=wall_impulse
    end if
  contains
    subroutine ordinary_transport_energy(rotation,grid,first,active,field,good)
      type(rotating_frame_3d), optional, intent(in) :: rotation
      type(eb_geometry_3d), intent(in) :: grid
      integer, intent(in) :: first(3)
      logical, intent(in) :: active(:,:,:)
      real(dp), intent(inout) :: field(:,:,:,:)
      logical, intent(out) :: good
      real(dp) :: position(3),energy
      integer :: x,y,z
      good=.true.
      if(.not.present(rotation)) return
      good=.false.
      if(.not.rotation%is_valid().or.any(shape(active)/=shape(field(1,:,:,:)))) return
      do z=1,size(active,3)
      do y=1,size(active,2)
      do x=1,size(active,1)
        if(.not.active(x,y,z)) cycle
        position=[grid%x_lower,grid%y_lower,grid%z_lower]+ &
          (real(first+[x,y,z]-1,dp)-.5_dp)*[grid%dx,grid%dy,grid%dz]
        energy=field(5,x,y,z)
        call convert_rotating_energy_3d(rotation,position,field(1,x,y,z),field(5,x,y,z),.false.,energy,good)
        if(.not.good) return
        field(5,x,y,z)=energy
      end do
      end do
      end do
      good=.true.
    end subroutine ordinary_transport_energy
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine

  ! Prepared same-resolution transport RHS before per-level StateRedist.
  ! Only owner-local fields are stored; publication is collectively transactional.
  subroutine compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,tree,owners,geometry,extended_geometry, &
      state,temperature,dt,flags,rate,flux,source_halos,old_halos,minimum_theta,comm,ok, &
      boundary_faces,inflow_state,inflow_temperature,wall_temperature,target_depth,ancestor_start,ancestor_end, &
      temperature_start,temperature_end,ancestor_fraction,frame,embedded_wall_temperature,embedded_noslip, &
      embedded_exchange,embedded_frequency,inflow_mach,inflow_pressure, &
      registers,max_order,target_volume_fraction, &
      constant_coefficients,internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end,coarse_linear, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions,primitive_from_total_energy)
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical, optional, intent(in) :: coarse_linear
    ! Normal face transport may share the ordinary MOL total-energy primitive recovery.
    ! Held halos and wall-source recovery retain their existing contracts.
    logical, optional, intent(in) :: primitive_from_total_energy
    logical :: linear_halos,use_total_energy
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), allocatable :: held_ghost
    real(dp), allocatable :: patch_held(:,:,:),ghost_held(:,:,:)
    logical :: constant_agreed
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: registers(:)
    integer, optional, intent(in) :: max_order
    real(dp), optional, intent(in) :: target_volume_fraction
    type(state_redist_register_3d), allocatable :: pending_registers(:)
    integer :: first_depth,last_depth,stage_depth
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: dt
    logical, intent(in) :: flags(4)
    type(amr_tree_state_3d), intent(inout) :: rate,source_halos,old_halos
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: minimum_theta
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(6),wall_temperature(6), &
      inflow_mach(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: embedded_noslip
    real(dp), optional, intent(inout) :: embedded_exchange(:),embedded_frequency
    type(amr_tree_state_3d) :: rates,pending_halos,pending_old,metrics,ghost,gt,vf,body
    type(eb_integrated_flux_3d), allocatable :: pending_flux(:)
    type(eb_geometry_3d) :: center
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: ordinary(:,:,:,:),patch_halo(:,:,:,:),body_halo(:,:,:,:),valid_body(:,:,:,:)
    real(dp) :: controls(21),low_controls(21),high_controls(21),theta,patch_theta,frequency, &
      patch_frequency,global_theta,global_frequency,wall_t,stamp, &
      exchange(reactive_nvar(size(species))),ignored_exchange(reactive_nvar(size(species))),ignored_frequency
    integer :: nv,depth,rank,status,p,n(3),level,a,ios
    logical :: valid,ready,noslip,wall_active,faces(6)
    character(len=32) :: labels(6)
    ok=.false.;nv=reactive_nvar(size(species));depth=-1;theta=1;frequency=0;wall_t=0;noslip=.false.
    use_total_energy=.false.
    if(present(primitive_from_total_energy)) use_total_energy=primitive_from_total_energy
    linear_halos=.false.
    if(present(coarse_linear)) linear_halos=coarse_linear
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed,species=species)
    if(.not.constant_agreed) return
    if(present(target_depth)) depth=target_depth
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    valid=ieee_is_finite(dt).and.dt>0.and.size(transport)==size(species).and.depth>=-1
    if(present(target_depth)) valid=valid.and.depth>=0
    if(present(ancestor_start).or.present(ancestor_end).or.present(temperature_start).or. &
      present(temperature_end).or.present(ancestor_fraction)) valid=valid.and.depth>=0
    if(noslip) valid=valid.and.flags(1)
    valid=valid.and.(wall_active.eqv.present(embedded_exchange)).and.(wall_active.eqv.present(embedded_frequency))
    if(present(embedded_exchange)) valid=valid.and.size(embedded_exchange)==nv
    if(present(embedded_wall_temperature)) then
      wall_t=embedded_wall_temperature
      valid=valid.and.flags(2).and.ieee_is_finite(wall_t).and.wall_t>0
    endif
    valid=valid.and.(present(registers).or.(.not.present(max_order).and..not.present(target_volume_fraction)))
    if(present(target_volume_fraction)) valid=valid.and.ieee_is_finite(target_volume_fraction)
    valid=valid.and.(.not.linear_halos.or.present(internal_energy))
    valid=valid.and.(.not.use_total_energy.or..not.present(frame))
    if(present(separate_momentum)) valid=valid.and.(.not.use_total_energy.or..not.separate_momentum)
    controls=0;controls(21)=merge(1._dp,0._dp,use_total_energy);controls(19)=merge(1._dp,0._dp,linear_halos)
    if(present(separate_momentum)) then
      controls(20)=merge(1._dp,0._dp,separate_momentum)
      valid=valid.and.(.not.separate_momentum.or.linear_halos)
    endif
    if(present(consistent_velocity)) controls(20)=controls(20)+merge(2._dp,0._dp,consistent_velocity)
    if(present(paired_pressure)) controls(20)=controls(20)+merge(4._dp,0._dp,paired_pressure)
    if(present(connected_directions)) controls(20)=controls(20)+merge(8._dp,0._dp,connected_directions)
    valid=valid.and.iand(int(controls(20)),3)/=2
    if(btest(int(controls(20)),3)) valid=valid.and.btest(int(controls(20)),0)
    if(btest(int(controls(20)),2)) valid=valid.and.iand(int(controls(20)),7)==7
    controls(16)=merge(1._dp,0._dp,present(registers));controls(18)=.5_dp
    if(present(max_order)) controls(17)=real(max_order,dp)
    if(present(target_volume_fraction)) controls(18)=target_volume_fraction
    controls(1)=dt;controls(2:5)=real(merge(1,0,flags),dp);controls(6)=real(depth,dp)
    controls(7)=merge(1._dp,0._dp,noslip);controls(8)=merge(1._dp,0._dp,present(embedded_wall_temperature))
    controls(9)=wall_t
    if(present(wall_temperature)) then
      valid=valid.and.present(boundary_faces).and.all(ieee_is_finite(wall_temperature)).and.all(wall_temperature>0)
      controls(10:15)=wall_temperature
    endif
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) valid=valid.and.flags(2).and.present(wall_temperature)
    endif
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low_controls,21,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high_controls,21,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low_controls/=high_controls)) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
    if(.not.ready) return
    call agree(depth<tree%level_count(),ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if(.not.ready) return
    allocate(rates%fields(size(owners)),pending_halos%fields(size(owners)),pending_old%fields(size(owners)), &
      metrics%fields(size(owners)),body%fields(size(owners)),pending_flux(size(owners)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    do p=1,size(owners)
      n=tree%nodes(p)%shape
      call extract_eb_geometry_3d(extended_geometry(p),[7,7,7],n+6,center,valid)
      if(valid) valid=matching_eb_metrics_3d(center,geometry(p))
      if(.not.valid) exit
      if(owners(p)/=rank) cycle
      allocate(metrics%fields(p)%values(1,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
        rates%fields(p)%values(nv,n(1),n(2),n(3)),body%fields(p)%values(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      metrics%fields(p)%values(1,:,:,:)=extended_geometry(p)%volume_fraction
      rates%fields(p)%values=0;body%fields(p)%values=0
    enddo
    call agree(valid,ready)
    if(.not.ready) return
    if(present(internal_energy)) allocate(held_ghost)
    if(linear_halos) then
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
      inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
      halo_width=6,target_volume_fraction=metrics,restrict_snapshot=.false.,ancestor_start=ancestor_start, &
      ancestor_end=ancestor_end,temperature_start=temperature_start,temperature_end=temperature_end, &
      ancestor_fraction=ancestor_fraction,frame=frame, &
      internal_energy=internal_energy,ghost_internal_energy=held_ghost, &
      ancestor_internal_energy_start=ancestor_internal_energy_start, &
      ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_geometry=extended_geometry, &
          separate_momentum=separate_momentum, &
        consistent_velocity=consistent_velocity,paired_pressure=paired_pressure,pressure_level=target_depth, &
      connected_directions=connected_directions,preserve_temperature_seed=use_total_energy)
    else
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
      inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
      halo_width=6,target_volume_fraction=metrics,restrict_snapshot=.false.,ancestor_start=ancestor_start, &
      ancestor_end=ancestor_end,temperature_start=temperature_start,temperature_end=temperature_end, &
      ancestor_fraction=ancestor_fraction,frame=frame, &
      internal_energy=internal_energy,ghost_internal_energy=held_ghost, &
      ancestor_internal_energy_start=ancestor_internal_energy_start, &
      ancestor_internal_energy_end=ancestor_internal_energy_end)
    endif
    if(.not.ready) return
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      n=tree%nodes(p)%shape;stamp=tree_work_start()
      allocate(pending_old%fields(p)%values(nv,n(1)+12,n(2)+12,n(3)+12),stat=ios)
      valid=ios==0
      if(.not.valid) exit
      pending_old%fields(p)%values=ghost%fields(p)%values
      if(present(internal_energy)) ghost_held=held_ghost%fields(p)%values(1,:,:,:)
      ordinary=ghost%fields(p)%values
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,geometry(p),[-5,-5,-5],vf%fields(p)%values(1,:,:,:)>0,ordinary,valid)
      endif
      if(.not.valid) exit
      level=tree%nodes(p)%level+1
      faces=[origin(1,p)>0,origin(1,p)+n(1)<domain(1,level),origin(2,p)>0,origin(2,p)+n(2)<domain(2,level), &
        origin(3,p)>0,origin(3,p)+n(3)<domain(3,level)]
      if(present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      labels='outflow'
      if(present(boundary_faces)) then
        where(.not.faces) labels=boundary_faces
        faces=.true.
        do a=1,3
          if((labels(2*a-1)=='periodic').neqv.(labels(2*a)=='periodic')) labels(2*a-1:2*a)='outflow'
        enddo
      endif
      if(use_total_energy) then
      call compute_sdc_eb_transport_halo_3d(species,transport,ordinary,gt%fields(p)%values(1,:,:,:), &
        extended_geometry(p),dt,flags,rates%fields(p)%values,pending_flux(p),patch_halo,patch_theta,valid, &
        faces,labels,wall_temperature, &
        constant_coefficients=constant_coefficients)
      else
      call compute_sdc_eb_transport_halo_3d(species,transport,ordinary,gt%fields(p)%values(1,:,:,:), &
        extended_geometry(p),dt,flags,rates%fields(p)%values,pending_flux(p),patch_halo,patch_theta,valid, &
        faces,labels,wall_temperature, &
        constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
      endif
      if(.not.valid) exit
      theta=min(theta,patch_theta)
      if(wall_active) then
        allocate(body_halo(nv,n(1)+6,n(2)+6,n(3)+6),valid_body(nv,n(1),n(2),n(3)),stat=ios)
        valid=ios==0
        if(.not.valid) exit
        call compute_eb_wall_transport_source_3d(species,transport,ordinary,gt%fields(p)%values(1,:,:,:), &
          extended_geometry(p),[4,4,4],n+9,present(embedded_wall_temperature),noslip,wall_t, &
          body_halo,ignored_exchange,ignored_frequency,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
        if(.not.valid) exit
        call compute_eb_wall_transport_source_3d(species,transport,ordinary,gt%fields(p)%values(1,:,:,:), &
          extended_geometry(p),[7,7,7],n+6,present(embedded_wall_temperature),noslip,wall_t, &
          valid_body,ignored_exchange,patch_frequency,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
        if(.not.valid) exit
        patch_halo(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)=patch_halo(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)+body_halo
        rates%fields(p)%values=rates%fields(p)%values+valid_body
        body%fields(p)%values=valid_body;frequency=max(frequency,patch_frequency)
        deallocate(body_halo,valid_body)
      endif
      valid=all(ieee_is_finite(patch_halo)).and.all(ieee_is_finite(rates%fields(p)%values))
      call tree_work_stop(p,stamp)
      if(.not.valid) exit
      call move_alloc(patch_halo,pending_halos%fields(p)%values)
    enddo
    call agree(valid,ready)
    if(.not.ready) return
    call refresh_peer_source_halos(ready)
    if(.not.ready) return
    call MPI_Allreduce(theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS.and.ieee_is_finite(global_theta),ready)
    if(.not.ready) return
    if(wall_active) then
      call integrate_mpi_eb_tree_3d(tree,owners,geometry,body,exchange,comm,ready)
      if(.not.ready) return
      exchange=dt*exchange
      call MPI_Allreduce(frequency,global_frequency,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(exchange)).and.ieee_is_finite(global_frequency),ready)
      if(.not.ready) return
    endif
    if(present(registers)) then
      first_depth=max(0,depth);last_depth=first_depth
      if(depth<0) last_depth=tree%level_count()-1
      ! Non-MOL Diffusion.cpp ML StateRedist uses factor1. Stage flux
      ! acceptance (.5 for old/final) remains the outer SDC driver's job.
      do stage_depth=first_depth,last_depth
        call compute_state_redist_mpi_level_stage_3d(tree,owners,extended_geometry,pending_old,pending_halos, &
          nv,stage_depth,dt,rates,pending_registers,comm,ready,boundary_faces,target_volume_fraction,max_order, &
          register_scale=1._dp)
        if(.not.ready) return
      enddo
    endif
    if(present(registers)) call move_alloc(pending_registers,registers)
    call move_alloc(rates%fields,rate%fields);call move_alloc(pending_flux,flux)
    call move_alloc(pending_halos%fields,source_halos%fields);call move_alloc(pending_old%fields,old_halos%fields)
    minimum_theta=global_theta
    if(wall_active) then
      embedded_exchange=exchange;embedded_frequency=global_frequency
    endif
    ok=.true.
  contains
    subroutine refresh_peer_source_halos(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:),requests(:,:),destinations(:,:)
      real(dp), allocatable :: values(:,:)
      integer(int64) :: cell(3),total
      integer :: p,q,level,i,j,k,index(3),a,m(3),stage,r,allocation_status
      character(len=32) :: faces(6)
      logical :: inside,good
      faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
      call agree(good,success)
      if(.not.success) return
      total=0
      do stage=1,2
        r=0;good=.true.
        do p=1,size(owners)
          if(.not.allocated(pending_halos%fields(p)%values)) cycle
          level=tree%nodes(p)%level;m=tree%nodes(p)%shape
          do k=4,m(3)+9
          do j=4,m(2)+9
          do i=4,m(1)+9
            cell=origin(:,p)+int([i,j,k]-7,int64);inside=.true.
            do a=1,3
              if(faces(2*a-1)=='periodic') then
                cell(a)=modulo(cell(a),domain(a,level+1))
              else if(cell(a)<0.or.cell(a)>=domain(a,level+1)) then
                inside=.false.
              endif
            enddo
            if(.not.inside) cycle
            do q=1,size(owners)
              if(tree%nodes(q)%level/=level) cycle
              if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%nodes(q)%shape,int64))) cycle
              index=int(cell-origin(:,q))+1
              if(owners(q)==rank) then
                if(stage==2) pending_halos%fields(p)%values(:,i,j,k)=rates%fields(q)%values(:,index(1),index(2),index(3))
              else if(stage==1) then
                total=total+1_int64
              else
                r=r+1;requests(:,r)=[q,index];destinations(:,r)=[p,i,j,k]
              endif
              exit
            enddo
          enddo
          enddo
          enddo
        enddo
        if(stage==1) then
          call agree(total<=int(huge(1)/max(4,nv),int64),success)
          if(.not.success) return
          allocate(requests(4,int(total)),destinations(4,int(total)),stat=allocation_status)
          call agree(allocation_status==0,success)
          if(.not.success) return
        endif
      enddo
      call fetch_mpi_amr_tree_cells_3d(tree,owners,rates,nv,requests,values,comm,success)
      if(.not.success) return
      do r=1,size(requests,2)
        p=destinations(1,r);index=destinations(2:4,r)
        pending_halos%fields(p)%values(:,index(1),index(2),index(3))=values(:,r)
      enddo
    end subroutine refresh_peer_source_halos
    subroutine ordinary_transport_energy(rotation,grid,first,active,field,good)
      type(rotating_frame_3d), optional, intent(in) :: rotation
      type(eb_geometry_3d), intent(in) :: grid
      integer, intent(in) :: first(3)
      logical, intent(in) :: active(:,:,:)
      real(dp), intent(inout) :: field(:,:,:,:)
      logical, intent(out) :: good
      real(dp) :: position(3),energy
      integer :: x,y,z
      good=.true.
      if(.not.present(rotation)) return
      good=.false.
      if(.not.rotation%is_valid().or.any(shape(active)/=shape(field(1,:,:,:)))) return
      do z=1,size(active,3)
      do y=1,size(active,2)
      do x=1,size(active,1)
        if(.not.active(x,y,z)) cycle
        position=[grid%x_lower,grid%y_lower,grid%z_lower]+ &
          (real(first+[x,y,z]-1,dp)-.5_dp)*[grid%dx,grid%dy,grid%dz]
        energy=field(5,x,y,z)
        call convert_rotating_energy_3d(rotation,position,field(1,x,y,z),field(5,x,y,z),.false.,energy,good)
        if(.not.good) return
        field(5,x,y,z)=energy
      end do
      end do
      end do
      good=.true.
    end subroutine ordinary_transport_energy
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
  end subroutine compute_sdc_mpi_eb_tree_transport_halos_3d
end module mpi_eb_tree_transport_3d_mod
