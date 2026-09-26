module mpi_amr_patch_tree_hydro_3d_mod
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: compute_reactive_face_fluxes_3d,compute_reactive_plm_face_fluxes_3d, &
    compute_reactive_ppm_halo_ctu_fluxes_3d,adjust_reactive_artificial_viscosity_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
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
  public :: advance_mpi_amr_tree_hydro_3d,compute_sdc_mpi_amr_tree_hydro_3d
contains
  subroutine compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,state,temperature,predictor_source,root_spacing, &
      dt,faces,hydro_rate,integrated_flux,comm,ok,inflow_state,inflow_temperature,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag, &
      target_depth,coarse_start,coarse_end,coarse_fractions,source_start,source_end,source_fractions)
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
    type(amr_tree_state_3d) :: local_state,local_temperature
    real(dp) :: budget(5+size(species)),maximum
    local_state=state;local_temperature=temperature
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local_state,local_temperature,root_spacing,dt,faces, &
      'characteristic_ppm','mc','pelec',budget,maximum,comm,ok,inflow_state,inflow_temperature, &
      outlet_pressure,inflow_mach,inflow_pressure,difmag,predictor_source,hydro_rate,integrated_flux, &
      target_depth,coarse_start,coarse_end,coarse_fractions,source_start,source_end,source_fractions)
  end subroutine compute_sdc_mpi_amr_tree_hydro_3d

  subroutine advance_mpi_amr_tree_hydro_3d(tree,owners,species,state,temperature,root_spacing,dt,faces, &
      reconstruction,limiter,solver,boundary_change,maximum_correction,comm,ok,inflow_state,inflow_temperature, &
      outlet_pressure,inflow_mach,inflow_pressure,difmag,predictor_source,hydro_rate,integrated_flux, &
      target_depth,coarse_start,coarse_end,coarse_fractions,source_start,source_end,source_fractions)
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
    real(dp) :: pressure(6),mach(6),pin(6),tin(6),viscosity
    type(amr_tree_state_3d) :: initial,candidate,temps,stage_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)),first_max,second_max
    integer :: node,nvar,rank,status,ios,parent,width,depth,selection(7),selection_low(7),selection_high(7)
    real(dp) :: controls(29),low(29),high(29)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(10),reference_options(10)
    character(len=:), allocatable :: donor_limiter
    logical :: valid,stage_only
    ok=.false.
    stage_only=present(predictor_source)
    nvar=reactive_nvar(size(species))
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
      merge(1,0,present(source_end)),merge(1,0,present(source_fractions))]
    call MPI_Allreduce(selection,selection_low,7,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selection,selection_high,7,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(selection_low/=selection_high)) return
    valid=depth>=-1.and.depth<tree%level_count()
    if(present(target_depth)) valid=valid.and.depth>=0.and.stage_only
    valid=valid.and.all(selection(2:4)==selection(2)).and.all(selection(5:7)==selection(5))
    if(any(selection(2:)>0)) valid=valid.and.depth>=0.and.stage_only
    call agree(valid)
    if(.not.valid) return
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
    controls=[root_spacing,dt,pressure,mach,tin,pin,viscosity]
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
    call MPI_Allreduce(controls,low,29,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,29,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
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
      allocate(inlet(nvar+1,6),reference_inlet(nvar+1,6),stat=ios)
      valid=ios==0;call agree(valid)
      if (.not.valid) return
      inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
      reference_inlet=inlet
      call MPI_Bcast(reference_inlet,size(reference_inlet),MPI_DOUBLE_PRECISION,0,comm,status)
      if (status/=MPI_SUCCESS) return
      valid=all(inlet==reference_inlet);call agree(valid)
      if (.not.valid) return
    end if
    if (reconstruction/='pcm') donor_limiter=limiter
    width=2
    if (reconstruction=='characteristic_ppm') width=4
    allocate(h(3,size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    h(:,1)=root_spacing
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if (any(h(:,node)<=0)) return
    end do
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
      type(amr_tree_state_3d) :: packed,updated,ghosts,source_ghosts,predictor
      type(amr_tree_flux_3d), allocatable :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),source_values(:,:,:,:)
      character(len=32) :: source_faces(6)
      real(dp) :: correction(nvar),rhs(nvar),theta
      integer :: peer,n(3),dims(3),i,j,k,allocation_status
      success=.false.;budget=0;largest=0;correction=0
      allocate(packed%fields(size(tree%nodes)),flux(size(tree%nodes)),stat=allocation_status)
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(nvar+1,n(1),n(2),n(3)),stat=allocation_status)
        if (allocation_status/=0) exit
        packed%fields(peer)%values(:nvar,:,:,:)=u%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=t%fields(peer)%values(1,:,:,:)
      end do
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      call copy_local(u,updated,success)
      if (.not.success) return
      if(stage_only) then
        do peer=1,size(tree%nodes)
          if(owners(peer)==rank) updated%fields(peer)%values=0
        end do
        predictor=predictor_source
        if(present(coarse_start)) then
          call set_ancestors(packed,coarse_start,nvar+1,success)
          if(.not.success) return
        end if
        if(present(source_start)) then
          call set_ancestors(predictor,source_start,nvar,success)
          if(.not.success) return
        end if
      end if
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,nvar+1,faces,ghosts,comm,success, &
        inflow_state=inlet,coarse_limiter=donor_limiter,species_count=size(species),ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        depth=target_depth,coarse_future=coarse_end,coarse_fractions=coarse_fractions)
      if (.not.success) return
      if(stage_only) then
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
        if(stage_only) source_values=source_ghosts%fields(peer)%values
        work_stamp=tree_work_start()
        associate(ghost=>ghosts%fields(peer)%values)
          call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,pressure,ghost,success,ghost_width=width)
          if (.not.success) exit
          allocate(fx(nvar,dims(1),dims(2),dims(3)),fy(nvar,dims(1),dims(2),dims(3)), &
            fz(nvar,dims(1),dims(2),dims(3)),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          if (reconstruction=='characteristic_ppm') then
            block
              logical :: inflow_mask(6),wall_mask(6),reflecting_mask(6),outflow_mask(6)
              call amr_tree_inflow_faces_3d(tree,peer,faces,inflow_mask,success,wall_mask,reflecting_mask,outflow_mask)
              if (.not.success) exit
              fx=0;fy=0;fz=0
            call compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              n,h(:,peer),dt,.false.,.true.,solver,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
              fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),theta,success, &
              inflow_faces=inflow_mask,no_slip_faces=wall_mask,difmag=viscosity, &
              reflecting_faces=reflecting_mask,outflow_faces=outflow_mask,outlet_pressure=outlet_pressure, &
              conserved_source=source_values)
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
          allocate(flux(peer)%x,source=dt*fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%y,source=dt*fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%z,source=dt*fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),stat=allocation_status)
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
          if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
          if(.not.all(ieee_is_finite(flux(peer)%x)).or..not.all(ieee_is_finite(flux(peer)%y)).or. &
            .not.all(ieee_is_finite(flux(peer)%z))) success=.false.
        end do
        call agree(success)
        if(.not.success) return
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
