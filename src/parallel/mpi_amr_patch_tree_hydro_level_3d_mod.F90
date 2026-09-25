module mpi_amr_patch_tree_hydro_level_3d_mod
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: compute_reactive_face_fluxes_3d,compute_reactive_plm_face_fluxes_3d, &
    compute_reactive_ppm_halo_ctu_fluxes_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_ghost_3d_mod, only: apply_amr_tree_pressure_ghosts_3d,amr_tree_inflow_faces_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_hydro_level_3d
contains
  subroutine advance_mpi_amr_tree_hydro_level_3d(tree,owners,species,state,temperature,depth,root_spacing,dt,faces, &
      reconstruction,limiter,solver,integrated_flux,comm,ok,coarse_start,coarse_end,fractions_start,fractions_end, &
      inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure,difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    integer, intent(in) :: depth
    type(amr_tree_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    type(amr_tree_flux_3d), allocatable :: first_flux(:),second_flux(:)
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: pressure(6),mach(6),pin(6),tin(6),viscosity
    type(amr_tree_state_3d) :: initial,candidate,temps
    real(dp), allocatable :: h(:,:),inlet(:,:)
    integer :: flags(5),minimum_flags(5),maximum_flags(5)
    integer :: node,nvar,rank,status,ios,parent,width
    real(dp) :: controls(29),low(29),high(29)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(9),reference_options(9)
    character(len=:), allocatable :: donor_limiter
    logical :: valid
    ok=.false.
    nvar=reactive_nvar(size(species))
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    viscosity=0
    if(reconstruction=='characteristic_ppm'.and.solver=='pelec') viscosity=.1_dp
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
    valid=valid.and.viscosity>=0.and.(viscosity==0.or.reconstruction=='characteristic_ppm')
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    valid=valid.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.(reconstruction=='pcm'.or.reconstruction=='characteristic_plm'.or. &
      reconstruction=='characteristic_ppm')
    valid=valid.and.(limiter=='mc'.or.limiter=='minmod').and.len_trim(solver)<=64
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
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
    options(:6)=faces;options(7)=reconstruction;options(8)=limiter;options(9)=solver
    reference_options=options
    call MPI_Bcast(reference_options,64*9,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(options==reference_options)
    call agree(valid)
    if (.not.valid) return
    flags=[depth,merge(1,0,present(coarse_start)),merge(1,0,present(coarse_end)), &
      merge(1,0,present(fractions_start)),merge(1,0,present(fractions_end))]
    call MPI_Allreduce(flags,minimum_flags,5,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,maximum_flags,5,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(minimum_flags/=maximum_flags)) return
    if (depth<0.or.depth>=tree%level_count().or.any(flags(2:)/=flags(2))) return
    if (depth>0.and.flags(2)==0) return
    if (present(coarse_start)) then
      valid=size(fractions_start)==depth.and.size(fractions_end)==depth
      if (valid) valid=all(ieee_is_finite(fractions_start)).and.all(ieee_is_finite(fractions_end))
      if (valid) valid=all(fractions_start>=0).and.all(fractions_end<=1).and.all(fractions_end>=fractions_start)
      call agree(valid)
      if (.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_start,nvar+1,comm,valid,maximum_depth=depth-1)
      if (.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_end,nvar+1,comm,valid,maximum_depth=depth-1)
      if (.not.valid) return
      do node=1,size(tree%nodes)
        if (owners(node)/=rank.or.tree%nodes(node)%level>=depth) cycle
        if (.not.all(ieee_is_finite(coarse_start%fields(node)%values))) valid=.false.
        if (.not.all(ieee_is_finite(coarse_end%fields(node)%values))) valid=.false.
      end do
      call agree(valid)
      if (.not.valid) return
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
    call copy_local(initial,candidate,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,first_flux,1,valid)
    if (.not.valid) return
    if (reconstruction/='characteristic_ppm') then
      call euler_stage(candidate,temps,second_flux,2,valid)
      if (.not.valid) return
    end if
    do node=1,size(tree%nodes)
      if (owners(node)/=rank.or.tree%nodes(node)%level/=depth) cycle
      if (reconstruction/='characteristic_ppm') then
        first_flux(node)%x=first_flux(node)%x+second_flux(node)%x
        first_flux(node)%y=first_flux(node)%y+second_flux(node)%y
        first_flux(node)%z=first_flux(node)%z+second_flux(node)%z
      end if
      if (.not.all(ieee_is_finite(first_flux(node)%x)).or..not.all(ieee_is_finite(first_flux(node)%y)).or. &
          .not.all(ieee_is_finite(first_flux(node)%z))) valid=.false.
    end do
    call agree(valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    call move_alloc(first_flux,integrated_flux)
    ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success,depth=depth)
    end subroutine
    subroutine euler_stage(u,t,flux,stage,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      integer, intent(in) :: stage
      real(dp) :: work_stamp
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated,ghosts
      type(amr_tree_flux_3d), allocatable, intent(out) :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),fractions(:)
      real(dp) :: rhs(nvar),theta,flux_weight
      integer :: peer,n(3),dims(3),i,j,k,allocation_status
      success=.false.;flux_weight=.5_dp*dt
      if (reconstruction=='characteristic_ppm') flux_weight=dt
      if (present(fractions_start)) then
        if (stage==1) then
          fractions=fractions_start
        else
          fractions=fractions_end
        end if
      end if
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
        if (tree%nodes(peer)%level<depth) packed%fields(peer)%values=coarse_start%fields(peer)%values
      end do
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      call copy_local(u,updated,success)
      if (.not.success) return
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,nvar+1,faces,ghosts,comm,success, &
        inflow_state=inlet,coarse_limiter=donor_limiter,species_count=size(species), &
        coarse_future=coarse_end,coarse_fractions=fractions,depth=depth,ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank.or.tree%nodes(peer)%level/=depth) cycle
        n=tree%nodes(peer)%shape;dims=n+4
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
              logical :: inflow_mask(6),wall_mask(6)
              call amr_tree_inflow_faces_3d(tree,peer,faces,inflow_mask,success,wall_mask)
              if (.not.success) exit
              fx=0;fy=0;fz=0
            call compute_reactive_ppm_halo_ctu_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              n,h(:,peer),dt,.false.,.true.,solver,fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2), &
              fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),theta,success, &
              inflow_faces=inflow_mask,no_slip_faces=wall_mask,difmag=viscosity)
            end block
          else if (reconstruction=='pcm') then
            call compute_reactive_face_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              dims(1),dims(2),dims(3),solver,fx,fy,fz,success)
          else
            call compute_reactive_plm_face_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              dims(1),dims(2),dims(3),limiter,solver,fx,fy,fz,success)
          end if
          if (.not.success) exit
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                rhs=-(fx(:,i+2,j+2,k+2)-fx(:,i+1,j+2,k+2))/h(1,peer) &
                  -(fy(:,i+2,j+2,k+2)-fy(:,i+2,j+1,k+2))/h(2,peer) &
                  -(fz(:,i+2,j+2,k+2)-fz(:,i+2,j+2,k+1))/h(3,peer)
                updated%fields(peer)%values(:,i,j,k)=u%fields(peer)%values(:,i,j,k)+dt*rhs
              end do
            end do
          end do
          allocate(flux(peer)%x,source=flux_weight*fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%y,source=flux_weight*fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%z,source=flux_weight*fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          deallocate(fx,fy,fz)
        end associate
        call tree_work_stop(peer,work_stamp)
      end do
      call agree(success)
      if (.not.success) return
      if (stage==2) then
        do peer=1,size(tree%nodes)
          if (owners(peer)/=rank.or.tree%nodes(peer)%level/=depth) cycle
          updated%fields(peer)%values=.5_dp*initial%fields(peer)%values+.5_dp*updated%fields(peer)%values
        end do
      end if
      call recover_tree(updated,t,success)
      if (.not.success) return
      call move_alloc(updated%fields,u%fields)
    end subroutine
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
  end subroutine advance_mpi_amr_tree_hydro_level_3d
end module
