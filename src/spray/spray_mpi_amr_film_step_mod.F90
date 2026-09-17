module spray_mpi_amr_film_step_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod, only: spray_liquid,spray_options,valid_spray_options,valid_spray_liquid
  use spray_mixture_mod, only: mixture_trace_fraction
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options, &
    valid_wall_component_inventory,wall_liquid_conductivity_at_temperature,mixture_wall_film_rates, &
    advance_mixture_wall_film_exchange
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_gas
  use eb_geometry_3d_mod, only: eb_geometry_3d
  implicit none
  private
  public :: advance_mpi_amr_film_coupling
contains
  ! Common film substeps, owner-local gas feedback and whole-interval publication.
  ! All model tables/options are identical across ranks by caller contract.
  subroutine advance_mpi_amr_film_coupling(comm,species,transport,components,options,wall_options,state,temperature, &
      tree,owners,films,film_mass,lower,spacing,interval,faces,substeps,boundary_change,ok,eb_geometry,nodal_phi)
    type(MPI_Comm), intent(in) :: comm
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(spray_wall_options), intent(in) :: wall_options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_patch_tree_3d), intent(in) :: tree
    type(spray_wall_film), allocatable, intent(inout) :: films(:)
    real(dp), allocatable, intent(inout) :: film_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    integer, intent(in) :: owners(:)
    character(len=*), intent(in) :: faces(6)
    integer, intent(out) :: substeps
    real(dp), intent(inout) :: boundary_change(:)
    logical, intent(out) :: ok
    type(spray_wall_film), allocatable :: patches(:),trial(:),kept(:)
    type(amr_tree_state_3d) :: gas,temp,next_gas,next_temp,field
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    logical :: periodic(3)
    character(len=32) :: reference_faces(6)
    real(dp) :: upper(3),minimum_h(3)
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:),change(:,:)
    real(dp), allocatable :: rates(:,:),hv(:,:),qg(:),qw(:)
    real(dp) :: primitive(reactive_nvar(size(species))),tg,sound,conductivity,tdot,elapsed,dt,local_dt,wall_energy
    real(dp) :: ledger,trial_ledger,global_ledger,controls(3),low(3),high(3)
    integer :: np,nc,nvar,p,c,axis,face,attempt,count_steps,ierr,node,shape3(3)
    logical :: valid,agreed,accepted,eb
    integer :: flag,lowflag,highflag
    real(dp) :: probe_offset,wall_temperature
    ok=.false.;substeps=0
    flag=merge(1,0,present(eb_geometry))+2*merge(1,0,present(nodal_phi))
    call MPI_Allreduce(flag,lowflag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,highflag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.lowflag/=highflag) return
    if (flag/=0.and.flag/=3) return
    eb=flag==3
    valid=allocated(films).and.allocated(film_mass).and.valid_spray_options(options)
    reference_faces=faces
    if (eb) reference_faces='slip_wall'
    valid=valid.and.valid_wall_options(wall_options,reference_faces).and.ieee_is_finite(interval)
    if (eb) valid=valid.and.all(wall_options%temperature==wall_options%temperature(1)).and. &
      all(faces/='periodic')
    valid=valid.and.interval>=0.and.interval<=1e6_dp.and.size(transport)==size(species)
    call agree()
    if (.not.agreed) return
    np=size(films);nc=size(components);nvar=reactive_nvar(size(species))
    valid=nc>=1.and.nc<=min(16,size(species)).and.np<=min(125000,wall_options%maximum_films)
    valid=valid.and.all(shape(film_mass)==[nc,np]).and.size(boundary_change)==nvar
    call agree()
    if (.not.agreed) return
    controls=[interval,real(nc,dp),merge(1.0_dp,0.0_dp,wall_options%evolve)]
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.any(low/=high)) return
    do p=1,np
      valid=valid.and.valid_wall_component_inventory(films(p),film_mass(:,p))
    end do
    do c=1,nc
      valid=valid.and.valid_spray_liquid(components(c)).and.components(c)%vapor_index<=size(species)
    end do
    call agree()
    if (.not.agreed) return
    reference_faces=faces
    call MPI_Bcast(reference_faces,32*6,MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(faces==reference_faces).and.all(len_trim(faces)<=32)
    call agree();if (.not.agreed) return
    periodic=faces([1,3,5])=='periodic'
    gas=state;temp=temperature;patches=films;masses=film_mass
    ! Validate the tree and closed film positions before any stencil lookup.
    call migrate_mpi_amr_wall_films(comm,tree,owners,patches,masses,lower,spacing,wall_options%maximum_films,valid,nodal_phi)
    if (.not.valid) return
    np=size(patches)
    valid=np<=125000
    call agree();if (.not.agreed) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree();if (.not.agreed) return
    upper=lower+spacing*real(tree%nodes(1)%shape,dp)
    minimum_h=spacing*real(domain(:,1),dp)/real(domain(:,size(ratios)),dp)

    allocate(positions(3,np),change(nvar,np),rates(nc,np),hv(nc,np),qg(np),qw(np))
    if (eb) then
      ! Migration already checks surface position/normal and rank agreement.
      probe_offset=max(1e-8_dp*minval((upper-lower)/real(shape(nodal_phi)-1,dp)), &
        128*epsilon(1._dp)*maxval(abs([lower,upper])))
      do p=1,np
        positions(:,p)=patches(p)%position+probe_offset*patches(p)%normal
      end do
    else
    do p=1,np
      axis=maxloc(abs(patches(p)%normal),dim=1)
      face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
      valid=wall_options%enabled(face).and.abs(patches(p)%normal(axis))==1
      if (mod(face,2)==1) then
        valid=valid.and.patches(p)%position(axis)==lower(axis)
      else
        valid=valid.and.patches(p)%position(axis)==upper(axis)
      end if
      if (.not.valid) exit
      ! Probe adjacent fluid cells without moving the exact stationary wall patch.
      positions(:,p)=max(lower+.5_dp*minimum_h,min(upper-.5_dp*minimum_h,patches(p)%position))
    end do
    end if
    call agree();if (.not.agreed) return
    change=0
    call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,gas,temp,positions,change, &
      lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi)
    if (.not.valid) return
    elapsed=0;ledger=0;count_steps=0
    do while(wall_options%evolve.and.elapsed<interval)
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      if (allocated(field%fields)) deallocate(field%fields)
      allocate(field%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        if (.not.allocated(gas%fields(node)%values)) cycle
        shape3=tree%nodes(node)%shape
        allocate(field%fields(node)%values(nvar+1,shape3(1),shape3(2),shape3(3)))
        field%fields(node)%values(:nvar,:,:,:)=gas%fields(node)%values
        field%fields(node)%values(nvar+1,:,:,:)=temp%fields(node)%values(1,:,:,:)
      end do
      call sample_mpi_amr_tree_spray_field(comm,tree,owners,field,nvar+1,positions, &
        lower,spacing,periodic,samples,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi)
      if (.not.valid) return
      local_dt=interval-elapsed;rates=0;hv=0;qg=0;qw=0
      do p=1,np
        if (patches(p)%mass<=0) cycle
        if (eb) then
          wall_temperature=wall_options%temperature(1)
        else
          axis=maxloc(abs(patches(p)%normal),dim=1)
          face=2*axis-merge(1,0,patches(p)%normal(axis)>0)
          valid=wall_options%enabled(face)
          if (.not.valid) exit
          wall_temperature=wall_options%temperature(face)
        end if
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive,tg,sound,valid)
        if (.not.valid) exit
        conductivity=wall_liquid_conductivity_at_temperature(wall_options,patches(p)%material_id,patches(p)%temperature)
        call mixture_wall_film_rates(patches(p),components,masses(:,p),species,transport,primitive,tg,options, &
          conductivity,wall_temperature,rates(:,p),hv(:,p),qg(p),qw(p),tdot,valid)
        if (.not.valid) exit
        do c=1,nc
          if (rates(c,p)>0.and.masses(c,p)>mixture_trace_fraction*patches(p)%mass) &
            local_dt=min(local_dt,options%max_fractional_change*masses(c,p)/rates(c,p))
        end do
        if (abs(tdot)>0) local_dt=min(local_dt,min(2.0_dp,options%max_fractional_change*patches(p)%temperature)/abs(tdot))
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=patches;next_mass=masses;change=0;trial_ledger=0;valid=.true.
        do p=1,np
          if (patches(p)%mass<=0) cycle
          call advance_mixture_wall_film_exchange(trial(p),components,next_mass(:,p),rates(:,p)*dt,hv(:,p), &
            qg(p)*dt,qw(p)*dt,wall_options%minimum_height,change(:,p),wall_energy,valid,deplete_traces=.true.)
          if (valid.and.trial(p)%mass>0) valid= &
            wall_liquid_conductivity_at_temperature(wall_options,trial(p)%material_id,trial(p)%temperature)>0
          if (.not.valid) exit
          trial_ledger=trial_ledger-wall_energy
        end do
        call agree()
        if (agreed) then
          next_gas=gas;next_temp=temp
          call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,next_gas,next_temp,positions,change, &
            lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi)
          if (valid) then
            accepted=.true.;exit
          end if
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      gas=next_gas;temp=next_temp;patches=trial;masses=next_mass;ledger=ledger+trial_ledger
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    kept=pack(patches,patches%mass>0)
    allocate(kept_mass(nc,size(kept)))
    do c=1,nc
      kept_mass(c,:)=pack(masses(c,:),patches%mass>0)
    end do
    call migrate_mpi_amr_wall_films(comm,tree,owners,kept,kept_mass,lower,spacing,wall_options%maximum_films,valid,nodal_phi)
    if (.not.valid) return
    call MPI_Allreduce(ledger,global_ledger,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (.not.ieee_is_finite(global_ledger)) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temp%fields,temperature%fields)
    call move_alloc(kept,films);call move_alloc(kept_mass,film_mass)
    boundary_change=0;boundary_change(5)=global_ledger;substeps=count_steps;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
