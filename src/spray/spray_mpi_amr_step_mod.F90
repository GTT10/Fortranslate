module spray_mpi_amr_step_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties, &
    advance_mixture_frozen_gas,mixture_transfer_rates,mixture_trace_fraction,build_mixture_film
  use spray_breakup_mod, only: advance_spray_breakup
  use spray_coupling_3d_mod, only: valid_spray_boundary_faces,apply_parcel_boundaries
  use spray_mpi_mod, only: advance_mpi_khrt_shedding,impact_mpi_planar_walls
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,valid_wall_options
  use spray_mpi_amr_film_mod, only: migrate_mpi_amr_wall_films
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_gas
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels,migrate_eb_tree_spray_parcels
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use spray_eb_drift_mod, only: drift_spray_level_set
  use spray_eb_intersection_mod, only: intersect_spray_level_set
  use spray_eb_impact_mod, only: spray_eb_wall_options
  use spray_mpi_eb_impact_mod, only: advance_mpi_level_set_impacts
  implicit none
  private
  public :: advance_mpi_amr_mixture_coupling
contains
  ! Regular or nodal EB AMR: common particle substeps and whole-interval rollback.
  ! Species, transport, liquids and options must be identical on every rank.
  ! Optional gas_change sums accepted deposition increments per fluid volume;
  ! rejected trials and wall/outflow inventories do not enter this source.
  subroutine advance_mpi_amr_mixture_coupling(comm,species,transport,components,options,state,temperature, &
      tree,owners,parcels,component_mass,lower,spacing,interval,substeps,ok,last_id,boundary_faces,boundary_change, &
      wall_films,film_mass,film_settings,eb_geometry,nodal_phi,model_rejected,eb_outflow_change,gas_change, &
      motion_mode,maximum_level)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(amr_tree_state_3d), optional, intent(inout) :: gas_change
    type(amr_tree_state_3d) :: integrated_change
    type(amr_tree_state_3d), allocatable :: trial_change
    integer :: export_flag,export_low,export_high
    integer, optional, intent(in) :: maximum_level
    character(len=*), optional, intent(in) :: motion_mode
    character(len=16) :: motion
    integer :: motion_code,motion_low,motion_high
    real(dp) :: kick_fraction,kick_dt,motion_dt
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), allocatable, intent(inout) :: component_mass(:,:)
    real(dp), intent(in) :: lower(3),spacing(3),interval
    integer, intent(out) :: substeps
    logical, intent(out) :: ok
    integer(int64), optional, intent(inout) :: last_id
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(inout) :: boundary_change(:)
    type(spray_wall_film), allocatable, optional, intent(inout) :: wall_films(:)
    real(dp), allocatable, optional, intent(inout) :: film_mass(:,:)
    type(spray_wall_options), optional, intent(in) :: film_settings
    type(eb_geometry_3d), optional, intent(in) :: eb_geometry(:)
    real(dp), optional, intent(in) :: nodal_phi(0:,0:,0:)
    real(dp), optional, intent(inout) :: eb_outflow_change(:)
    real(dp), allocatable :: impact_exterior(:),liquid_energies(:)
    character(len=16) :: outflow_faces(6)
    real(dp) :: exterior(reactive_nvar(size(species))),finish(3),delta(3),flight_fraction
    integer :: axis
    logical, optional, intent(out) :: model_rejected
    type(spray_eb_wall_options) :: eb_walls
    logical :: rejected
    integer :: eb_hits
    type(spray_wall_film), allocatable :: patches(:),trial_patches(:)
    real(dp), allocatable :: patch_mass(:,:),trial_patch_mass(:,:),impact_previous(:,:),impact_pressure(:)
    real(dp) :: mw(size(components)),impact_wall(5),wall_ledger(5)
    logical, allocatable :: hit(:),eligible(:)
    type(spray_parcel) :: preview(1)
    logical :: preview_contact(1)
    real(dp) :: preview_change(5),preview_exterior(5)
    integer :: preview_hits
    logical :: crossing
    real(dp) :: crossing_fraction,crossing_position(3),crossing_normal(3)
    integer :: crossing_cell(3)
    integer :: film_flag
    real(dp) :: ledger(reactive_nvar(size(species))),trial_ledger(reactive_nvar(size(species)))
    real(dp) :: global_ledger(reactive_nvar(size(species))),exchange(5)
    character(len=32) :: faces(6),reference_faces(6)
    integer(int64) :: accepted_id,trial_id
    real(dp), allocatable :: reynolds(:),child_energy(:)
    integer :: child_flag,low_flag,high_flag
    type(spray_parcel), allocatable :: drops(:),trial(:),survivors(:)
    type(spray_liquid) :: effective
    type(spray_film) :: film
    real(dp) :: old_mass
    type(amr_tree_state_3d) :: gas,temp,next_gas,next_temp,field
    real(dp), allocatable :: masses(:,:),next_mass(:,:),kept_mass(:,:),positions(:,:),samples(:,:)
    real(dp), allocatable :: primitive(:,:),tg(:),change(:,:)
    real(dp) :: elapsed,dt,local_dt,min_interval,max_interval,sound,invtau,heat,rho,cp,e,tdot,speed
    real(dp) :: y(size(components)),rates(size(components)),hv(size(components)),length(3)
    integer :: p,n,np,nvar,nc,ierr,attempt,count_steps,kept,min_nc,max_nc,node,shape3(3)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: minimum_spacing
    logical :: periodic(3)
    logical :: valid,agreed,accepted
    ok=.false.;substeps=0;faces='periodic';ledger=0;wall_ledger=0
    exterior=0
    if(present(model_rejected)) model_rejected=.false.
    export_flag=merge(1,0,present(gas_change))
    call MPI_Allreduce(export_flag,export_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(export_flag,export_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.export_low/=export_high) return
    motion='coupled';kick_fraction=1
    if(present(motion_mode)) motion=motion_mode
    select case(motion)
    case('coupled')
      motion_code=0
    case('kick_drift')
      motion_code=1;kick_fraction=.5_dp
    case('kick')
      motion_code=2;kick_fraction=.5_dp
    case default
      motion_code=-1
    end select
    call MPI_Allreduce(motion_code,motion_low,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(motion_code,motion_high,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.motion_low<0.or.motion_low/=motion_high) return
    if(present(eb_geometry)) faces="slip_wall"
    if (present(boundary_faces)) faces=boundary_faces
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
    valid=valid.and.(present(eb_geometry).eqv.present(nodal_phi))
    if(present(eb_geometry).and.present(boundary_faces)) then
      valid=valid.and.all(faces=='periodic'.or.faces=='outflow'.or.faces=='slip_wall'.or.faces=='no_slip_wall')
      valid=valid.and.present(eb_outflow_change)
    end if
    valid=valid.and.ieee_is_finite(interval)
    valid=valid.and.interval>=0.and.interval<=1e6_dp
    valid=valid.and.(.not.options%khrt_children.or.present(last_id))
    valid=valid.and.size(components)>=1.and.size(components)<=16.and.valid_spray_boundary_faces(faces)
    if (present(boundary_faces)) then
      do n=1,6
        valid=valid.and.len_trim(boundary_faces(n))<=len(faces(n))
      end do
    end if
    if (present(boundary_change)) valid=valid.and.size(boundary_change)==size(ledger)
    if(present(eb_outflow_change)) then
      valid=valid.and.present(eb_geometry)
      valid=valid.and.size(eb_outflow_change)==size(ledger)
    end if
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.allocated(wall_films).and.present(last_id).and. &
      (present(boundary_faces).or.present(eb_geometry))
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    if (present(film_settings)) then
      if(present(eb_geometry)) then
        outflow_faces='slip_wall'
        valid=valid.and.valid_wall_options(film_settings,outflow_faces)
      else
        valid=valid.and.valid_wall_options(film_settings,faces)
      end if
      if (present(eb_geometry)) then
        valid=valid.and.all(film_settings%temperature==film_settings%temperature(1))
        eb_walls%temperature=film_settings%temperature(1)
        eb_walls%contact_angle=film_settings%contact_angle
        eb_walls%maximum_films=film_settings%maximum_films
      end if
    end if
    call agree()
    if (.not.agreed) return
    film_flag=merge(1,0,present(wall_films))+2*merge(1,0,present(eb_geometry))+4*merge(1,0,present(nodal_phi))
    film_flag=film_flag+8*merge(1,0,present(eb_outflow_change))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    if(present(eb_outflow_change)) allocate(impact_exterior(5+size(components)))
    reference_faces=faces
    call MPI_Bcast(reference_faces,32*6,MPI_CHARACTER,0,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(faces==reference_faces)
    call agree()
    if (.not.agreed) return
    periodic=[faces(1)=='periodic',faces(3)=='periodic',faces(5)=='periodic']
    child_flag=merge(1,0,options%khrt_children)
    call MPI_Allreduce(child_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(child_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
    accepted_id=0
    if (present(last_id)) accepted_id=last_id
    np=size(parcels);nc=size(components);nvar=reactive_nvar(size(species))
    valid=size(component_mass,1)==nc.and.size(component_mass,2)==np.and.np<=125000
    valid=valid.and.size(transport)==size(species)
    call MPI_Allreduce(interval,min_interval,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(interval,max_interval,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,min_nc,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,max_nc,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=valid.and.min_interval==max_interval.and.min_nc==max_nc
    call agree()
    if (.not.agreed) return
    do p=1,np
      valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
    end do
    do n=1,nc
      valid=valid.and.valid_spray_liquid(components(n)).and.components(n)%vapor_index<=size(species)
    end do
    call agree()
    if (.not.agreed) return
    gas=state;temp=temperature;drops=parcels;masses=component_mass
    if (present(wall_films)) then
      patches=wall_films;patch_mass=film_mass
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
        film_settings%maximum_films,valid,nodal_phi,maximum_level)
      if (.not.valid) return
      do n=1,nc
        mw(n)=species(components(n)%vapor_index)%molecular_weight
      end do
    end if
    allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np))
    do p=1,np
      positions(:,p)=drops(p)%position
    end do
    change=0
    ! This collective also validates geometry and the original gas on every rank.
    call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,gas,temp,positions,change, &
      lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree()
    if (.not.agreed) return
    if(present(gas_change)) then
      integrated_change=gas
      do node=1,size(owners)
        if(.not.allocated(integrated_change%fields(node)%values)) cycle
        integrated_change%fields(node)%values=0
      end do
      allocate(trial_change)
    end if
    length=spacing*real(tree%nodes(1)%shape,dp)
    node=size(domain,2)
    if(present(maximum_level)) node=maximum_level+1
    minimum_spacing=minval(length/real(domain(:,node),dp))
    elapsed=0;count_steps=0
    do while(elapsed<interval)
      valid=count_steps<options%maximum_substeps
      call agree()
      if (.not.agreed) return
      np=size(drops)
      deallocate(positions,change,primitive,tg)
      if (allocated(reynolds)) deallocate(reynolds)
      allocate(positions(3,np),change(nvar,np),primitive(nvar,np),tg(np),reynolds(np))
      primitive=0;tg=0
      ! Sample the temperature guess alongside conserved variables in one exchange.
      if (allocated(field%fields)) deallocate(field%fields)
      allocate(field%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        if (.not.allocated(gas%fields(node)%values)) cycle
        shape3=tree%nodes(node)%shape
        allocate(field%fields(node)%values(nvar+1,shape3(1),shape3(2),shape3(3)))
        field%fields(node)%values(:nvar,:,:,:)=gas%fields(node)%values
        field%fields(node)%values(nvar+1,:,:,:)=temp%fields(node)%values(1,:,:,:)
      end do
      do p=1,np
        positions(:,p)=drops(p)%position
      end do
      call sample_mpi_amr_tree_spray_field(comm,tree,owners,field,nvar+1,positions, &
        lower,spacing,periodic,samples,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi,maximum_level=maximum_level)
      if (.not.valid) return
      local_dt=interval-elapsed
      do p=1,np
        if (drops(p)%mass<=0) cycle
        call reactive_conserved_to_primitive(species,samples(:nvar,p),samples(nvar+1,p),primitive(:,p),tg(p),sound,valid)
        if (.not.valid) exit
        y=masses(:,p)/drops(p)%mass
        call mixture_transfer_rates(species,transport,primitive(:,p),tg(p),drops(p),components,y,options, &
          invtau,heat,rates,hv,valid)
        if (.not.valid) exit
        call liquid_mixture_properties(components,y,drops(p)%temperature,rho,cp,e,valid)
        if (.not.valid) exit
        if (invtau>0) local_dt=min(local_dt,options%max_fractional_change/invtau/kick_fraction)
        tdot=heat*(tg(p)-drops(p)%temperature)
        do n=1,nc
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            local_dt=min(local_dt,options%max_fractional_change*masses(n,p)/rates(n)/kick_fraction)
          tdot=tdot-rates(n)*(hv(n)-liquid_energy(components(n),drops(p)%temperature))
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) local_dt=min(local_dt, &
          min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot)/kick_fraction)
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4,p))))
        if (speed>0.and.motion/='kick') local_dt=min(local_dt,options%particle_cfl*minimum_spacing/speed)
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        kick_dt=kick_fraction*dt;motion_dt=dt
        if(motion=='kick') motion_dt=0
        trial=drops;next_mass=masses;change=0;valid=.true.;trial_id=accepted_id;reynolds=0;trial_ledger=0
        do p=1,np
          if (drops(p)%mass<=0) cycle
          call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
            next_mass(:,p),options,kick_dt,change(:,p),valid)
          if (.not.valid) exit
          if (trial(p)%mass>0.and.options%breakup_model/='none') then
            old_mass=trial(p)%mass;y=next_mass(:,p)/old_mass
            call liquid_mixture_properties(components,y,trial(p)%temperature,rho,cp,e,valid)
            if (.not.valid) exit
            effective=components(1);effective%density=rho;effective%density_model='constant'
            effective%density_coefficients=0;effective%density_temperature_range=0
            effective%heat_capacity=cp;effective%heat_capacity_model='constant'
            effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
            effective%reference_temperature=trial(p)%temperature;effective%reference_energy=e
            effective%boiling_model='reject';effective%critical_temperature=0
            effective%normal_boiling_temperature=0
            call build_mixture_film(species,transport,primitive(:,p),tg(p),trial(p),components,y,options, &
              film,rates,hv,valid)
            if (.not.valid) exit
            if (options%khrt_children) then
              reynolds(p)=film%density*norm2(primitive(2:4,p)-trial(p)%velocity)* &
                parcel_diameter(trial(p),effective)/film%viscosity
            else
              call advance_spray_breakup(trial(p),effective,options,primitive(1,p),primitive(2:4,p), &
                film%density,film%viscosity,kick_dt,valid)
              if (.not.valid) exit
              next_mass(:,p)=next_mass(:,p)*(trial(p)%mass/old_mass)
            end if
          end if
          if(motion/='coupled') then
            trial(p)%position=drops(p)%position
            if(motion=='kick_drift'.and..not.present(eb_geometry)) &
              trial(p)%position=trial(p)%position+motion_dt*trial(p)%velocity
          end if
          if(present(eb_geometry)) then
            trial(p)%position=drops(p)%position
          else if(.not.present(boundary_faces)) then
            trial(p)%position=lower+modulo(trial(p)%position-lower,length)
          end if
        end do
        call agree()
        if (agreed.and.options%khrt_children) then
          if(allocated(eligible)) deallocate(eligible)
          allocate(eligible(np));eligible=.true.
          if(present(eb_geometry).and.motion_dt>0) then
            ! A parent crossing the EB takes the wall path, with no KH release
            ! on that substep. Unhit parents and their children then drift together.
            do p=1,np
              if(trial(p)%mass<=0) cycle
              if(present(boundary_faces)) then
                if(any(faces=='periodic')) then
                  ! Use the same segmented flight as the eventual drift, including
                  ! EB or reflecting-wall contacts after a periodic crossing.
                  preview(1)=trial(p);preview_contact=.false.
                  preview_change=0;preview_exterior=0;preview_hits=0
                  call drift_spray_level_set(nodal_phi,lower,lower+length,motion_dt,32,preview, &
                    preview_change,preview_hits,valid,faces,[0._dp],preview_exterior,preview_contact)
                  if(.not.valid) exit
                  eligible(p)=.not.preview_contact(1)
                  cycle
                end if
              end if
              crossing=.false.;crossing_fraction=0;crossing_position=0;crossing_normal=0;crossing_cell=0
              delta=motion_dt*trial(p)%velocity;finish=trial(p)%position+delta
              if(present(eb_outflow_change)) then
                flight_fraction=1
                do axis=1,3
                  if(delta(axis)>0.and.finish(axis)>=lower(axis)+length(axis)) then
                    flight_fraction=min(flight_fraction,(lower(axis)+length(axis)-trial(p)%position(axis))/delta(axis))
                  else if(delta(axis)<0.and.finish(axis)<=lower(axis)) then
                    flight_fraction=min(flight_fraction,(lower(axis)-trial(p)%position(axis))/delta(axis))
                  end if
                end do
                finish=max(lower,min(lower+length,trial(p)%position+flight_fraction*delta))
              end if
              call intersect_spray_level_set(nodal_phi,lower,lower+length,trial(p)%position, &
                finish,crossing,crossing_fraction, &
                crossing_position,crossing_normal,crossing_cell,valid)
              if(.not.valid) exit
              eligible(p)=.not.crossing
              if(present(boundary_faces)) then
                ! Domain-wall impact also takes precedence over KH release.
                do axis=1,3
                  if(delta(axis)<0.and.trial(p)%position(axis)+delta(axis)<=lower(axis)) then
                    if(faces(2*axis-1)=='slip_wall'.or.faces(2*axis-1)=='no_slip_wall') eligible(p)=.false.
                  else if(delta(axis)>0.and.trial(p)%position(axis)+delta(axis)>=lower(axis)+length(axis)) then
                    if(faces(2*axis)=='slip_wall'.or.faces(2*axis)=='no_slip_wall') eligible(p)=.false.
                  end if
                end do
              end if
            end do
            call agree()
          end if
          if(agreed) call advance_mpi_khrt_shedding(comm,components,options,trial,next_mass,primitive(1,:), &
            primitive(2:4,:),reynolds,kick_dt,trial_id,child_energy,valid,eligible)
          if (valid) valid=size(trial)<=125000
          call agree()
          if (agreed) change(5,:)=change(5,:)+child_energy
        end if
        if(allocated(impact_exterior)) impact_exterior=0
        if(agreed.and.present(eb_geometry).and..not.present(wall_films).and.motion_dt>0) then
          exchange=0;eb_hits=0
          if(allocated(impact_exterior)) then
            if(allocated(liquid_energies)) deallocate(liquid_energies)
            allocate(liquid_energies(size(trial)));liquid_energies=0;outflow_faces='outflow'
            if(present(boundary_faces)) outflow_faces=faces
            do p=1,size(trial)
              if(trial(p)%mass<=0) cycle
              call liquid_mixture_properties(components,next_mass(:,p)/trial(p)%mass,trial(p)%temperature,rho,cp,e,valid)
              if(.not.valid) exit
              liquid_energies(p)=e
            end do
            if(valid) call drift_spray_level_set(nodal_phi,lower,lower+length,motion_dt,32,trial,exchange,eb_hits,valid, &
              outflow_faces,liquid_energies,impact_exterior(:5))
            if(valid) then
              do p=1,size(trial)
                if(trial(p)%mass>0) cycle
                impact_exterior(6:)=impact_exterior(6:)-next_mass(:,p)*trial(p)%multiplicity
                next_mass(:,p)=0
              end do
            end if
            call agree()
            if(agreed) then
              call MPI_Allreduce(MPI_IN_PLACE,impact_exterior,size(impact_exterior),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
              valid=ierr==MPI_SUCCESS
            end if
          else
            call drift_spray_level_set(nodal_phi,lower,lower+length,motion_dt,32,trial,exchange,eb_hits,valid)
          end if
          if(valid) trial_ledger(:5)=exchange
          call agree()
        end if
        impact_wall=0
        if (agreed.and.present(wall_films)) then
          trial_patches=patches;trial_patch_mass=patch_mass
          call prepare_impacts()
          call agree()
          if (agreed.and.motion_dt>0) then
            if (present(eb_geometry)) then
              eb_hits=0
              call advance_mpi_level_set_impacts(comm,components,options,eb_walls,nodal_phi,mw,impact_pressure, &
                lower,lower+length,motion_dt,trial,next_mass,trial_patches,trial_patch_mass,trial_id,impact_wall, &
                eb_hits,valid,rejected,impact_exterior,boundary_faces)
              if (rejected) then
                if (present(model_rejected)) model_rejected=.true.
                return
              end if
            else
              call impact_mpi_planar_walls(comm,components,options,film_settings,faces,mw,impact_pressure,impact_previous, &
                lower,lower+length,trial,next_mass,trial_patches,trial_patch_mass,trial_id,hit,impact_wall,valid)
            end if
            if (valid) valid=size(trial)<=125000
            call agree()
          end if
        end if
        if (agreed) then
          if (present(boundary_faces).and.motion_dt>0) call apply_boundaries()
          call agree()
        end if
        if (agreed) then
          next_gas=gas;next_temp=temp
          call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,next_gas,next_temp,positions,change, &
            lower,spacing,periodic,valid,eb_geometry=eb_geometry,nodal_phi=nodal_phi, &
            gas_change=trial_change,maximum_level=maximum_level)
          if (valid) then
            accepted=.true.;exit
          end if
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      if(present(gas_change)) then
        do node=1,size(owners)
          if(.not.allocated(integrated_change%fields(node)%values)) cycle
          associate(a=>integrated_change%fields(node)%values,b=>trial_change%fields(node)%values)
            if(any(a>huge(1._dp)-max(b,0._dp)).or.any(a< -huge(1._dp)-min(b,0._dp))) then
              valid=.false.;exit
            end if
            a=a+b
          end associate
        end do
        call agree()
        if(.not.agreed) return
      end if
      drops=trial;masses=next_mass;gas=next_gas;temp=next_temp;accepted_id=trial_id
      ledger=ledger+trial_ledger
      if(allocated(impact_exterior)) then
        exterior(:5)=exterior(:5)+impact_exterior(:5)
        do n=1,nc
          exterior(5+components(n)%vapor_index)=exterior(5+components(n)%vapor_index)+impact_exterior(5+n)
        end do
      end if
      if (present(wall_films)) then
        patches=trial_patches;patch_mass=trial_patch_mass;wall_ledger=wall_ledger+impact_wall
      end if
      if(present(wall_films).or.allocated(impact_exterior)) then
        kept=count(drops%mass>0);survivors=pack(drops,drops%mass>0)
        allocate(kept_mass(nc,kept));kept=0
        do p=1,size(drops)
          if (drops(p)%mass<=0) cycle
          kept=kept+1;kept_mass(:,kept)=masses(:,p)
        end do
        call move_alloc(survivors,drops);call move_alloc(kept_mass,masses)
      end if
      elapsed=elapsed+dt;count_steps=count_steps+1
    end do
    np=size(drops)
    kept=count(drops%mass>0);allocate(survivors(kept),kept_mass(nc,kept));kept=0
    do p=1,np
      if (drops(p)%mass<=0) cycle
      kept=kept+1;survivors(kept)=drops(p);kept_mass(:,kept)=masses(:,p)
    end do
    if(present(eb_geometry)) then
      call migrate_eb_tree_spray_parcels(comm,tree,owners,nodal_phi,survivors,lower,spacing, &
        options%maximum_parcels,valid,kept_mass,maximum_level=maximum_level)
    else
      call migrate_amr_tree_spray_parcels(comm,tree,owners,survivors,lower,spacing,periodic, &
        options%maximum_parcels,valid,kept_mass,maximum_level=maximum_level)
    end if
    if (.not.valid) return
    if (present(wall_films)) then
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing, &
        film_settings%maximum_films,valid,nodal_phi,maximum_level)
      if (.not.valid) return
    end if
    call MPI_Allreduce(ledger,global_ledger,size(ledger),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    global_ledger(:5)=global_ledger(:5)-wall_ledger
    if (.not.all(ieee_is_finite(global_ledger))) return
    if(present(gas_change)) call move_alloc(integrated_change%fields,gas_change%fields)
    state=gas;temperature=temp
    call move_alloc(survivors,parcels);call move_alloc(kept_mass,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(last_id)) last_id=accepted_id
    if (present(boundary_change)) boundary_change=global_ledger
    if(present(eb_outflow_change)) eb_outflow_change=exterior
    substeps=count_steps;ok=.true.
  contains
    subroutine prepare_impacts()
      integer :: q,parent,lo,hi,mid
      if (allocated(impact_previous)) deallocate(impact_previous,impact_pressure)
      allocate(impact_previous(3,size(trial)),impact_pressure(size(trial)))
      valid=.true.
      do q=1,size(trial)
        parent=q
        if (q>np) then
          lo=1;hi=np;parent=0
          do while(lo<=hi)
            mid=lo+(hi-lo)/2
            if (drops(mid)%id==trial(q)%parent_id) then
              parent=mid;exit
            else if (drops(mid)%id<trial(q)%parent_id) then
              lo=mid+1
            else
              hi=mid-1
            end if
          end do
          if (parent==0) then
            valid=.false.;return
          end if
        end if
        impact_previous(:,q)=positions(:,parent);impact_pressure(q)=1
        if (drops(parent)%mass>0) impact_pressure(q)=primitive(5,parent)
      end do
    end subroutine
    subroutine apply_boundaries()
      integer :: q,c
      valid=.true.
      do q=1,size(trial)
        if (trial(q)%mass<=0) cycle
        y=next_mass(:,q)/trial(q)%mass
        call liquid_mixture_properties(components,y,trial(q)%temperature,rho,cp,e,valid)
        if (.not.valid) return
        effective=components(1);effective%density=rho;effective%density_model='constant'
        effective%density_coefficients=0;effective%density_temperature_range=0
        effective%heat_capacity=cp;effective%heat_capacity_model='constant'
        effective%heat_capacity_coefficients=0;effective%heat_capacity_temperature_range=0
        effective%reference_temperature=trial(q)%temperature;effective%reference_energy=e
        effective%boiling_model='reject';effective%critical_temperature=0
        effective%normal_boiling_temperature=0
        call apply_parcel_boundaries(trial(q),effective,lower,lower+length,faces,exchange,valid)
        if (.not.valid) return
        trial_ledger(:5)=trial_ledger(:5)+exchange
        do c=1,nc
          trial_ledger(5+components(c)%vapor_index)=trial_ledger(5+components(c)%vapor_index)+exchange(1)*y(c)
        end do
        if (trial(q)%mass<=0) next_mass(:,q)=0
      end do
      valid=all(ieee_is_finite(trial_ledger))
    end subroutine
    subroutine agree()
      call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) agreed=.false.
    end subroutine
  end subroutine
end module
