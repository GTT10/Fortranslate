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
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  implicit none
  private
  public :: advance_mpi_amr_mixture_coupling
contains
  ! Regular AMR: common particle substeps and whole-interval rollback.
  ! Species, transport, liquids and options must be identical on every rank.
  subroutine advance_mpi_amr_mixture_coupling(comm,species,transport,components,options,state,temperature, &
      tree,owners,parcels,component_mass,lower,spacing,interval,substeps,ok,last_id,boundary_faces,boundary_change, &
      wall_films,film_mass,film_settings)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: components(:)
    type(spray_options), intent(in) :: options
    type(amr_tree_state_3d), intent(inout) :: state,temperature
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
    type(spray_wall_film), allocatable :: patches(:),trial_patches(:)
    real(dp), allocatable :: patch_mass(:,:),trial_patch_mass(:,:),impact_previous(:,:),impact_pressure(:)
    real(dp) :: mw(size(components)),impact_wall(5),wall_ledger(5)
    logical, allocatable :: hit(:)
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
    if (present(boundary_faces)) faces=boundary_faces
    valid=allocated(parcels).and.allocated(component_mass).and.valid_spray_options(options)
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
    valid=valid.and.(present(wall_films).eqv.present(film_mass))
    valid=valid.and.(present(wall_films).eqv.present(film_settings))
    if (present(wall_films)) valid=valid.and.allocated(wall_films).and.present(last_id).and.present(boundary_faces)
    if (present(film_mass)) valid=valid.and.allocated(film_mass)
    if (present(film_settings)) valid=valid.and.valid_wall_options(film_settings,faces)
    call agree()
    if (.not.agreed) return
    film_flag=merge(1,0,present(wall_films))
    call MPI_Allreduce(film_flag,low_flag,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(film_flag,high_flag,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.low_flag/=high_flag) return
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
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing,film_settings%maximum_films,valid)
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
      lower,spacing,periodic,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree()
    if (.not.agreed) return
    length=spacing*real(tree%nodes(1)%shape,dp)
    minimum_spacing=minval(length/real(domain(:,size(domain,2)),dp))
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
        lower,spacing,periodic,samples,valid)
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
        if (invtau>0) local_dt=min(local_dt,options%max_fractional_change/invtau)
        tdot=heat*(tg(p)-drops(p)%temperature)
        do n=1,nc
          if (rates(n)>0.and.y(n)>mixture_trace_fraction) &
            local_dt=min(local_dt,options%max_fractional_change*masses(n,p)/rates(n))
          tdot=tdot-rates(n)*(hv(n)-liquid_energy(components(n),drops(p)%temperature))
        end do
        tdot=tdot/(drops(p)%mass*cp)
        if (abs(tdot)>0) local_dt=min(local_dt,min(2.0_dp,options%max_fractional_change*drops(p)%temperature)/abs(tdot))
        speed=max(maxval(abs(drops(p)%velocity)),maxval(abs(primitive(2:4,p))))
        if (speed>0) local_dt=min(local_dt,options%particle_cfl*minimum_spacing/speed)
      end do
      call agree()
      if (.not.agreed) return
      call MPI_Allreduce(local_dt,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      accepted=.false.
      do attempt=1,30
        if (dt<=0.or.elapsed+dt<=elapsed) exit
        trial=drops;next_mass=masses;change=0;valid=.true.;trial_id=accepted_id;reynolds=0;trial_ledger=0
        do p=1,np
          if (drops(p)%mass<=0) cycle
          call advance_mixture_frozen_gas(species,transport,primitive(:,p),tg(p),trial(p),components, &
            next_mass(:,p),options,dt,change(:,p),valid)
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
            call build_mixture_film(species,transport,primitive(:,p),tg(p),trial(p),components,y,options, &
              film,rates,hv,valid)
            if (.not.valid) exit
            if (options%khrt_children) then
              reynolds(p)=film%density*norm2(primitive(2:4,p)-trial(p)%velocity)* &
                parcel_diameter(trial(p),effective)/film%viscosity
            else
              call advance_spray_breakup(trial(p),effective,options,primitive(1,p),primitive(2:4,p), &
                film%density,film%viscosity,dt,valid)
              if (.not.valid) exit
              next_mass(:,p)=next_mass(:,p)*(trial(p)%mass/old_mass)
            end if
          end if
          if (.not.present(boundary_faces)) trial(p)%position=lower+modulo(trial(p)%position-lower,length)
        end do
        call agree()
        if (agreed.and.options%khrt_children) then
          call advance_mpi_khrt_shedding(comm,components,options,trial,next_mass,primitive(1,:), &
            primitive(2:4,:),reynolds,dt,trial_id,child_energy,valid)
          if (valid) valid=size(trial)<=125000
          call agree()
          if (agreed) change(5,:)=change(5,:)+child_energy
        end if
        impact_wall=0
        if (agreed.and.present(wall_films)) then
          trial_patches=patches;trial_patch_mass=patch_mass
          call prepare_impacts()
          call agree()
          if (agreed) then
            call impact_mpi_planar_walls(comm,components,options,film_settings,faces,mw,impact_pressure,impact_previous, &
              lower,lower+length,trial,next_mass,trial_patches,trial_patch_mass,trial_id,hit,impact_wall,valid)
            if (valid) valid=size(trial)<=125000
            call agree()
          end if
        end if
        if (agreed) then
          if (present(boundary_faces)) call apply_boundaries()
          call agree()
        end if
        if (agreed) then
          next_gas=gas;next_temp=temp
          call deposit_mpi_amr_tree_spray_gas(comm,species,tree,owners,next_gas,next_temp,positions,change, &
            lower,spacing,periodic,valid)
          if (valid) then
            accepted=.true.;exit
          end if
        end if
        dt=.5_dp*dt
      end do
      if (.not.accepted) return
      drops=trial;masses=next_mass;gas=next_gas;temp=next_temp;accepted_id=trial_id
      ledger=ledger+trial_ledger
      if (present(wall_films)) then
        patches=trial_patches;patch_mass=trial_patch_mass;wall_ledger=wall_ledger+impact_wall
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
    call migrate_amr_tree_spray_parcels(comm,tree,owners,survivors,lower,spacing,periodic, &
      options%maximum_parcels,valid,kept_mass)
    if (.not.valid) return
    if (present(wall_films)) then
      call migrate_mpi_amr_wall_films(comm,tree,owners,patches,patch_mass,lower,spacing,film_settings%maximum_films,valid)
      if (.not.valid) return
    end if
    call MPI_Allreduce(ledger,global_ledger,size(ledger),MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    global_ledger(:5)=global_ledger(:5)-wall_ledger
    if (.not.all(ieee_is_finite(global_ledger))) return
    state=gas;temperature=temp
    call move_alloc(survivors,parcels);call move_alloc(kept_mass,component_mass)
    if (present(wall_films)) then
      call move_alloc(patches,wall_films);call move_alloc(patch_mass,film_mass)
    end if
    if (present(last_id)) last_id=accepted_id
    if (present(boundary_change)) boundary_change=global_ledger
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
