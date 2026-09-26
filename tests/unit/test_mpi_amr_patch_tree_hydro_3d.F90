program test_mpi_amr_patch_tree_hydro_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d, &
    compute_sdc_amr_patch_tree_hydro_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d,compute_sdc_mpi_amr_tree_hydro_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_ghost_3d_mod, only: amr_tree_inflow_faces_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  real(dp) :: pressure(6),mach(6),pin(6)
  real(dp), allocatable :: serial_coefficient,mpi_coefficient
  real(dp) :: coefficient
  character(len=32) :: faces(6),reconstruction,arg,solver
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth
  logical :: ok,ppm,inflow_mask(6)
  call get_command_argument(1,arg)
  ppm=trim(arg)=='ppm'
  solver="hllc"
  if (ppm.or.trim(arg)=='viscosity') solver="pelec"
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,3,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  faces='inflow'
  call amr_tree_inflow_faces_3d(tree,1,faces,inflow_mask,ok)
  call require(ok.and.all(inflow_mask),'root touches all physical faces')
  call amr_tree_inflow_faces_3d(tree,2,faces,inflow_mask,ok)
  call require(ok.and.all(inflow_mask.eqv.[.true.,.false.,.true.,.false.,.true.,.false.]), &
    'child touches only lower physical faces')
  call amr_tree_inflow_faces_3d(tree,3,faces,inflow_mask,ok)
  call require(ok.and..not.any(inflow_mask),'nested interior patch has no physical inflow')

  allocate(global%fields(3),reference%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(global%fields(node)%values(12,n(1),n(2),n(3)),reference%fields(node)%values(1,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=101325*(1+.02_dp*node)
          q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
          q(6)=.025_dp+.003_dp*i;q(9)=.22_dp;q(12)=1-q(6)-q(9)
          call reactive_primitive_to_conserved(species,q,global%fields(node)%values(:,i,j,k),exact_t,sound,ok)
          call require(ok,'initial conserved state')
        end do
      end do
    end do
    call recover_reactive_temperatures_3d(species,global%fields(node)%values, &
      spread(spread(spread(800.0_dp,1,n(1)),2,n(2)),3,n(3)),n(1),n(2),n(3),reference%fields(node)%values(1,:,:,:),ok)
    call require(ok,'serial recovery oracle')
  end do
  spacing=[.001_dp,.0015_dp,.002_dp]
  call compute_amr_patch_tree_cfl_3d(tree,species,global,reference,spacing,.15_dp,expected_dt,ok)
  call require(ok,'serial CFL')
  inlet=spread(global%fields(1)%values(:,1,1,1),2,6)
  inlet_t=reference%fields(1)%values(1,1,1,1)
  if(trim(arg)=='sdc') then
    call check_sdc()
    if(rank==0) print *, 'PASS: MPI SDC hydro stage'
    call MPI_Finalize(ierr)
    stop
  end if
  do mode=1,8
    if (ppm.and.(mode==2.or.mode==6)) cycle
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='slip_wall'
    end if
    if (mode==4) faces(1)='inflow'
    pressure=0;mach=-1;pin=0
    if (mode>=5) pressure([1,6])=102325
    if(mode>=7) then
      pressure=0;faces(1)='inflow';mach(1)=.2_dp
      if(mode==8) pin(1)=120000
    end if
    if (ppm) reconstruction='characteristic_ppm'
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    coefficient=0
    if(ppm.or.trim(arg)=='viscosity') coefficient=.037_dp
    if(allocated(serial_coefficient)) deallocate(serial_coefficient)
    if(allocated(mpi_coefficient)) deallocate(mpi_coefficient)
    if(trim(arg)=='default') then
      if(mod(mode,2)==0) then
        serial_coefficient=.1_dp
      else
        mpi_coefficient=.1_dp
      end if
    else
      serial_coefficient=coefficient;mpi_coefficient=coefficient
    end if
    oracle=global;oracle_t=reference;dt=.1_dp*expected_dt
    call advance_amr_patch_tree_hydro_3d(tree,species,oracle,oracle_t,spacing,dt,faces,reconstruction,'mc',solver, &
      serial_budget,serial_largest,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin,difmag=serial_coefficient)
    call require(ok,'serial hydro oracle')
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc',solver, &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin,difmag=mpi_coefficient)
    call require(ok,'distributed hydro')
    do node=1,3
      call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owned state only')
      call require(allocated(temperature%fields(node)%values).eqv.(owners(node)==rank),'owned temperature only')
      if (owners(node)/=rank) cycle
      call require(all(local%fields(node)%values==oracle%fields(node)%values),'exact serial conserved state')
      call require(all(temperature%fields(node)%values==oracle_t%fields(node)%values),'exact serial temperature')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest,'exact boundary and reflux diagnostics')
    saved=local;saved_t=temperature
    if(ppm.or.trim(arg)=='viscosity') then
      coefficient=.037_dp
      if(rank==nranks-1) coefficient=.05_dp
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc',solver, &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin,difmag=coefficient)
      call require(.not.ok,'mismatched coefficient rejected')
      call unchanged()
      coefficient=.037_dp
      if(rank==nranks-1) coefficient=-1._dp
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc',solver, &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin,difmag=coefficient)
      call require(.not.ok,'invalid coefficient rejected')
      call unchanged()
      coefficient=.037_dp
    end if
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc','bad', &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'flux solver failure rejects collectively')
    call unchanged()
    if (rank==nranks-1) dt=2*dt
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc',solver, &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'mismatched timestep rejected')
    call unchanged()
    if (mode>=5.and.mode<=6) then
      dt=.1_dp*expected_dt
      if (rank==nranks-1) pressure(1)=103325
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc',solver, &
        budget,largest,MPI_COMM_WORLD,ok,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'mismatched pressure rejected')
      call unchanged()
      pressure(1)=102325
      if (rank==nranks-1) pressure(3)=102325
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc',solver, &
        budget,largest,MPI_COMM_WORLD,ok,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'rank-local wall pressure rejected')
      call unchanged()
    end if
    if(mode>=7) then
      dt=.1_dp*expected_dt
      if(rank==nranks-1) mach(1)=.3_dp
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
        reconstruction,'mc',solver,budget,largest,MPI_COMM_WORLD,ok,inflow_temperature=inlet_t, &
        inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'rank-disagreeing Mach rejected')
      call unchanged()
    end if
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: distributed tree hydro, physical boundaries, exact serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_sdc()
    type(amr_tree_state_3d) :: source,owned_source,rates,serial_rates,prior,bad_state
    type(amr_tree_flux_3d), allocatable :: flux(:),serial_flux(:),prior_flux(:)
    type(amr_tree_state_3d) :: start,finish,owned_start,owned_end,source_end,owned_end_source
    real(dp) :: fractions(2)
    integer :: peer,boundary_case,n(3)
    owners=mod([0,1,0],nranks);dt=.01_dp*expected_dt
    source=global
    allocate(local%fields(3),temperature%fields(3),owned_source%fields(3))
    do peer=1,3
      source%fields(peer)%values=0
      source%fields(peer)%values(2,:,:,:)=20000*peer
      source%fields(peer)%values(3,:,:,:)=-10000*peer
      if(owners(peer)/=rank) cycle
      local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
      owned_source%fields(peer)=source%fields(peer)
    end do
    do boundary_case=1,2
      faces='periodic'
      if(boundary_case==2) faces='outflow'
      call compute_sdc_amr_patch_tree_hydro_3d(tree,species,global,reference,source,spacing,dt,faces, &
        serial_rates,serial_flux,ok)
      call require(ok,'serial SDC source-aware hydro')
      call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,owned_source,spacing,dt,faces, &
        rates,flux,MPI_COMM_WORLD,ok)
      call require(ok,'MPI SDC source-aware hydro')
      do peer=1,3
        call require(allocated(rates%fields(peer)%values).eqv.(owners(peer)==rank),'owned hydro rate only')
        call require(allocated(flux(peer)%x).eqv.(owners(peer)==rank),'owned hydro flux only')
        if(owners(peer)/=rank) cycle
        call require(all(rates%fields(peer)%values==serial_rates%fields(peer)%values),'serial/MPI SDC rate exact')
        call require(all(flux(peer)%x==serial_flux(peer)%x).and.all(flux(peer)%y==serial_flux(peer)%y).and. &
          all(flux(peer)%z==serial_flux(peer)%z),'serial/MPI SDC flux exact')
        call require(all(local%fields(peer)%values==global%fields(peer)%values),'hydro stage input state unchanged')
        call require(all(temperature%fields(peer)%values==reference%fields(peer)%values),'hydro stage input T unchanged')
      end do
    end do
    prior=rates;prior_flux=flux;bad_state=local
    if(rank==owners(3)) bad_state%fields(3)%values(5,1,1,1)=-1e20_dp
    call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,bad_state,temperature,owned_source,spacing,dt,faces, &
      rates,flux,MPI_COMM_WORLD,ok)
    call require(.not.ok,'one-owner hydro failure rejected collectively')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'MPI SDC rate rollback')
      call require(all(flux(peer)%x==prior_flux(peer)%x).and.all(flux(peer)%y==prior_flux(peer)%y).and. &
        all(flux(peer)%z==prior_flux(peer)%z),'MPI SDC flux rollback')
    end do
    allocate(start%fields(3),owned_start%fields(3),owned_end%fields(3),owned_end_source%fields(3))
    do peer=1,3
      n=tree%nodes(peer)%shape
      allocate(start%fields(peer)%values(size(q)+1,n(1),n(2),n(3)))
      start%fields(peer)%values(:size(q),:,:,:)=global%fields(peer)%values
      start%fields(peer)%values(size(q)+1,:,:,:)=reference%fields(peer)%values(1,:,:,:)
    end do
    finish=start;source_end=source
    do peer=1,2
      ! Change density and energy together, retaining temperature and velocity.
      finish%fields(peer)%values(:size(q),:,:,:)=1.02_dp*start%fields(peer)%values(:size(q),:,:,:)
      source_end%fields(peer)%values=-2*source%fields(peer)%values
      if(owners(peer)/=rank) cycle
      owned_start%fields(peer)=start%fields(peer);owned_end%fields(peer)=finish%fields(peer)
      owned_end_source%fields(peer)=source_end%fields(peer)
    end do
    faces='periodic';fractions=[.25_dp,.5_dp]
    call compute_sdc_amr_patch_tree_hydro_3d(tree,species,global,reference,source,spacing,dt,faces, &
      serial_rates,serial_flux,ok,target_depth=2,coarse_start=start,coarse_end=finish,coarse_fractions=fractions, &
      source_start=source,source_end=source_end,source_fractions=[.75_dp,.5_dp])
    call require(ok,'serial selected hydro with temporal donors')
    call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,owned_source,spacing,dt,faces, &
      rates,flux,MPI_COMM_WORLD,ok,target_depth=2,coarse_start=owned_start,coarse_end=owned_end, &
      coarse_fractions=fractions,source_start=owned_source,source_end=owned_end_source,source_fractions=[.75_dp,.5_dp])
    call require(ok,'MPI selected hydro with temporal donors')
    do peer=1,3
      call require(allocated(flux(peer)%x).eqv.(owners(peer)==rank.and.peer==3),'selected owned flux only')
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==serial_rates%fields(peer)%values),'temporal donor rates exact')
      if(peer/=3) cycle
      call require(all(flux(peer)%x==serial_flux(peer)%x).and.all(flux(peer)%y==serial_flux(peer)%y).and. &
        all(flux(peer)%z==serial_flux(peer)%z),'temporal donor flux exact')
    end do
    prior=rates;prior_flux=flux
    if(rank==nranks-1) fractions(1)=.3_dp
    call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,owned_source,spacing,dt,faces, &
      rates,flux,MPI_COMM_WORLD,ok,target_depth=2,coarse_start=owned_start,coarse_end=owned_end, &
      coarse_fractions=fractions,source_start=owned_source,source_end=owned_end_source,source_fractions=[.75_dp,.5_dp])
    call require(.not.ok,'idle-rank temporal fraction mismatch rejected')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'temporal failure rate rollback')
      if(peer/=3) cycle
      call require(all(flux(peer)%x==prior_flux(peer)%x),'temporal failure flux rollback')
    end do
  end subroutine check_sdc
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest,'diagnostic rollback')
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
