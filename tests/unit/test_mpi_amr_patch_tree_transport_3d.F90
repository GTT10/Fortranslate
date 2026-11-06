program test_mpi_amr_patch_tree_transport_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use amr_patch_tree_transport_3d_mod, only: advance_amr_patch_tree_transport_3d, &
    compute_amr_patch_tree_transport_timestep_3d,compute_sdc_amr_patch_tree_transport_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: advance_mpi_amr_tree_transport_3d, &
    compute_mpi_amr_tree_transport_timestep_3d,compute_sdc_mpi_amr_tree_transport_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:),bad(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  real(dp) :: pressure(6),mach(6),pin(6)
  character(len=32) :: faces(6),reconstruction,arg
  character(len=:), allocatable :: coarse_limiter
  real(dp) :: theta,serial_theta,wall_t(6),wall_v(3,6)
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth
  real(dp) :: transport_dt,serial_dt,diffusivity,serial_diffusivity
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,3,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
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
  call get_command_argument(1,arg)
  if(trim(arg)=='sdc'.or.trim(arg)=='constant_sdc') then
    call check_sdc()
    if(rank==0) print *, 'PASS: MPI SDC transport stage'
    call MPI_Finalize(ierr)
    stop
  end if
  do mode=1,8
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    wall_t=1100;wall_v=0
    if (allocated(coarse_limiter)) deallocate(coarse_limiter)
    if (mod(mode,2)==0) coarse_limiter='mc'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
    end if
    if (mode==4) faces(1)='inflow'
    pressure=0;mach=-1;pin=0
    if (mode>=5) pressure([1,6])=102325
    if(mode>=7) then
      pressure=0;faces(1)='inflow';mach(1)=.2_dp
      if(mode==8) pin(1)=120000
    end if
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,global,reference,spacing,.15_dp, &
      .true.,.true.,.true.,serial_dt,serial_diffusivity,ok,faces,wall_t)
    call require(ok,'serial diffusion and wall timestep')
    call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,local,temperature,spacing,.15_dp, &
      .true.,.true.,.true.,transport_dt,diffusivity,MPI_COMM_WORLD,ok,faces,wall_t)
    call require(ok.and.transport_dt==serial_dt.and.diffusivity==serial_diffusivity,'exact transport timestep bounds')
    oracle=global;oracle_t=reference;dt=1e-9_dp
    call advance_amr_patch_tree_transport_3d(tree,species,transport,oracle,oracle_t,spacing,dt,faces, &
      .true.,.true.,.true.,.true.,serial_budget,serial_largest,serial_theta,ok,inlet,inlet_t,wall_t,wall_v,coarse_limiter, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'serial transport RK2 oracle')
    call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
      .true.,.true.,.true.,.true.,budget,largest,theta,MPI_COMM_WORLD,ok,inlet,inlet_t,wall_t,wall_v,coarse_limiter, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'distributed transport RK2')
    do node=1,3
      call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owned state only')
      call require(allocated(temperature%fields(node)%values).eqv.(owners(node)==rank),'owned temperature only')
      if (owners(node)/=rank) cycle
      call require(all(local%fields(node)%values==oracle%fields(node)%values),'exact serial conserved state')
      call require(all(temperature%fields(node)%values==oracle_t%fields(node)%values),'exact serial temperature')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest,'exact boundary and reflux diagnostics')
    call require(theta==serial_theta,'exact global transport limiter')
    saved=local;saved_t=temperature
    bad=transport
    if (rank==owners(3)) bad(1)%diameter=-1
    call advance_mpi_amr_tree_transport_3d(tree,owners,species,bad,local,temperature,spacing,dt,faces, &
      .true.,.true.,.true.,.true.,budget,largest,theta,MPI_COMM_WORLD,ok,inlet,inlet_t,wall_t,wall_v,coarse_limiter, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'one owner transport flux failure rejects collectively')
    call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,bad,local,temperature,spacing,.15_dp, &
      .true.,.true.,.true.,transport_dt,diffusivity,MPI_COMM_WORLD,ok,faces,wall_t)
    call require(.not.ok.and.transport_dt==serial_dt.and.diffusivity==serial_diffusivity,'timestep rollback')
    call unchanged()
    if (rank==nranks-1) dt=2*dt
    call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
      .true.,.true.,.true.,.true.,budget,largest,theta,MPI_COMM_WORLD,ok,inlet,inlet_t,wall_t,wall_v,coarse_limiter, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'mismatched timestep rejected')
    call unchanged()
    if (mode>=5.and.mode<=6) then
      dt=1e-9_dp
      if (rank==nranks-1) pressure(1)=103325
      call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
        .true.,.true.,.true.,.true.,budget,largest,theta,MPI_COMM_WORLD,ok, &
        wall_temperature=wall_t,wall_velocity=wall_v,outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'mismatched pressure rejected')
      call unchanged()
    end if
    if(mode>=7) then
      if(rank==nranks-1) pin(1)=130000
      ! Restore the common step after the timestep-disagreement guard.
      dt=.1_dp*expected_dt
      call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
        .true.,.true.,.true.,.true.,budget,largest,theta,MPI_COMM_WORLD,ok,inflow_temperature=inlet_t, &
        inflow_mach=mach,inflow_pressure=pin,wall_temperature=wall_t,wall_velocity=wall_v)
      call require(.not.ok,'rank-disagreeing Mach inlet pressure rejected')
      call unchanged()
    end if
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: owner transport RK2, PCM/PLM, moving isothermal walls, serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_sdc()
    type(amr_tree_state_3d) :: rates,serial_rates,prior
    type(amr_tree_flux_3d), allocatable :: flux(:),serial_flux(:),prior_flux(:)
    type(amr_tree_state_3d) :: start,finish,owned_start,owned_end,global_held,local_held
    real(dp) :: fractions(2)
    real(dp), allocatable :: coefficients(:)
    integer :: peer,boundary_case,n(3)
    logical :: thermal_mask(6),bad_mask(6)
    real(dp) :: thermal_dt,serial_thermal_dt,thermal_diff,serial_thermal_diff
    real(dp) :: retained_theta
    if(trim(arg)=="constant_sdc") coefficients=[1e-5_dp,.02_dp,3e-5_dp,2e-5_dp]
    owners=mod([0,1,0],nranks);dt=.01_dp*expected_dt
    allocate(local%fields(3),temperature%fields(3))
    do peer=1,3
      if(owners(peer)/=rank) cycle
      local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
    end do
    global_held=reference;local_held=temperature
    do peer=1,3
      global_held%fields(peer)%values(1,:,:,:)=1.001_dp*global%fields(peer)%values(5,:,:,:)
      if(owners(peer)==rank) local_held%fields(peer)=global_held%fields(peer)
    enddo
    do boundary_case=0,3
      faces='periodic';wall_t=1100;wall_v=0;thermal_mask=.false.
      if(boundary_case==0) then
        faces='outflow';faces(3:4)='slip_wall';thermal_mask(3:4)=.true.
      endif
      if(boundary_case==2) then
        faces='outflow';faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
      end if
      call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,global,reference,spacing,dt,faces, &
        [.true.,.true.,.true.,.true.],serial_rates,serial_flux,serial_theta,ok, &
        wall_temperature=wall_t,wall_velocity=wall_v, &
        constant_coefficients=coefficients,internal_energy=global_held, &
        total_energy_primitives=boundary_case==3,isothermal_slip_faces=thermal_mask)
      call require(ok,'serial SDC transport stage')
      call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
        [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        constant_coefficients=coefficients,internal_energy=local_held, &
        total_energy_primitives=boundary_case==3,isothermal_slip_faces=thermal_mask)
      call require(ok.and.theta==serial_theta,'MPI SDC transport stage and limiter')
      do peer=1,3
        call require(allocated(rates%fields(peer)%values).eqv.(owners(peer)==rank),'owned transport rate only')
        call require(allocated(flux(peer)%x).eqv.(owners(peer)==rank),'owned transport flux only')
        if(owners(peer)/=rank) cycle
        call require(all(rates%fields(peer)%values==serial_rates%fields(peer)%values),'serial/MPI transport rate exact')
        call require(all(flux(peer)%x==serial_flux(peer)%x).and.all(flux(peer)%y==serial_flux(peer)%y).and. &
          all(flux(peer)%z==serial_flux(peer)%z),'serial/MPI transport flux exact')
        call require(all(local%fields(peer)%values==global%fields(peer)%values),'transport input state unchanged')
        call require(all(temperature%fields(peer)%values==reference%fields(peer)%values),'transport input T unchanged')
      end do
      if(boundary_case==0) then
        call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,global,reference,spacing, &
          .15_dp,.true.,.true.,.true.,serial_thermal_dt,serial_thermal_diff,ok,faces,wall_t, &
          constant_coefficients=coefficients,isothermal_slip_faces=thermal_mask)
        call require(ok,'serial thermal slip timestep')
        call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,local,temperature,spacing, &
          .15_dp,.true.,.true.,.true.,thermal_dt,thermal_diff,MPI_COMM_WORLD,ok,faces,wall_t, &
          constant_coefficients=coefficients,isothermal_slip_faces=thermal_mask)
        call require(ok.and.thermal_dt==serial_thermal_dt.and.thermal_diff==serial_thermal_diff, &
          'thermal slip timestep serial/MPI exact')
        prior=rates;prior_flux=flux;retained_theta=theta;bad_mask=thermal_mask
        if(rank==0) bad_mask(1)=.true.
        call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
          [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok,wall_temperature=wall_t, &
          isothermal_slip_faces=bad_mask)
        call require(.not.ok.and.theta==retained_theta,'rank thermal mask mismatch rollback')
        do peer=1,3
          if(owners(peer)/=rank) cycle
          call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'thermal mask rate rollback')
          call require(all(flux(peer)%x==prior_flux(peer)%x).and.all(flux(peer)%y==prior_flux(peer)%y).and. &
            all(flux(peer)%z==prior_flux(peer)%z),'thermal mask flux rollback')
        enddo
        thermal_dt=-17;thermal_diff=-19
        call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,local,temperature,spacing, &
          .15_dp,.true.,.true.,.true.,thermal_dt,thermal_diff,MPI_COMM_WORLD,ok,faces,wall_t, &
          isothermal_slip_faces=bad_mask)
        call require(.not.ok.and.thermal_dt==-17.and.thermal_diff==-19,'rank thermal timestep mismatch rollback')
      endif
    end do
    prior=rates;prior_flux=flux;retained_theta=theta;bad=transport
    call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok, &
      internal_energy=local_held,total_energy_primitives=rank==0)
    call require(.not.ok.and.theta==retained_theta,'rank-inconsistent primitive energy mode rejected')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'energy mode mismatch rate rollback')
    enddo
    if(allocated(coefficients)) then
      if(rank==nranks-1) coefficients(3)=2*coefficients(3)
      call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
        [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
        constant_coefficients=coefficients)
      call require(.not.ok.and.theta==retained_theta,'idle-rank constant coefficient mismatch rejected')
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'coefficient mismatch rate rollback')
        call require(all(flux(peer)%x==prior_flux(peer)%x).and.all(flux(peer)%y==prior_flux(peer)%y).and. &
          all(flux(peer)%z==prior_flux(peer)%z),'coefficient mismatch flux rollback')
      enddo
      coefficients=[1e-5_dp,.02_dp,3e-5_dp,2e-5_dp]
    endif
    if(rank==owners(3)) bad(1)%diameter=-1
    call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,bad,local,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok,wall_temperature=wall_t,wall_velocity=wall_v)
    call require(.not.ok.and.theta==retained_theta,'one-owner transport failure rejected collectively')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'MPI transport rate rollback')
      call require(all(flux(peer)%x==prior_flux(peer)%x).and.all(flux(peer)%y==prior_flux(peer)%y).and. &
        all(flux(peer)%z==prior_flux(peer)%z),'MPI transport flux rollback')
    end do
    allocate(start%fields(3),owned_start%fields(3),owned_end%fields(3))
    do peer=1,3
      n=tree%nodes(peer)%shape
      allocate(start%fields(peer)%values(size(q)+1,n(1),n(2),n(3)))
      start%fields(peer)%values(:size(q),:,:,:)=global%fields(peer)%values
      start%fields(peer)%values(size(q)+1,:,:,:)=reference%fields(peer)%values(1,:,:,:)
    end do
    finish=start
    do peer=1,2
      finish%fields(peer)%values(size(q)+1,:,:,:)=1.02_dp*start%fields(peer)%values(size(q)+1,:,:,:)
      if(owners(peer)/=rank) cycle
      owned_start%fields(peer)=start%fields(peer);owned_end%fields(peer)=finish%fields(peer)
    end do
    faces='periodic';fractions=[.25_dp,.5_dp]
    call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,global,reference,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],serial_rates,serial_flux,serial_theta,ok,target_depth=2, &
      coarse_start=start,coarse_end=finish,coarse_fractions=fractions, &
        constant_coefficients=coefficients)
    call require(ok,'serial level transport temporal donors')
    call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok,target_depth=2, &
      coarse_start=owned_start,coarse_end=owned_end,coarse_fractions=fractions, &
        constant_coefficients=coefficients)
    call require(ok.and.theta==serial_theta,'MPI level transport temporal donors and limiter')
    do peer=1,3
      call require(allocated(flux(peer)%x).eqv.(owners(peer)==rank.and.peer==3),'selected owned transport flux only')
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==serial_rates%fields(peer)%values),'temporal transport rates exact')
      if(peer/=3) cycle
      call require(all(flux(peer)%x==serial_flux(peer)%x).and.all(flux(peer)%y==serial_flux(peer)%y).and. &
        all(flux(peer)%z==serial_flux(peer)%z),'temporal transport flux exact')
    end do
    prior=rates;prior_flux=flux;retained_theta=theta
    if(rank==nranks-1) fractions(1)=.3_dp
    call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,local,temperature,spacing,dt,faces, &
      [.true.,.true.,.true.,.true.],rates,flux,theta,MPI_COMM_WORLD,ok,target_depth=2, &
      coarse_start=owned_start,coarse_end=owned_end,coarse_fractions=fractions, &
        constant_coefficients=coefficients)
    call require(.not.ok.and.theta==retained_theta,'idle-rank transport fraction mismatch rejected')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(rates%fields(peer)%values==prior%fields(peer)%values),'temporal transport rate rollback')
      if(peer/=3) cycle
      call require(all(flux(peer)%x==prior_flux(peer)%x),'temporal transport flux rollback')
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
    call require(theta==serial_theta,'limiter diagnostic rollback')
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
