program test_mpi_amr_patch_tree_chemistry_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use amr_patch_tree_chemistry_3d_mod, only: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d, &
    react_sdc_amr_patch_tree_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: advance_mpi_amr_tree_chemistry_3d,advance_mpi_amr_tree_strang_3d, &
    react_sdc_mpi_amr_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:),bad(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  integer :: first_mode,last_mode
  character(len=32) :: test_mode
  real(dp) :: pressure(6)
  character(len=32) :: faces(6),reconstruction,integrator
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'mechanism')
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
          q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=202650*(1+.02_dp*node)
          q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
          q(6)=.025_dp+.003_dp*i;q(7)=1e-5_dp;q(9)=.22_dp;q(12)=1-q(6)-q(7)-q(9)
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
  first_mode=1;last_mode=3
  call get_command_argument(1,test_mode)
  if (trim(test_mode)=='--sdc-only') then
    call check_sdc()
    if(rank==0) print *, 'PASS: MPI SDC owner reaction and collective rollback'
    call MPI_Finalize(ierr)
    stop
  end if
  if (trim(test_mode)=='--pressure-only') then
    first_mode=4;last_mode=4
  end if
  do mode=first_mode,last_mode
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='slip_wall'
    end if
    if (mode==3) faces(1)='inflow'
    integrator='explicit'
    if (mod(mode,2)==0) integrator='implicit'
    pressure=0
    if (mode==4) pressure([1,6])=203650
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    oracle=global;oracle_t=reference;dt=1e-9_dp
    if (mode==4) dt=1e-10_dp
    call advance_amr_patch_tree_chemistry_3d(tree,species,reactions,oracle,oracle_t,dt,1e-6_dp,1e-14_dp,ok,integrator)
    call require(ok,'serial exposed-cell chemistry')
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,local,temperature,dt,1e-6_dp,1e-14_dp, &
      MPI_COMM_WORLD,ok,integrator)
    call require(ok,'distributed exposed-cell chemistry')
    call compare_state()
    call require(any(oracle%fields(3)%values(6:12,:,:,:)/=global%fields(3)%values(6:12,:,:,:)), &
      'leaf chemistry is active')
    oracle=global;oracle_t=reference
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    call advance_amr_patch_tree_strang_3d(tree,species,reactions,oracle,oracle_t,spacing,dt,faces,reconstruction,'mc','hllc', &
      1e-6_dp,1e-14_dp,serial_budget,serial_largest,ok,integrator,inlet,inlet_t, &
      outlet_pressure=pressure)
    call require(ok,'serial reaction-hydro composition')
    call advance_mpi_amr_tree_strang_3d(tree,owners,species,reactions,local,temperature,spacing,dt,faces, &
      reconstruction,'mc','hllc',1e-6_dp,1e-14_dp,budget,largest,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t, &
      outlet_pressure=pressure)
    call require(ok,'distributed reaction-hydro composition')
    call compare_state()
    call require(all(budget==serial_budget).and.largest==serial_largest,'exact coupled diagnostics')
    saved=local;saved_t=temperature
    call advance_mpi_amr_tree_strang_3d(tree,owners,species,reactions,local,temperature,spacing,dt,faces, &
      reconstruction,'mc','bad',1e-6_dp,1e-14_dp,budget,largest,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t, &
      outlet_pressure=pressure)
    call require(.not.ok,'hydro failure after chemistry rolls back whole split')
    call unchanged()
    bad=reactions
    if (rank==owners(3)) bad(1)%reactant_stoich(1)=-1
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,bad,local,temperature,dt,1e-6_dp,1e-14_dp, &
      MPI_COMM_WORLD,ok,integrator)
    call require(.not.ok,'one owner invalid reaction rejects collectively')
    call unchanged()
    if (rank==nranks-1) dt=2*dt
    call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,local,temperature,dt,1e-6_dp,1e-14_dp, &
      MPI_COMM_WORLD,ok,integrator)
    call require(.not.ok,'inconsistent reaction interval')
    call unchanged()
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: owner chemistry and Strang composition, explicit/implicit parity and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_sdc()
    type(amr_tree_state_3d) :: rates,local_rates,old,old_t,ir,serial_ir,retained_ir
    type(amr_tree_state_3d) :: serial_aux,local_aux,retained_aux
    real(dp), allocatable :: forced(:,:,:,:)
    integer :: peer
    owners=mod([0,1,0],nranks);dt=1e-9_dp
    rates=global;serial_aux=reference
    do peer=1,3
      serial_aux%fields(peer)%values=10000
    enddo
    do peer=1,3
      rates%fields(peer)%values=.1_dp*global%fields(peer)%values
    end do
    call react_sdc_amr_patch_tree_3d(tree,species,reactions,global,reference,rates,dt,1e-6_dp,1e-14_dp, &
      oracle,oracle_t,serial_ir,ok,internal_energy=serial_aux)
    call require(ok,'serial SDC reaction reference')
    do peer=1,3
      forced=global%fields(peer)%values+dt*rates%fields(peer)%values
      call require(all(serial_aux%fields(peer)%values(1,:,:,:)== &
        forced(5,:,:,:)-.5_dp*sum(forced(2:4,:,:,:)**2,dim=1)/forced(1,:,:,:)), &
        'AMR reactor exports exact forced internal energy endpoint')
    enddo
    allocate(old%fields(3),old_t%fields(3),local_rates%fields(3),local_aux%fields(3))
    do peer=1,3
      if(owners(peer)/=rank) cycle
      old%fields(peer)=global%fields(peer);old_t%fields(peer)=reference%fields(peer)
      local_rates%fields(peer)=rates%fields(peer)
      local_aux%fields(peer)=reference%fields(peer);local_aux%fields(peer)%values=10000
    end do
    call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,old_t,local_rates,dt,1e-6_dp,1e-14_dp, &
      local,temperature,ir,MPI_COMM_WORLD,ok,internal_energy=local_aux)
    call require(ok,'owner-local SDC reaction')
    call compare_state()
    do peer=1,3
      call require(allocated(ir%fields(peer)%values).eqv.(owners(peer)==rank),'owner-local reaction estimates only')
      call require(allocated(local_aux%fields(peer)%values).eqv.(owners(peer)==rank),'reaction scalar owned only')
      if(owners(peer)/=rank) cycle
      call require(all(local_aux%fields(peer)%values==serial_aux%fields(peer)%values),'serial MPI reaction scalar exact')
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'serial SDC reaction estimates exact')
      if(tree%nodes(peer)%level==1) cycle
      local%fields(peer)%values(2,:,:,:)=local%fields(peer)%values(2,:,:,:)+1
      temperature%fields(peer)%values=temperature%fields(peer)%values+1
      ir%fields(peer)%values=ir%fields(peer)%values+1
      local_aux%fields(peer)%values=local_aux%fields(peer)%values+1
    end do
    saved=local;saved_t=temperature;retained_ir=ir;retained_aux=local_aux
    call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,old_t,local_rates,dt,1e-6_dp,1e-14_dp, &
      local,temperature,ir,MPI_COMM_WORLD,ok,target_depth=1,internal_energy=local_aux)
    call require(ok,'MPI level-selective SDC reaction')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'other-level state retained')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'other-level temperature retained')
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'other-level source retained')
      call require(all(local_aux%fields(peer)%values==retained_aux%fields(peer)%values),'other-level scalar retained')
    end do
    if(rank==owners(3)) old%fields(3)%values(1,1,1,1)=-1
    call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,old_t,local_rates,dt,1e-6_dp,1e-14_dp, &
      local,temperature,ir,MPI_COMM_WORLD,ok,internal_energy=local_aux)
    call require(.not.ok,'one-owner SDC failure rejected collectively')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'collective SDC state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'collective SDC temperature rollback')
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'collective SDC source rollback')
      call require(all(local_aux%fields(peer)%values==retained_aux%fields(peer)%values),'collective SDC scalar rollback')
      old%fields(peer)=global%fields(peer)
    end do
    if(rank==nranks-1) dt=2*dt
    call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,old_t,local_rates,dt,1e-6_dp,1e-14_dp, &
      local,temperature,ir,MPI_COMM_WORLD,ok,internal_energy=local_aux)
    call require(.not.ok,'idle-rank interval disagreement rejected')
    dt=1e-9_dp
    if(rank==nranks-1) then
      call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,old_t,local_rates,dt,1e-6_dp,1e-14_dp, &
        local,temperature,ir,MPI_COMM_WORLD,ok)
    else
      call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,old_t,local_rates,dt,1e-6_dp,1e-14_dp, &
        local,temperature,ir,MPI_COMM_WORLD,ok,internal_energy=local_aux)
    endif
    call require(.not.ok,'rank-disagreeing reactor scalar presence rejected')
  end subroutine check_sdc
  subroutine compare_state()
    integer :: peer
    do peer=1,3
      call require(allocated(local%fields(peer)%values).eqv.(owners(peer)==rank),'owned state only')
      call require(allocated(temperature%fields(peer)%values).eqv.(owners(peer)==rank),'owned temperature only')
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==oracle%fields(peer)%values),'exact serial conserved state')
      call require(all(temperature%fields(peer)%values==oracle_t%fields(peer)%values),'exact serial temperature')
    end do
  end subroutine
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
