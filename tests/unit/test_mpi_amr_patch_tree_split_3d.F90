program test_mpi_amr_patch_tree_split_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_nan
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use dual_energy_mod, only: dual_energy_policy
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    composite_integrals_amr_patch_tree_3d,average_down_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d, &
    convert_amr_patch_tree_rotating_energy_3d
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d, &
    advance_amr_patch_tree_sdc_subcycle_3d
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d,advance_mpi_amr_tree_sdc_3d, &
    advance_mpi_amr_tree_sdc_subcycle_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use amr_patch_tree_chemistry_3d_mod, only: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: advance_mpi_amr_tree_chemistry_3d,advance_mpi_amr_tree_strang_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_checkpoint_3d_mod, only: amr_tree_checkpoint_3d,read_amr_tree_checkpoint_3d,write_amr_tree_checkpoint_3d
  use mpi_amr_patch_tree_checkpoint_3d_mod, only: read_mpi_amr_tree_checkpoint_3d,write_mpi_amr_tree_checkpoint_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:),bad(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  real(dp) :: pressure(6),mach(6),pin(6)
  character(len=32) :: faces(6),reconstruction,integrator,arg,solver
  real(dp) :: coefficient,reflux_specific_energy
  logical :: diffusion_failure=.false.
  integer :: diffusion_calls=0,reactor_calls=0
  logical :: reactor_failure=.false.
  integer :: recursive_auxiliary_calls(3)=0
  integer :: auxiliary_calls=0
  integer :: ordered_mol_calls(3)=0
  logical :: auxiliary_failure=.false.,serial_auxiliary_failure=.false.
  logical :: ppm,stored_energy_failure=.false.
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth
  logical :: ok,enabled,bad_enabled
  real(dp) :: theta,serial_theta,wall_t(6),wall_v(3,6)
  call get_command_argument(1,arg)
  ppm=trim(arg)=='ppm'
  coefficient=0;solver='hllc'
  if(ppm) then
    coefficient=.037_dp;solver='pelec'
  end if
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'mechanism')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(2,2,2,2,2,2,2,2,2,3,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([2,2,2],[1,2],patches,tree,ok)
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
  if(trim(arg)=='repair_restart') then
    call check_reflux_held_restart()
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='reflux_budget') then
    call check_reflux_budget()
    if(rank==0) print *, 'PASS: serial/MPI split reflux repair and refusal, common/recursive'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='stored_energy') then
    call check_old_admission()
    call check_stored_energy()
    if(rank==0) print *, 'PASS: MPI retained auxiliary MOL common/recursive time integration'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='mol') then
    call check_mol()
    call check_ordered_mol_auxiliary()
    if(rank==0) print *, 'PASS: common/recursive regular MOL state, transport, chemistry and rollback'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='mol_conserve') then
    call check_mol_conservation()
    if(rank==0) print *, 'PASS: regular AMR nonreacting MOL composite conservation'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='rotating_mol'.or.trim(arg)=='rotating_boundaries') then
    call check_rotating_mol(trim(arg)=='rotating_boundaries')
    if(rank==0) print *, 'PASS: MPI rotating regular AMR common/recursive MOL conservation and serial parity'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='recursive_auxiliary') then
    call check_recursive_auxiliary()
    if(rank==0) print *, 'PASS: recursive regular AMR auxiliary stages, times, inventory and rollback'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='auxiliary') then
    call check_auxiliary()
    if(rank==0) print *, 'PASS: regular AMR SDC auxiliary source stages, inventory and rollback'
    call MPI_Finalize(ierr)
    stop
  endif
  if(trim(arg)=='sdc') then
    call check_sdc()
    if(rank==0) print *, 'PASS: coupled MPI SDC parity and rollback'
    call MPI_Finalize(ierr)
    stop
  end if
  do mode=1,6
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if(ppm) reconstruction='characteristic_ppm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='slip_wall'
    end if
    if (mode==3) faces(1)='inflow'
    integrator='explicit'
    if (mod(mode,2)==0) integrator='implicit'
    wall_t=1100;wall_v=0;enabled=mode/=3
    if (mod(mode,2)==0) then
      faces='outflow';faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
    end if
    pressure=0;mach=-1;pin=0
    if (mode==4) pressure([1,6])=203650
    if(mode>=5) then
      faces(1)='inflow';mach(1)=.2_dp
      if(mode==6) pin(1)=210000
    end if
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    oracle=global;oracle_t=reference;dt=1e-10_dp
    call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
      reconstruction,'mc',solver,enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      serial_budget,serial_largest,serial_theta,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6,difmag=coefficient)
    call require(ok,'serial full operator split')
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc',solver,enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6,difmag=coefficient)
    call require(ok,'distributed full operator split')
    call compare_state()
    call require(all(budget==serial_budget).and.largest==serial_largest.and.theta==serial_theta,'exact coupled diagnostics')
    saved=local;saved_t=temperature
    if(ppm) then
      if(rank==nranks-1) coefficient=.05_dp
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc',solver,enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6,difmag=coefficient)
      call require(.not.ok,'coefficient mismatch rolls back chemistry and transport')
      call unchanged()
      coefficient=.037_dp
    end if
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc','bad',enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6)
    call require(.not.ok,'hydro failure rolls back preceding chemistry and transport')
    call unchanged()
    bad_enabled=enabled
    if (rank==nranks-1) bad_enabled=.not.enabled
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc',solver,bad_enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6)
    call require(.not.ok,'inconsistent enabled operators rejected')
    call unchanged()
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: distributed reaction transport hydro split, disabled operators, serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_reflux_held_restart()
    type(amr_tree_state_3d) :: sh,mh,initial_state,initial_temp,initial_held,control,control_t,control_h
    type(amr_tree_state_3d) :: saved_local,saved_temp,saved_held,saved_sh
    type(amr_tree_flux_3d), allocatable :: sf(:),mf(:),saf(:),maf(:),save_sf(:),save_mf(:),save_saf(:),save_maf(:)
    type(amr_tree_checkpoint_3d) :: serial_snapshot,mpi_snapshot
    type(dual_energy_policy) :: policy
    integer, allocatable :: restored_owners(:)
    real(dp) :: before(12),after(12),held_before(1),held_after(1),external(12),serial_external(12)
    real(dp) :: se,me,spw,mpw,auxdiff,gasdiff,control_budget(12),control_external(12),control_max,control_theta,ce,cpw
    real(dp) :: accumulated_boundary(12),initial_inventory(12),held_residual,gas_residual,total_energy,total_pressure
    real(dp) :: total_maximum,total_minimum
    integer :: p,nn(3),ii,jj,kk,step,steps_to_run,completed
    logical :: cycling
    character(len=32) :: clock_mode,operation
    character(len=1024) :: checkpoint_path,output_path
    character(len=:), allocatable :: context
    call get_command_argument(2,clock_mode)
    call get_command_argument(3,operation)
    call get_command_argument(4,checkpoint_path)
    call get_command_argument(5,output_path)
    call require(clock_mode=='common'.or.clock_mode=='recursive','explicit restart clock')
    call require(operation=='write'.or.operation=='direct'.or.operation=='continue','explicit restart operation')
    cycling=clock_mode=='recursive'
    context='manufactured-reflux-held-v1 '//trim(clock_mode)//' eta2=1 update_E=false'//achar(10)
    policy=dual_energy_policy();policy%eta2=1;policy%update_total_from_internal=.false.
    policy%allow_negative_energy=.false.
    species=spread(species(size(species)),1,size(species))
    call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,2,patches(1),ok)
    call require(ok,'held repair child')
    call initialize_amr_patch_3d(2,2,2,2,2,2,2,2,2,2,patches(2),ok)
    call require(ok,'held repair grandchild')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
    call require(ok,'held repair tree')
    deallocate(global%fields,reference%fields)
    allocate(global%fields(3),reference%fields(3),sh%fields(3))
    owners=mod([0,1,2],nranks);spacing=1;dt=1.e-8_dp;faces='periodic'
    do p=1,3
      nn=tree%nodes(p)%shape
      allocate(global%fields(p)%values(12,nn(1),nn(2),nn(3)),reference%fields(p)%values(1,nn(1),nn(2),nn(3)), &
        sh%fields(p)%values(1,nn(1),nn(2),nn(3)))
      ! This pressure keeps NASA internal energy positive for the explicit policy.
      q=0;q(1)=1;q(5)=2.e5_dp;q(2)=.2_dp;q(12)=1
      if(p>1) then
        q(6)=.1_dp;q(12)=.9_dp
        q(5)=.98_dp*q(5)
      endif
      do kk=1,nn(3)
        do jj=1,nn(2)
          do ii=1,nn(1)
            call reactive_primitive_to_conserved(species,q,global%fields(p)%values(:,ii,jj,kk),exact_t,sound,ok)
            call require(ok,'held repair initial state')
            reference%fields(p)%values(1,ii,jj,kk)=exact_t
            sh%fields(p)%values(1,ii,jj,kk)=global%fields(p)%values(5,ii,jj,kk)-.5_dp*q(1)*sum(q(2:4)**2)
          enddo
        enddo
      enddo
    enddo
    ! Initial covered coarse values must represent the same composite field.
    call average_down_amr_patch_tree_3d(tree,global,ok)
    call require(ok,'initial gas average down')
    call average_down_amr_patch_tree_3d(tree,sh,ok)
    call require(ok,'initial held average down')
    initial_temp=reference
    do p=1,3
      nn=tree%nodes(p)%shape
      call recover_reactive_temperatures_3d(species,global%fields(p)%values,initial_temp%fields(p)%values(1,:,:,:), &
        nn(1),nn(2),nn(3),reference%fields(p)%values(1,:,:,:),ok,internal_energy_density=sh%fields(p)%values(1,:,:,:))
      call require(ok,'initial held recovery')
    enddo
    reflux_specific_energy=global%fields(1)%values(5,1,1,1)
    oracle=global;oracle_t=reference
    call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,initial_inventory,ok)
    call require(ok,'held initial inventory')
    se=0;me=0;spw=0;mpw=0;completed=0;accumulated_boundary=0;total_energy=0;total_pressure=0
    total_maximum=0;total_minimum=1
    if(operation=='continue') then
      call read_amr_tree_checkpoint_3d(trim(checkpoint_path),species,context,[4,4,4],serial_snapshot,ok, &
        allow_internal_energy=.true.,mol_energy_reset=.true.)
      call require(ok,'serial reads active repaired checkpoint')
      call read_mpi_amr_tree_checkpoint_3d(trim(checkpoint_path),species,context,[4,4,4],'weighted',2, &
        mpi_snapshot,restored_owners,MPI_COMM_WORLD,ok,allow_internal_energy=.true.,mol_energy_reset=.true.)
      call require(ok,'changed-rank reads active repaired checkpoint')
      call require(serial_snapshot%steps==1.and.serial_snapshot%time==dt.and. &
        serial_snapshot%previous_timestep==dt,'repaired checkpoint accepted clock')
      call require(all(serial_snapshot%initial==initial_inventory),'repaired checkpoint initial inventory identity')
      owners=restored_owners
      oracle=serial_snapshot%state;oracle_t=serial_snapshot%temperature;sh=serial_snapshot%internal_energy
      local=mpi_snapshot%state;temperature=mpi_snapshot%temperature;mh=mpi_snapshot%internal_energy
      se=serial_snapshot%energy_correction;me=mpi_snapshot%energy_correction
      spw=serial_snapshot%pressure_work_change;mpw=mpi_snapshot%pressure_work_change
      total_energy=se;total_pressure=spw
      total_maximum=serial_snapshot%maximum_reflux;total_minimum=serial_snapshot%minimum_theta
      completed=serial_snapshot%steps;accumulated_boundary=serial_snapshot%boundary
      call require(all(mpi_snapshot%initial==serial_snapshot%initial).and. &
        all(mpi_snapshot%boundary==serial_snapshot%boundary).and.me==se.and.mpw==spw,'restored ledger parity')
    else
      allocate(local%fields(3),temperature%fields(3),mh%fields(3))
      do p=1,3
        if(owners(p)/=rank) cycle
        local%fields(p)=oracle%fields(p);temperature%fields(p)=oracle_t%fields(p);mh%fields(p)=sh%fields(p)
      enddo
    endif
    if(rank==0) print *, 'REPAIR_RESTART ',trim(clock_mode),' ',trim(operation),' ranks=',nranks,' owners=',owners
    call compare_state()
    do p=1,3
      if(owners(p)/=rank) cycle
      call require(all(mh%fields(p)%values==sh%fields(p)%values),'initial/restored held parity')
    enddo
    steps_to_run=1
    if(operation=='direct') steps_to_run=2
    do step=1,steps_to_run
      initial_state=oracle;initial_temp=oracle_t;initial_held=sh
      call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,before,ok)
      call require(ok,'pre-step gas inventory')
      call composite_integrals_amr_patch_tree_3d(tree,sh,spacing,held_before,ok)
      call require(ok,'pre-step held inventory')
      se=0;me=0;spw=0;mpw=0
      if(cycling) then
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,subcycling=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=sf,internal_energy=sh,energy_correction=se, &
          energy_policy=policy,internal_energy_integrated_flux=saf,internal_energy_pressure_work_change=spw)
      else
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=sf,internal_energy=sh,energy_correction=se, &
          energy_policy=policy,internal_energy_integrated_flux=saf,internal_energy_pressure_work_change=spw)
      endif
      call require(ok,'active held serial full split accepted')
      if(cycling) then
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,subcycling=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=mf,internal_energy=mh,energy_correction=me, &
          energy_policy=policy,internal_energy_integrated_flux=maf,internal_energy_pressure_work_change=mpw)
      else
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=mf,internal_energy=mh,energy_correction=me, &
          energy_policy=policy,internal_energy_integrated_flux=maf,internal_energy_pressure_work_change=mpw)
      endif
      call require(ok,'active held MPI full split accepted')
      call compare_state()
      call require(all(budget==serial_budget).and.all(external==serial_external).and. &
        largest==serial_largest.and.theta==serial_theta.and.me==se.and.mpw==spw,'active held full split diagnostics')
      do p=1,3
        call require(all(oracle%fields(p)%values(6:12,:,:,:)>=0),'active held species nonnegative')
        if(owners(p)/=rank) cycle
        call require(all(sh%fields(p)%values==mh%fields(p)%values),'active held exact serial/MPI scalar')
        call require(all(sf(p)%x==mf(p)%x).and.all(sf(p)%y==mf(p)%y).and.all(sf(p)%z==mf(p)%z), &
          'active held gas flux parity')
        call require(all(saf(p)%x==maf(p)%x).and.all(saf(p)%y==maf(p)%y).and.all(saf(p)%z==maf(p)%z), &
          'active held scalar flux parity')
      enddo
      call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,after,ok)
      call require(ok,'post-step gas inventory')
      after(5)=after(5)-se
      gas_residual=maxval(abs(after-before-serial_budget-serial_external)/max(1._dp,abs(before)))
      call require(gas_residual<2.e-12_dp,'active held gas conservation with correction ledger')
      call composite_integrals_amr_patch_tree_3d(tree,sh,spacing,held_after,ok)
      call require(ok,'post-step held inventory')
      held_residual=abs(held_after(1)-held_before(1)-spw)/max(1._dp,abs(held_before(1)))
      call require(held_residual<2.e-12_dp,'active held periodic inventory with pressure-work ledger')
      ! Same thermo/flow and fine cells, but give exposed coarse cells fuel:
      ! this control removes the exhausted-species trigger without changing EOS.
      if(completed==0) then
        control=initial_state;control_t=initial_temp;control_h=initial_held
        do p=1,3
          control%fields(p)%values(12,:,:,:)=control%fields(p)%values(12,:,:,:)+ &
            control%fields(p)%values(6,:,:,:)-.1_dp*control%fields(p)%values(1,:,:,:)
          control%fields(p)%values(6,:,:,:)=.1_dp*control%fields(p)%values(1,:,:,:)
        enddo
        ce=0;cpw=0
        if(cycling) then
          call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,control,control_t,spacing,dt,faces, &
            1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
            control_budget,control_external,control_max,control_theta,ok,mol_mode=.true.,subcycling=.true., &
            evaluate_diffusion=reflux_diffusion,internal_energy=control_h,energy_correction=ce,energy_policy=policy, &
            internal_energy_pressure_work_change=cpw)
        else
          call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,control,control_t,spacing,dt,faces, &
            1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
            control_budget,control_external,control_max,control_theta,ok,mol_mode=.true., &
            evaluate_diffusion=reflux_diffusion,internal_energy=control_h,energy_correction=ce,energy_policy=policy, &
            internal_energy_pressure_work_change=cpw)
        endif
        call require(ok,'nonexhausted matched thermo control')
        gasdiff=oracle%fields(1)%values(1,1,2,2)-control%fields(1)%values(1,1,2,2)
        auxdiff=maxval(abs(sh%fields(1)%values-control_h%fields(1)%values))
        if(rank==0) print *, 'ACTIVE_REPAIR',trim(clock_mode),'rho correction=',gasdiff,'held correction=',auxdiff
        call require(gasdiff>.005_dp,'repair activation restores full density vector')
        call require(auxdiff>100*epsilon(1._dp)*maxval(abs(initial_held%fields(1)%values)), &
          'repair activation changes auxiliary held energy above roundoff')
      endif
      completed=completed+1;accumulated_boundary=accumulated_boundary+serial_budget
      total_energy=total_energy+se;total_pressure=total_pressure+spw
      total_maximum=max(total_maximum,serial_largest);total_minimum=min(total_minimum,serial_theta)
      if(rank==0) print *, 'ACCEPTED_STEP',completed,'gas residual=',gas_residual,'held residual=',held_residual
    enddo
    serial_snapshot%tree=tree;serial_snapshot%state=oracle;serial_snapshot%temperature=oracle_t
    serial_snapshot%internal_energy=sh;serial_snapshot%energy_correction=total_energy
    serial_snapshot%pressure_work_change=total_pressure
    serial_snapshot%initial=initial_inventory;serial_snapshot%boundary=accumulated_boundary
    serial_snapshot%time=completed*dt;serial_snapshot%steps=completed;serial_snapshot%previous_timestep=dt
    serial_snapshot%maximum_reflux=total_maximum;serial_snapshot%minimum_theta=total_minimum
    mpi_snapshot=serial_snapshot;mpi_snapshot%state=local;mpi_snapshot%temperature=temperature
    mpi_snapshot%internal_energy=mh
    if(operation=='write') output_path=checkpoint_path
    call write_mpi_amr_tree_checkpoint_3d(trim(output_path),species,context,mpi_snapshot,owners,local,temperature, &
      MPI_COMM_WORLD,ok,indexed=.true.,mol_energy_reset=.true.)
    call require(ok,'publish legitimately repaired MPI checkpoint')
    if(rank==0) then
      call write_amr_tree_checkpoint_3d(trim(output_path)//'.serial',species,context,serial_snapshot,ok, &
        indexed=.true.,mol_energy_reset=.true.)
      call require(ok,'publish serial authority repaired checkpoint')
    endif
    ! Failure after a successful trajectory: remove all fuel but retain gas,
    ! held state and preexisting flux/ledger outputs as rollback sentinels.
    do p=1,3
      oracle%fields(p)%values(12,:,:,:)=oracle%fields(p)%values(12,:,:,:)+oracle%fields(p)%values(6,:,:,:)
      oracle%fields(p)%values(6,:,:,:)=0
      if(owners(p)/=rank) cycle
      local%fields(p)%values=oracle%fields(p)%values
    enddo
    initial_state=oracle;initial_temp=oracle_t;saved_sh=sh
    saved_local=local;saved_temp=temperature;saved_held=mh
    save_sf=sf;save_mf=mf;save_saf=saf;save_maf=maf
    serial_budget=17;serial_external=18;serial_largest=19;serial_theta=20;se=21;spw=22
    budget=17;external=18;largest=19;theta=20;me=21;mpw=22
    if(cycling) then
      call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
        serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,subcycling=.true., &
        evaluate_diffusion=reflux_diffusion,internal_energy=sh,energy_correction=se,energy_policy=policy, &
        integrated_flux=sf,internal_energy_integrated_flux=saf,internal_energy_pressure_work_change=spw)
    else
      call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
        serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
        evaluate_diffusion=reflux_diffusion,internal_energy=sh,energy_correction=se,energy_policy=policy, &
        integrated_flux=sf,internal_energy_integrated_flux=saf,internal_energy_pressure_work_change=spw)
    endif
    call require(.not.ok,'active held serial unfunded repair refuses')
    if(cycling) then
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
        budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,subcycling=.true., &
        evaluate_diffusion=reflux_diffusion,internal_energy=mh,energy_correction=me,energy_policy=policy, &
        integrated_flux=mf,internal_energy_integrated_flux=maf,internal_energy_pressure_work_change=mpw)
    else
      call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
        budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
        evaluate_diffusion=reflux_diffusion,internal_energy=mh,energy_correction=me,energy_policy=policy, &
        integrated_flux=mf,internal_energy_integrated_flux=maf,internal_energy_pressure_work_change=mpw)
    endif
    call require(.not.ok,'active held MPI unfunded repair refuses')
    call require(all(serial_budget==17).and.all(serial_external==18).and.serial_largest==19.and.serial_theta==20.and. &
      se==21.and.spw==22,'active held serial diagnostic rollback')
    call require(all(budget==17).and.all(external==18).and.largest==19.and.theta==20.and.me==21.and.mpw==22, &
      'active held MPI diagnostic rollback')
    do p=1,3
      call require(all(oracle%fields(p)%values==initial_state%fields(p)%values).and. &
        all(oracle_t%fields(p)%values==initial_temp%fields(p)%values).and. &
        all(sh%fields(p)%values==saved_sh%fields(p)%values),'active held serial state rollback')
      call require(all(sf(p)%x==save_sf(p)%x).and.all(sf(p)%y==save_sf(p)%y).and.all(sf(p)%z==save_sf(p)%z).and. &
        all(saf(p)%x==save_saf(p)%x).and.all(saf(p)%y==save_saf(p)%y).and.all(saf(p)%z==save_saf(p)%z), &
        'active held serial flux rollback')
      if(owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==saved_local%fields(p)%values).and. &
        all(temperature%fields(p)%values==saved_temp%fields(p)%values).and. &
        all(mh%fields(p)%values==saved_held%fields(p)%values),'active held MPI state rollback')
      call require(all(mf(p)%x==save_mf(p)%x).and.all(mf(p)%y==save_mf(p)%y).and.all(mf(p)%z==save_mf(p)%z).and. &
        all(maf(p)%x==save_maf(p)%x).and.all(maf(p)%y==save_maf(p)%y).and.all(maf(p)%z==save_maf(p)%z), &
        'active held MPI flux rollback')
    enddo
    if(rank==0) print *, 'PASS_ACTIVE_HELD_REPAIR_RESTART',trim(clock_mode),trim(operation),completed
  end subroutine check_reflux_held_restart

  subroutine check_reflux_budget()
    type(amr_tree_state_3d) :: serial_initial,serial_temp_initial,mpi_initial,mpi_temp_initial
    type(amr_tree_flux_3d), allocatable :: sf(:),mf(:)
    real(dp) :: before(12),after(12),external(12),serial_external(12)
    integer :: p,nn(3),ii,jj,kk,variant
    logical :: cycling,failure
    ! Equal caloric species isolate reflux from composition-dependent hydro.
    species=spread(species(size(species)),1,size(species))
    call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,2,patches(1),ok)
    call require(ok,'repair child')
    call initialize_amr_patch_3d(2,2,2,2,2,2,2,2,2,2,patches(2),ok)
    call require(ok,'repair grandchild')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
    call require(ok,'repair tree')
    deallocate(global%fields,reference%fields)
    allocate(global%fields(3),reference%fields(3))
    do p=1,3
      nn=tree%nodes(p)%shape
      allocate(global%fields(p)%values(12,nn(1),nn(2),nn(3)),reference%fields(p)%values(1,nn(1),nn(2),nn(3)))
    enddo
    owners=mod([0,1,2],nranks);spacing=1;dt=1.e-8_dp;faces='periodic'
    do variant=1,4
      cycling=mod(variant,2)==0;failure=variant>=3
      do p=1,3
        nn=tree%nodes(p)%shape
        q=0;q(1)=1;q(5)=1.e5_dp;q(12)=1
        if(p>1.and..not.failure) then
          q(6)=.1_dp;q(12)=.9_dp
        endif
        do kk=1,nn(3)
          do jj=1,nn(2)
            do ii=1,nn(1)
              call reactive_primitive_to_conserved(species,q,global%fields(p)%values(:,ii,jj,kk),exact_t,sound,ok)
              call require(ok,'reflux initial state')
              reference%fields(p)%values(1,ii,jj,kk)=exact_t
            enddo
          enddo
        enddo
      enddo
      reflux_specific_energy=global%fields(1)%values(5,1,1,1)
      oracle=global;oracle_t=reference
      if(allocated(local%fields)) deallocate(local%fields)
      if(allocated(temperature%fields)) deallocate(temperature%fields)
      allocate(local%fields(3),temperature%fields(3))
      do p=1,3
        if(owners(p)/=rank) cycle
        local%fields(p)=global%fields(p);temperature%fields(p)=reference%fields(p)
      enddo
      serial_initial=oracle;serial_temp_initial=oracle_t;mpi_initial=local;mpi_temp_initial=temperature
      call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,before,ok)
      call require(ok,'initial repair inventory')
      serial_budget=17;serial_external=18;serial_largest=19;serial_theta=20
      budget=17;external=18;largest=19;theta=20
      if(cycling) then
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,subcycling=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=sf)
      else
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=sf)
      endif
      call require(ok.eqv.(.not.failure),'serial full split repair acceptance')
      if(cycling) then
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,subcycling=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=mf)
      else
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          evaluate_diffusion=reflux_diffusion,integrated_flux=mf)
      endif
      call require(ok.eqv.(.not.failure),'MPI full split repair acceptance')
      call compare_state()
      call require(all(budget==serial_budget).and.all(external==serial_external).and. &
        largest==serial_largest.and.theta==serial_theta,'full split repair diagnostics')
      if(failure) then
        call require(all(serial_budget==17).and.all(serial_external==18).and.serial_largest==19.and.serial_theta==20, &
          'failed repair diagnostic rollback')
        do p=1,3
          call require(all(oracle%fields(p)%values==serial_initial%fields(p)%values).and. &
            all(oracle_t%fields(p)%values==serial_temp_initial%fields(p)%values),'serial split atomic rollback')
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==mpi_initial%fields(p)%values).and. &
            all(temperature%fields(p)%values==mpi_temp_initial%fields(p)%values),'MPI split atomic rollback')
        enddo
      else
        call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,after,ok)
        call require(ok.and.maxval(abs(after-before)/max(1._dp,abs(before)))<2.e-12_dp,'full split conservation')
        call require(oracle%fields(1)%values(6,3,2,2)>0,'nonzero coarse/fine synchronization exercised')
        do p=1,3
          call require(all(oracle%fields(p)%values(6:12,:,:,:)>=0),'full split species nonnegative')
          if(owners(p)/=rank) cycle
          call require(all(sf(p)%x==mf(p)%x).and.all(sf(p)%y==mf(p)%y).and.all(sf(p)%z==mf(p)%z), &
            'full split synchronized flux parity')
        enddo
      endif
    enddo
  end subroutine check_reflux_budget
  subroutine reflux_diffusion(state,temps,selected_depth,fractions,coarse_start,coarse_end,rates,flux,good, &
      internal_energy)
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end,internal_energy
    integer, intent(in) :: selected_depth
    real(dp), optional, intent(in) :: fractions(:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: good
    integer :: p,nn(3)
    rates=state
    if(allocated(flux)) deallocate(flux)
    allocate(flux(size(tree%nodes)))
    do p=1,size(tree%nodes)
      if(.not.allocated(state%fields(p)%values)) cycle
      rates%fields(p)%values=0
      if(selected_depth>=0.and.tree%nodes(p)%level/=selected_depth) cycle
      nn=tree%nodes(p)%shape
      allocate(flux(p)%x(12,nn(1)+1,nn(2),nn(3)),flux(p)%y(12,nn(1),nn(2)+1,nn(3)), &
        flux(p)%z(12,nn(1),nn(2),nn(3)+1))
      flux(p)%x=0;flux(p)%y=0;flux(p)%z=0
      ! A constant conservative fine flux has zero local divergence. Its
      ! mismatch to coarse zero flux exhausts the exposed coarse fuel only
      ! during reflux, after the stage states have already been admitted.
      if(p==1) cycle
      flux(p)%x(1,:,:,:)=1.e6_dp
      flux(p)%x(5,:,:,:)=1.e6_dp*reflux_specific_energy
      flux(p)%x(6,:,:,:)=1.e6_dp
    enddo
    good=.true.
  end subroutine reflux_diffusion
  subroutine check_auxiliary()
    real(dp) :: external(12),exchange(12),force(3),expected,serial_external(12),serial_exchange(12)
    type(amr_tree_state_3d) :: retained,retained_t
    integer :: peer,cells(3)
    owners=0;faces='periodic';dt=1e-8_dp;force=0
    do peer=1,3
      cells=tree%nodes(peer)%shape
      do k=1,cells(3)
      do j=1,cells(2)
      do i=1,cells(1)
        global%fields(peer)%values(:,i,j,k)=global%fields(1)%values(:,1,1,1)
        reference%fields(peer)%values(:,i,j,k)=reference%fields(1)%values(:,1,1,1)
      enddo
      enddo
      enddo
    enddo
    oracle=global;oracle_t=reference;serial_budget=0;serial_external=0;serial_exchange=0
    serial_largest=0;serial_theta=1
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external, &
      serial_largest,serial_theta,ok,auxiliary_source=serial_auxiliary,auxiliary_change=serial_exchange)
    call require(ok,'serial auxiliary oracle')
    allocate(local%fields(3),temperature%fields(3))
    if(rank==0) then
      local=global;temperature=reference
    endif
    budget=0;external=0;exchange=-9;largest=0;theta=1
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD, &
      ok,auxiliary_source=test_auxiliary,auxiliary_change=exchange)
    call require(ok.and.auxiliary_calls==3,'auxiliary stages 0,1,2')
    expected=3e6_dp*dt*product(spacing)*8
    call require(abs(exchange(5)-expected)<1e-14_dp*expected,'composite auxiliary inventory')
    call require(all(exchange(:4)==0).and.all(exchange(6:)==0).and.all(external==0),'separate source inventory')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'auxiliary serial MPI budget parity')
    call require(maxval(abs(exchange-serial_exchange))<1e-14_dp*expected,'auxiliary serial MPI source parity')
    if(rank==0) then
      do peer=1,3
        call require(maxval(abs(local%fields(peer)%values(5,:,:,:)- &
          global%fields(peer)%values(5,:,:,:)-3e6_dp*dt))<1e-10_dp,'accepted trapezoidal auxiliary energy')
      enddo
    endif
    retained=oracle;retained_t=oracle_t;serial_auxiliary_failure=.true.
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external, &
      serial_largest,serial_theta,ok,auxiliary_source=serial_auxiliary,auxiliary_change=serial_exchange)
    call require(.not.ok,'serial late auxiliary failure rejected')
    do peer=1,3
      call require(all(oracle%fields(peer)%values==retained%fields(peer)%values),'serial auxiliary state rollback')
      call require(all(oracle_t%fields(peer)%values==retained_t%fields(peer)%values),'serial auxiliary temperature rollback')
    enddo
    call require(all(budget==serial_budget).and.all(external==serial_external),'serial auxiliary budget rollback')
    call require(maxval(abs(exchange-serial_exchange))<1e-14_dp*expected,'serial auxiliary source rollback')
    serial_auxiliary_failure=.false.
    saved=local;saved_t=temperature;serial_budget=budget;serial_largest=largest;serial_theta=theta
    auxiliary_calls=0;auxiliary_failure=.true.
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD, &
      ok,auxiliary_source=test_auxiliary,auxiliary_change=exchange)
    call require(.not.ok,'late idle-rank auxiliary rejection')
    call unchanged()
    call require(abs(exchange(5)-expected)<1e-14_dp*expected.and.all(external==0),'auxiliary ledger rollback')
  end subroutine
  subroutine check_recursive_auxiliary()
    real(dp) :: external(12),exchange(12),force(3),expected,serial_external(12),serial_exchange(12)
    integer(int64) :: counts(3)
    type(amr_tree_state_3d) :: retained,retained_t
    integer :: peer,cells(3)
    owners=0;faces='periodic';dt=1e-8_dp;force=0
    do peer=1,3
      cells=tree%nodes(peer)%shape
      do k=1,cells(3)
      do j=1,cells(2)
      do i=1,cells(1)
        global%fields(peer)%values(:,i,j,k)=global%fields(1)%values(:,1,1,1)
        reference%fields(peer)%values(:,i,j,k)=reference%fields(1)%values(:,1,1,1)
      enddo
      enddo
      enddo
    enddo
    oracle=global;oracle_t=reference;serial_budget=0;serial_external=0;serial_exchange=0
    serial_largest=0;serial_theta=1
    call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external, &
      serial_largest,serial_theta,ok,auxiliary_source=serial_recursive_auxiliary,auxiliary_change=serial_exchange)
    call require(ok,'serial auxiliary oracle')
    allocate(local%fields(3),temperature%fields(3))
    if(rank==0) then
      local=global;temperature=reference
    endif
    budget=0;external=0;exchange=-9;largest=0;theta=1;counts=-1
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD, &
      ok,auxiliary_source=test_recursive_auxiliary,auxiliary_change=exchange,level_steps=counts)
    call require(ok.and.all(counts==[1_int64,2_int64,6_int64]),'recursive auxiliary interval counts')
    call require(all(recursive_auxiliary_calls==[3,6,18]),'recursive auxiliary stage counts')
    expected=2.5e6_dp*dt*product(spacing)*8
    call require(abs(exchange(5)-expected)<1e-14_dp*expected,'time-dependent composite auxiliary inventory')
    call require(all(exchange(:4)==0).and.all(exchange(6:)==0).and.all(external==0),'recursive source inventory split')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'auxiliary serial MPI budget parity')
    call require(maxval(abs(exchange-serial_exchange))<1e-14_dp*expected,'auxiliary serial MPI source parity')
    retained=oracle;retained_t=oracle_t;serial_auxiliary_failure=.true.
    call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external, &
      serial_largest,serial_theta,ok,auxiliary_source=serial_recursive_auxiliary,auxiliary_change=serial_exchange)
    call require(.not.ok,'serial late auxiliary failure rejected')
    do peer=1,3
      call require(all(oracle%fields(peer)%values==retained%fields(peer)%values),'serial auxiliary state rollback')
      call require(all(oracle_t%fields(peer)%values==retained_t%fields(peer)%values),'serial auxiliary temperature rollback')
    enddo
    call require(all(budget==serial_budget).and.all(external==serial_external),'serial auxiliary budget rollback')
    call require(maxval(abs(exchange-serial_exchange))<1e-14_dp*expected,'serial auxiliary source rollback')
    serial_auxiliary_failure=.false.
    saved=local;saved_t=temperature;serial_budget=budget;serial_largest=largest;serial_theta=theta
    recursive_auxiliary_calls=0;auxiliary_failure=.true.
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD, &
      ok,auxiliary_source=test_recursive_auxiliary,auxiliary_change=exchange,level_steps=counts)
    call require(.not.ok,'late fine-level idle-rank auxiliary rejection')
    call unchanged()
    call require(all(counts==[1_int64,2_int64,6_int64]),'recursive failure count rollback')
    call require(abs(exchange(5)-expected)<1e-14_dp*expected.and.all(external==0),'recursive auxiliary ledger rollback')
    auxiliary_failure=.false.
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2, &
      .false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD, &
      ok,level_steps=counts)
    call require(ok.and.all(recursive_auxiliary_calls==[3,6,18]),'absent callback after auxiliary-enabled call')
  end subroutine
  subroutine test_recursive_auxiliary(state,temps,interval,stage,selected_depth,time,rate,source_comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,selected_depth
    type(MPI_Comm), intent(in) :: source_comm
    logical, intent(out) :: good
    integer :: peer,ordinal,steps(3)
    real(dp) :: evaluation_time
    steps=[1,2,6];ordinal=recursive_auxiliary_calls(selected_depth+1)
    call require(stage==mod(ordinal,3),'recursive auxiliary stage ordering')
    call require(abs(interval-dt/steps(selected_depth+1))<1e-23_dp,'level interval duration')
    call require(abs(time-dt*(ordinal/3)/steps(selected_depth+1))<1e-23_dp,'level interval start time')
    recursive_auxiliary_calls(selected_depth+1)=ordinal+1
    evaluation_time=time
    if(stage>0) evaluation_time=time+interval
    do peer=1,3
      if(owners(peer)/=rank) cycle
      rate%fields(peer)%values=0
      if(tree%nodes(peer)%level/=selected_depth) cycle
      rate%fields(peer)%values(5,:,:,:)=2e6_dp+1e6_dp*evaluation_time/dt
    enddo
    good=.not.(auxiliary_failure.and.rank==nranks-1.and.selected_depth==2.and.stage==2.and.time>.8_dp*dt)
  end subroutine
  subroutine test_auxiliary(state,temps,interval,stage,selected_depth,rate,source_comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,selected_depth
    type(MPI_Comm), intent(in) :: source_comm
    logical, intent(out) :: good
    integer :: peer
    auxiliary_calls=auxiliary_calls+1
    call require(stage==auxiliary_calls-1.and.interval==dt.and.selected_depth==-1,'callback stage contract')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      rate%fields(peer)%values=0
      rate%fields(peer)%values(5,:,:,:)=merge(2e6_dp,4e6_dp,stage==0)
    enddo
    good=.not.(auxiliary_failure.and.rank==nranks-1.and.stage==2)
  end subroutine
  subroutine serial_auxiliary(state,temps,interval,stage,selected_depth,rate,good)
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,selected_depth
    logical, intent(out) :: good
    integer :: peer
    do peer=1,3
      rate%fields(peer)%values=0
      rate%fields(peer)%values(5,:,:,:)=merge(2e6_dp,4e6_dp,stage==0)
    enddo
    good=.not.(serial_auxiliary_failure.and.stage==2)
  end subroutine
  subroutine serial_recursive_auxiliary(state,temps,interval,stage,selected_depth,time,rate,good)
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,selected_depth
    logical, intent(out) :: good
    integer :: peer
    real(dp) :: evaluation_time
    evaluation_time=time
    if(stage>0) evaluation_time=time+interval
    do peer=1,3
      rate%fields(peer)%values=0
      if(tree%nodes(peer)%level/=selected_depth) cycle
      rate%fields(peer)%values(5,:,:,:)=2e6_dp+1e6_dp*evaluation_time/dt
    enddo
    good=.not.(serial_auxiliary_failure.and.stage==2.and.selected_depth==2.and.time>.8_dp*dt)
  end subroutine
  subroutine reactor_adapter(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    reactor_calls=reactor_calls+1
    if(reactor_failure.and.rank==0.and.reactor_calls==2) then
      species_density=-1;temperature=-1;completed_steps=0;ok=.false.
      return
    end if
    call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps=maximum_steps)
  end subroutine reactor_adapter
  subroutine check_old_admission()
    type(amr_tree_state_3d) :: u,t,h,prior,prior_t,prior_h
    type(dual_energy_policy) :: policy
    real(dp) :: ledger,bound(12),external(12),maximum,minimum
    integer :: peer,kind
    owners=mod([0,1,0],nranks)
    allocate(u%fields(3),t%fields(3),h%fields(3))
    do peer=1,3
      if(owners(peer)/=rank) cycle
      u%fields(peer)=global%fields(peer);t%fields(peer)=reference%fields(peer);h%fields(peer)=reference%fields(peer)
      h%fields(peer)%values(1,:,:,:)=u%fields(peer)%values(5,:,:,:)- &
        sum(u%fields(peer)%values(2:4,:,:,:)**2,dim=1)/(2*u%fields(peer)%values(1,:,:,:))
    enddo
    prior=u;prior_t=t;prior_h=h
    policy%eta2=2;policy%allow_negative_energy=.false.;faces='periodic'
    do kind=1,2
      u=prior;t=prior_t;h=prior_h
      if(rank==owners(2)) then
        if(kind==1) u%fields(2)%values(1,1,1,1)=-1
        if(kind==2) h%fields(2)%values(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      endif
      ledger=-7;bound=-3;external=-4;maximum=-5;minimum=-6
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,u,t,spacing,1e-8_dp,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
        bound,external,maximum,minimum,MPI_COMM_WORLD,ok,mol_mode=.true.,internal_energy=h, &
        energy_correction=ledger,energy_policy=policy)
      call require(.not.ok,'all ranks reject invalid owner old state')
      call require(ledger==-7.and.all(bound==-3).and.all(external==-4).and.maximum==-5.and.minimum==-6, &
        'MPI admission ledger rollback')
      if(rank==owners(2)) then
        if(kind==1) then
          call require(u%fields(2)%values(1,1,1,1)==-1,'MPI invalid density unchanged')
          u%fields(2)%values(1,1,1,1)=prior%fields(2)%values(1,1,1,1)
        else
          call require(ieee_is_nan(h%fields(2)%values(1,1,1,1)),'MPI invalid held unchanged')
          h%fields(2)=prior_h%fields(2)
        endif
      endif
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(u%fields(peer)%values==prior%fields(peer)%values).and. &
          all(t%fields(peer)%values==prior_t%fields(peer)%values).and. &
          all(h%fields(peer)%values==prior_h%fields(peer)%values),'MPI admission field rollback')
      enddo
    enddo
  end subroutine check_old_admission

  subroutine check_stored_energy()
    type(amr_tree_state_3d) :: held,serial_held,prior_held,prior,prior_t
    type(amr_tree_flux_3d), allocatable :: flux(:),aflux(:),sflux(:),saflux(:),prior_flux(:),prior_aflux(:)
    type(dual_energy_policy) :: policy,other_policy
    real(dp) :: ledger,serial_ledger,external(12),serial_external(12),exchange(12),serial_exchange(12)
    real(dp) :: prim(12),ct,c
    integer(int64) :: steps(3),serial_steps(3)
    integer :: peer,layout,cycle,test_case,dims(3)
    logical :: cycling,flags(4)
    procedure(stored_inert_reactor), pointer :: reactor
    dt=1e-8_dp;faces='periodic'
    prim=0;prim(1)=.3_dp;prim(5)=202650;prim(12)=1
    do layout=1,2
      owners=mod([0,1,0],nranks)
      if(layout==2) owners=0
      do cycle=-1,1
      cycling=cycle==1
      do test_case=1,4
        do peer=1,size(tree%nodes)
          call reactive_primitive_to_conserved(species,prim,global%fields(peer)%values(:,1,1,1),ct,c,ok)
          call require(ok,'stored MPI uniform N2 state')
          do i=1,12
            global%fields(peer)%values(i,:,:,:)=global%fields(peer)%values(i,1,1,1)
          enddo
          reference%fields(peer)%values=ct
        enddo
        if(test_case==1) then
          prim(5)=prim(5)*1.02_dp
          call reactive_primitive_to_conserved(species,prim,global%fields(3)%values(:,1,1,1), &
            reference%fields(3)%values(1,1,1,1),c,ok)
          call require(ok,'stored MPI finest pressure perturbation')
          prim(5)=202650
        endif
        oracle=global;oracle_t=reference;serial_held=reference
        do peer=1,size(tree%nodes)
          serial_held%fields(peer)%values(1,:,:,:)=global%fields(peer)%values(5,:,:,:)
          if(test_case==1) serial_held%fields(peer)%values=serial_held%fields(peer)%values-100
        enddo
        if(allocated(local%fields)) deallocate(local%fields,temperature%fields,held%fields)
        allocate(local%fields(3),temperature%fields(3),held%fields(3))
        do peer=1,3
          if(owners(peer)/=rank) cycle
          local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
          held%fields(peer)=serial_held%fields(peer)
        enddo
        policy=dual_energy_policy();policy%eta2=1;policy%allow_negative_energy=.false.
        policy%update_total_from_internal=test_case/=1
        ledger=0;serial_ledger=0;flags=cycling.and.test_case==2
        nullify(reactor)
        if(test_case>2) reactor=>stored_inert_reactor
        if(cycle==-1) then
          steps=1;serial_steps=1
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          merge(2,1,test_case==4),test_case>2,flags,[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          internal_energy=serial_held,energy_correction=serial_ledger,energy_policy=policy, &
          auxiliary_source=stored_serial_level_source,auxiliary_change=serial_exchange,evaluate_reaction=reactor, &
          integrated_flux=sflux,internal_energy_integrated_flux=saflux)
        call require(ok,'serial retained MOL time oracle')
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          merge(2,1,test_case==4),test_case>2,flags,[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          internal_energy=held,energy_correction=ledger,energy_policy=policy, &
          auxiliary_source=stored_mpi_level_source,auxiliary_change=exchange,evaluate_reaction=reactor, &
          integrated_flux=flux,internal_energy_integrated_flux=aflux)
        else
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          merge(2,1,test_case==4),test_case>2,flags,[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,subcycling=cycling, &
          level_steps=serial_steps,internal_energy=serial_held,energy_correction=serial_ledger,energy_policy=policy, &
          auxiliary_source=stored_serial_source,auxiliary_change=serial_exchange,evaluate_reaction=reactor, &
          integrated_flux=sflux,internal_energy_integrated_flux=saflux)
        call require(ok,'serial retained MOL time oracle')
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          merge(2,1,test_case==4),test_case>2,flags,[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,subcycling=cycling, &
          level_steps=steps,internal_energy=held,energy_correction=ledger,energy_policy=policy, &
          auxiliary_source=stored_mpi_source,auxiliary_change=exchange,evaluate_reaction=reactor, &
          integrated_flux=flux,internal_energy_integrated_flux=aflux)
        endif
        if(.not.ok.and.rank==0) print *, 'stored fixture layout/cycle/case:',layout,cycle,test_case
        call require(ok,'MPI retained MOL time update')
        if(rank==0) print *, 'retained comparison layout/cycle/case:',layout,cycle,test_case
        call compare_state()
        call require(ledger==serial_ledger.and.all(steps==serial_steps),'retained MPI correction and schedule exact')
        call require(all(budget==serial_budget).and.all(external==serial_external).and.all(exchange==serial_exchange), &
          'retained MPI boundary and source ledgers exact')
        do peer=1,3
          call require(allocated(held%fields(peer)%values).eqv.(owners(peer)==rank),'retained scalar owner-local only')
          call require(allocated(aflux(peer)%x).eqv.(owners(peer)==rank),'retained scalar flux owner-local only')
          if(owners(peer)/=rank) cycle
          call require(all(held%fields(peer)%values==serial_held%fields(peer)%values),'retained MPI scalar time state exact')
          call require(all(aflux(peer)%x==saflux(peer)%x).and.all(aflux(peer)%y==saflux(peer)%y).and. &
            all(aflux(peer)%z==saflux(peer)%z),'retained MPI synchronized scalar flux exact')
        enddo
      enddo
      enddo
      prior=local;prior_t=temperature;prior_held=held;prior_flux=flux;prior_aflux=aflux
      ledger=-7;budget=-3;external=-4;exchange=-8;largest=-5;theta=-6;steps=-9
      stored_energy_failure=.true.
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
        budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
        level_steps=steps,internal_energy=held,energy_correction=ledger,energy_policy=policy, &
        auxiliary_source=stored_mpi_source,auxiliary_change=exchange,integrated_flux=flux,internal_energy_integrated_flux=aflux)
      stored_energy_failure=.false.
      call require(.not.ok.and.ledger==-7.and.all(budget==-3).and.all(external==-4).and.all(exchange==-8).and. &
        largest==-5.and.theta==-6.and.all(steps==-9),'retained recursive collective ledger rollback')
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(local%fields(peer)%values==prior%fields(peer)%values).and. &
          all(temperature%fields(peer)%values==prior_t%fields(peer)%values).and. &
          all(held%fields(peer)%values==prior_held%fields(peer)%values),'retained recursive collective field rollback')
        call require(all(flux(peer)%x==prior_flux(peer)%x).and.all(aflux(peer)%x==prior_aflux(peer)%x), &
          'retained recursive collective flux rollback')
      enddo
      other_policy=policy
      if(rank==nranks-1) other_policy%eta2=.5_dp
      call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
        1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
        budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
        internal_energy=held,energy_correction=ledger,energy_policy=other_policy)
      call require(.not.ok.and.ledger==-7,'rank-disagreeing retained policy rejected')
      if(rank==nranks-1) then
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,internal_energy=held)
      else
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2e-7_dp,1e-12_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,internal_energy=held,energy_correction=ledger)
      endif
      call require(.not.ok.and.ledger==-7,'rank-disagreeing retained ledger presence rejected')
    enddo
  end subroutine check_stored_energy

  subroutine stored_serial_source(gas,temp,interval,stage,selected_depth,time,rate,good)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,selected_depth
    logical, intent(out) :: good
    integer :: peer
    rate=gas
    do peer=1,3
      rate%fields(peer)%values=0
      if(tree%nodes(peer)%level/=selected_depth) cycle
      rate%fields(peer)%values(5,:,:,:)=-10/dt
    enddo
    good=.true.
  end subroutine stored_serial_source

  subroutine stored_mpi_source(gas,temp,interval,stage,selected_depth,time,rate,source_comm,good)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,selected_depth
    type(MPI_Comm), intent(in) :: source_comm
    logical, intent(out) :: good
    integer :: peer
    rate=gas
    do peer=1,3
      if(owners(peer)/=rank) cycle
      rate%fields(peer)%values=0
      if(tree%nodes(peer)%level/=selected_depth) cycle
      rate%fields(peer)%values(5,:,:,:)=-10/dt
    enddo
    good=.true.
    if(stored_energy_failure.and.rank==owners(3).and.selected_depth==2.and.stage>0) good=.false.
  end subroutine stored_mpi_source

  subroutine stored_serial_level_source(gas,temp,interval,stage,selected_depth,rate,good)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,selected_depth
    logical, intent(out) :: good
    call stored_serial_source(gas,temp,interval,stage,selected_depth,0._dp,rate,good)
  end subroutine stored_serial_level_source

  subroutine stored_mpi_level_source(gas,temp,interval,stage,selected_depth,rate,source_comm,good)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,selected_depth
    type(MPI_Comm), intent(in) :: source_comm
    logical, intent(out) :: good
    call stored_mpi_source(gas,temp,interval,stage,selected_depth,0._dp,rate,source_comm,good)
  end subroutine stored_mpi_level_source

  subroutine stored_inert_reactor(sp,rx,rhoe,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: rhoe,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    densities=densities+duration*external_y;count=1;good=.true.
  end subroutine stored_inert_reactor

  subroutine check_rotating_mol(physical_boundaries)
    logical, intent(in) :: physical_boundaries
    type(rotating_frame_3d) :: frame,bad_frame
    type(amr_tree_state_3d) :: initial_state,held0,serial_held,mpi_held
    type(dual_energy_policy) :: policy
    real(dp) :: lower(3),initial(12),final(12),external(12),serial_external(12),residual,case_dt
    real(dp) :: saved_budget(12),saved_external(12)
    real(dp) :: serial_energy_change,mpi_energy_change
    logical :: recursive,reacting,case_transport(4)
    integer :: peer,variant,first_variant,last_variant
    frame%axis=3;frame%angular_speed=400._dp
    frame%axis_point=[.0004_dp,-.0007_dp,0._dp]
    lower=[-.002_dp,.003_dp,.001_dp]
    faces=['outflow ','outflow ','outflow ','outflow ','periodic','periodic']
    owners=mod([0,1,0],nranks)
    initial_state=global
    call convert_amr_patch_tree_rotating_energy_3d(tree,initial_state,lower,spacing,frame,.true.,ok)
    call require(ok,'rotating regular AMR initial state')
    call composite_integrals_amr_patch_tree_3d(tree,initial_state,spacing,initial,ok)
    call require(ok,'rotating regular AMR initial inventory')
    case_transport=.false.;case_dt=1e-8_dp
    first_variant=1;last_variant=12
    if(physical_boundaries) then
      first_variant=13;last_variant=20
    endif
    do variant=first_variant,last_variant
      recursive=mod(variant,2)==0
      reacting=(variant>=7.and.variant<=8).or.(variant>=11.and.variant<=12)
      case_transport=variant>=9
      case_dt=merge(1e-10_dp,1e-8_dp,variant>=7)
      faces=['outflow ','outflow ','outflow ','outflow ','periodic','periodic']
      if(variant>=3) faces(1)='inflow'
      mach=-1;pin=0;pressure=0;wall_t=900;wall_v=0
      if(variant>=5) mach(1)=.2_dp
      if(variant>=13) then
        faces(1)='outflow';mach=-1
        select case((variant-13)/2)
        case(0)
          faces(1:2)='slip_wall'
        case(1)
          faces(1:2)='no_slip_wall';wall_v(2,1:2)=3._dp
        case(2)
          faces(1:2)='isothermal_wall';wall_t(1:2)=850._dp
        case(3)
          pressure(2)=2.08e5_dp
        end select
      endif
      oracle=initial_state;oracle_t=reference
      if(allocated(local%fields)) deallocate(local%fields,temperature%fields)
      allocate(local%fields(3),temperature%fields(3))
      do peer=1,3
        if(owners(peer)/=rank) cycle
        local%fields(peer)=initial_state%fields(peer)
        temperature%fields(peer)=reference%fields(peer)
      enddo
      if(recursive) then
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          2,reacting,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          outlet_pressure=pressure,wall_temperature=wall_t,wall_velocity=wall_v)
        call require(ok,'serial recursive rotating regular AMR step')
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt, &
          faces,2,reacting,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          outlet_pressure=pressure,wall_temperature=wall_t,wall_velocity=wall_v)
      else
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          2,reacting,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          outlet_pressure=pressure,wall_temperature=wall_t,wall_velocity=wall_v)
        call require(ok,'serial common rotating regular AMR step')
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          2,reacting,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          outlet_pressure=pressure,wall_temperature=wall_t,wall_velocity=wall_v)
      endif
      call require(ok,'MPI rotating regular AMR step')
      call compare_state()
      call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,final,ok)
      call require(ok,'rotating regular AMR final inventory')
      residual=maxval(abs(final-initial-serial_budget-serial_external)/max(1.e-18_dp,abs(initial)))
      if(reacting) residual=maxval(abs(final(:5)-initial(:5)-serial_budget(:5)-serial_external(:5))/ &
        max(1.e-18_dp,abs(initial(:5))))
      if(rank==0) print *, 'rotating regular MOL variant, composite residual:',variant,residual
      call require(residual<3e-12_dp,'rotating regular AMR composite conservation')
      saved=local;saved_t=temperature;saved_budget=budget;saved_external=external
      bad_frame=frame
      if(rank==nranks-1) bad_frame%angular_speed=bad_frame%angular_speed+1._dp
      if(recursive) then
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt, &
          faces,2,reacting,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,frame=bad_frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          outlet_pressure=pressure,wall_temperature=wall_t,wall_velocity=wall_v)
      else
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          2,reacting,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,frame=bad_frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          outlet_pressure=pressure,wall_temperature=wall_t,wall_velocity=wall_v)
      endif
      call require(.not.ok,'rank-disagreeing rotating frame rejected')
      call require(all(budget==saved_budget).and.all(external==saved_external),'rotating rollback ledgers')
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(local%fields(peer)%values==saved%fields(peer)%values),'rotating rollback state')
        call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'rotating rollback temperature')
      enddo
    enddo
    if(physical_boundaries) return
    case_transport=.false.;case_dt=1e-8_dp
    held0=reference
    do peer=1,3
      held0%fields(peer)%values(1,:,:,:)=global%fields(peer)%values(5,:,:,:)- &
        .5_dp*sum(global%fields(peer)%values(2:4,:,:,:)**2,dim=1)/global%fields(peer)%values(1,:,:,:)
    enddo
    policy=dual_energy_policy();policy%eta2=1;policy%allow_negative_energy=.false.
    do variant=1,4
      case_transport=variant>=3
      case_dt=merge(1e-10_dp,1e-8_dp,variant>=3)
      faces=['outflow ','outflow ','outflow ','outflow ','periodic','periodic']
      if(mod(variant,2)==0) faces(1)='inflow'
      mach=-1;pin=0
      if(mod(variant,2)==0) mach(1)=.2_dp
      oracle=initial_state;oracle_t=reference;serial_held=held0
      if(allocated(local%fields)) deallocate(local%fields,temperature%fields)
      if(allocated(mpi_held%fields)) deallocate(mpi_held%fields)
      allocate(local%fields(3),temperature%fields(3),mpi_held%fields(3))
      do peer=1,3
        if(owners(peer)/=rank) cycle
        local%fields(peer)=initial_state%fields(peer)
        temperature%fields(peer)=reference%fields(peer)
        mpi_held%fields(peer)=held0%fields(peer)
      enddo
      serial_energy_change=0;mpi_energy_change=0
      if(mod(variant,2)==1) then
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          2,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          internal_energy=serial_held,energy_correction=serial_energy_change,energy_policy=policy)
        call require(ok,'serial rotating retained common step')
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          2,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          internal_energy=mpi_held,energy_correction=mpi_energy_change,energy_policy=policy)
      else
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          2,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          internal_energy=serial_held,energy_correction=serial_energy_change,energy_policy=policy)
        call require(ok,'serial rotating retained recursive step')
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing, &
          case_dt,faces,2,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,frame=frame,domain_lower=lower, &
          inflow_state=inlet,inflow_temperature=inlet_t,inflow_mach=mach,inflow_pressure=pin, &
          internal_energy=mpi_held,energy_correction=mpi_energy_change,energy_policy=policy)
      endif
      call require(ok,'MPI rotating retained step')
      call compare_state()
      call require(mpi_energy_change==serial_energy_change,'rotating retained energy correction parity')
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(mpi_held%fields(peer)%values==serial_held%fields(peer)%values), &
          'rotating retained internal energy parity')
      enddo
    enddo
  end subroutine check_rotating_mol

  subroutine check_mol_conservation()
    type(amr_tree_flux_3d), allocatable :: registered(:),physical(:),mpi_registered(:),mpi_physical(:)
    type(amr_tree_flux_3d), allocatable :: scalar_registered(:),scalar_physical(:), &
      mpi_scalar_registered(:),mpi_scalar_physical(:)
    type(amr_tree_state_3d) :: serial_held,mpi_held,initial_held
    type(dual_energy_policy) :: policy
    real(dp) :: initial(12),final(12),expected(12),residual,external(12),serial_external(12),case_dt
    real(dp) :: initial_scalar(1),final_scalar(1),serial_energy_change,mpi_energy_change
    integer :: peer,variant,mol_iterations,ncell(3)
    logical :: recursive,case_transport(4)
    owners=mod([0,1,0],nranks);faces='periodic';dt=1e-8_dp
    call composite_integrals_amr_patch_tree_3d(tree,global,spacing,initial,ok)
    call require(ok,'initial regular MOL composite inventory')
    mol_iterations=2
    do variant=1,3
      recursive=variant>=2
      case_transport=variant==3
      case_dt=merge(1e-10_dp,dt,variant==3)
      oracle=global;oracle_t=reference
      if(allocated(local%fields)) deallocate(local%fields,temperature%fields)
      allocate(local%fields(3),temperature%fields(3))
      do peer=1,3
        if(owners(peer)/=rank) cycle
        local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
      enddo
      if(recursive) then
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          integrated_flux=registered,physical_integrated_flux=physical)
        call require(ok,'serial recursive nonreacting MOL conservation step')
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          integrated_flux=mpi_registered,physical_integrated_flux=mpi_physical)
      else
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          integrated_flux=registered,physical_integrated_flux=physical)
        call require(ok,'serial common nonreacting MOL conservation step')
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          integrated_flux=mpi_registered,physical_integrated_flux=mpi_physical)
      endif
      call require(ok,'MPI nonreacting MOL conservation step')
      call compare_state()
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(mpi_registered(peer)%x==registered(peer)%x).and. &
          all(mpi_registered(peer)%y==registered(peer)%y).and. &
          all(mpi_registered(peer)%z==registered(peer)%z),'regular MOL registered serial MPI parity')
        call require(all(mpi_physical(peer)%x==physical(peer)%x).and. &
          all(mpi_physical(peer)%y==physical(peer)%y).and. &
          all(mpi_physical(peer)%z==physical(peer)%z),'regular MOL physical serial MPI parity')
      enddo
      if(rank==owners(1)) call require(any(physical(1)%x/=registered(1)%x).or. &
        any(physical(1)%y/=registered(1)%y).or.any(physical(1)%z/=registered(1)%z), &
        'regular MOL physical and registered transfers differ')
      call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,final,ok)
      call require(ok,'final regular MOL composite inventory')
      residual=maxval(abs(final-initial-serial_budget-serial_external)/max(1._dp,abs(initial)))
      if(rank==0) print *, 'regular MOL variant, composite residual:',variant,residual
      call require(residual<3e-14_dp,'nonreacting MOL composite conservation')
    enddo
    policy=dual_energy_policy();policy%allow_negative_energy=.false.
    policy%eta2=1;policy%update_total_from_internal=.false.
    allocate(serial_held%fields(3),mpi_held%fields(3))
    do peer=1,3
      ncell=tree%nodes(peer)%shape
      allocate(serial_held%fields(peer)%values(1,ncell(1),ncell(2),ncell(3)))
      serial_held%fields(peer)%values(1,:,:,:)=global%fields(peer)%values(5,:,:,:)- &
        .5_dp*sum(global%fields(peer)%values(2:4,:,:,:)**2,dim=1)/global%fields(peer)%values(1,:,:,:)
    enddo
    call composite_integrals_amr_patch_tree_3d(tree,serial_held,spacing,initial_scalar,ok)
    call require(ok,'initial retained MOL composite inventory')
    initial_held=serial_held
    do variant=1,3
      recursive=variant>=2
      case_transport=variant==3
      case_dt=merge(1e-10_dp,dt,case_transport(1))
      oracle=global;oracle_t=reference;serial_held=initial_held
      serial_energy_change=0;mpi_energy_change=0
      if(allocated(local%fields)) deallocate(local%fields,temperature%fields)
      allocate(local%fields(3),temperature%fields(3))
      do peer=1,3
        if(owners(peer)/=rank) cycle
        local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
        mpi_held%fields(peer)=serial_held%fields(peer)
      enddo
      if(recursive) then
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          internal_energy=serial_held,energy_correction=serial_energy_change,energy_policy=policy, &
          internal_energy_integrated_flux=scalar_registered, &
          physical_internal_energy_integrated_flux=scalar_physical)
        call require(ok,'serial recursive retained MOL conservation step')
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          internal_energy=mpi_held,energy_correction=mpi_energy_change,energy_policy=policy, &
          internal_energy_integrated_flux=mpi_scalar_registered, &
          physical_internal_energy_integrated_flux=mpi_scalar_physical)
      else
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
          internal_energy=serial_held,energy_correction=serial_energy_change,energy_policy=policy, &
          internal_energy_integrated_flux=scalar_registered, &
          physical_internal_energy_integrated_flux=scalar_physical)
        call require(ok,'serial common retained MOL conservation step')
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,case_dt,faces, &
          mol_iterations,.false.,case_transport,[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
          internal_energy=mpi_held,energy_correction=mpi_energy_change,energy_policy=policy, &
          internal_energy_integrated_flux=mpi_scalar_registered, &
          physical_internal_energy_integrated_flux=mpi_scalar_physical)
      endif
      call require(ok,'MPI retained MOL conservation step')
      call compare_state()
      call require(serial_energy_change==mpi_energy_change,'retained MOL serial MPI energy correction parity')
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(mpi_held%fields(peer)%values==serial_held%fields(peer)%values), &
          'retained MOL serial MPI state parity')
        call require(all(mpi_scalar_registered(peer)%x==scalar_registered(peer)%x).and. &
          all(mpi_scalar_registered(peer)%y==scalar_registered(peer)%y).and. &
          all(mpi_scalar_registered(peer)%z==scalar_registered(peer)%z), &
          'retained MOL serial MPI registered flux parity')
        call require(all(mpi_scalar_physical(peer)%x==scalar_physical(peer)%x).and. &
          all(mpi_scalar_physical(peer)%y==scalar_physical(peer)%y).and. &
          all(mpi_scalar_physical(peer)%z==scalar_physical(peer)%z), &
          'retained MOL serial MPI physical flux parity')
      enddo
      if(rank==owners(1)) call require(any(scalar_physical(1)%x/=scalar_registered(1)%x).or. &
        any(scalar_physical(1)%y/=scalar_registered(1)%y).or. &
        any(scalar_physical(1)%z/=scalar_registered(1)%z), &
        'retained MOL physical and registered transfers differ')
      call require(any(serial_held%fields(1)%values/=initial_held%fields(1)%values), &
        'retained MOL state changes during the conservation step')
      call composite_integrals_amr_patch_tree_3d(tree,serial_held,spacing,final_scalar,ok)
      call require(ok,'final retained MOL composite inventory')
      residual=abs(final_scalar(1)-initial_scalar(1))/max(1._dp,abs(initial_scalar(1)))
      if(rank==0) print *, 'retained MOL variant, composite residual:',variant,residual
      call require(residual<3e-14_dp,'retained MOL composite conservation')
      call composite_integrals_amr_patch_tree_3d(tree,oracle,spacing,final,ok)
      call require(ok,'final retained MOL gas inventory')
      expected=serial_budget+serial_external;expected(5)=expected(5)+serial_energy_change
      residual=maxval(abs(final-initial-expected)/max(1._dp,abs(initial)))
      call require(residual<3e-14_dp,'retained MOL gas composite conservation')
    enddo
  end subroutine check_mol_conservation

  subroutine check_mol()
    real(dp) :: external(12),serial_external(12),force(3)
    integer :: variant,lf_case,peer
    logical :: chemistry,recursive,lf
    owners=mod([0,1,0],nranks);dt=1e-10_dp;faces='periodic';force=[20000._dp,-30000._dp,40000._dp]
    do variant=1,4
      chemistry=mod(variant,2)==0;recursive=variant>2
      do lf_case=0,1
        lf=lf_case==1
        if(allocated(local%fields)) deallocate(local%fields,temperature%fields)
        allocate(local%fields(3),temperature%fields(3))
        do peer=1,3
          if(owners(peer)/=rank) cycle
          local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
        enddo
        oracle=global;oracle_t=reference
        if(recursive) then
          call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
            2,chemistry,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp, &
            serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,mol_lax_friedrichs=lf)
          call require(ok,'serial recursive MOL stage')
          call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
            2,chemistry,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp, &
            budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,mol_lax_friedrichs=lf)
        else
          call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
            2,chemistry,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp, &
            serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true.,mol_lax_friedrichs=lf)
          call require(ok,'serial common MOL stage')
          call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
            2,chemistry,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp, &
            budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,mol_lax_friedrichs=lf)
        endif
        call require(ok,'MPI coupled MOL update')
        call compare_state()
        call require(all(budget==serial_budget).and.all(external==serial_external),'MOL boundary/source inventories exact')
        if(largest/=serial_largest.or.theta/=serial_theta) &
          print *, 'MOL diagnostic mismatch:',rank,variant,lf_case,largest,serial_largest,theta,serial_theta
        call require(largest==serial_largest.and.theta==serial_theta,'MOL diagnostics exact')
        saved=local;saved_t=temperature
        lf=rank==nranks-1
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          2,chemistry,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
          mol_mode=.true.,mol_lax_friedrichs=lf)
        call require(.not.ok,'MOL idle-rank LF mismatch rejected')
        call unchanged()
      enddo
    enddo
  end subroutine check_mol
  subroutine check_ordered_mol_auxiliary()
    real(dp) :: external(12),serial_external(12),exchange(12),serial_exchange(12)
    real(dp) :: retained_budget(12),retained_external(12),retained_exchange(12)
    integer :: peer
    owners=mod([0,1,0],nranks);dt=1e-10_dp;faces='periodic'
    if(allocated(local%fields)) deallocate(local%fields,temperature%fields)
    allocate(local%fields(3),temperature%fields(3))
    do peer=1,3
      if(owners(peer)/=rank) cycle
      local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
    enddo
    oracle=global;oracle_t=reference
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
      1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
      serial_budget,serial_external,serial_largest,serial_theta,ok,mol_mode=.true., &
      auxiliary_source=serial_auxiliary,auxiliary_change=serial_exchange)
    call require(ok,'serial ordered MOL auxiliary adapter')
    ordered_mol_calls=0;auxiliary_failure=.false.
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
      budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
      auxiliary_source=ordered_mol_auxiliary,auxiliary_change=exchange)
    call require(ok.and.all(ordered_mol_calls==2),'common MOL two stages per level with equal dt')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'ordered MOL source budgets exact')
    call require(all(exchange==serial_exchange),'ordered MOL auxiliary inventory exact')
    saved=local;saved_t=temperature
    retained_budget=budget;retained_external=external;retained_exchange=exchange
    ordered_mol_calls=0;auxiliary_failure=.true.
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
      budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
      auxiliary_source=ordered_mol_auxiliary,auxiliary_change=exchange)
    call require(.not.ok,'ordered MOL late fine-level idle-rank callback rejection')
    call unchanged()
    call require(all(budget==retained_budget).and.all(external==retained_external), &
      'ordered MOL rejected step inventory rollback')
    call require(all(exchange==retained_exchange),'ordered MOL rejected auxiliary inventory rollback')
    auxiliary_failure=.false.
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      1,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
      budget,external,largest,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,subcycling=rank==nranks-1)
    call require(.not.ok,'rank-disagreeing level time schedule rejected')
    call unchanged()
  end subroutine check_ordered_mol_auxiliary
  subroutine ordered_mol_auxiliary(state,temps,interval,stage,selected_depth,rate,source_comm,good)
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,selected_depth
    type(MPI_Comm), intent(in) :: source_comm
    logical, intent(out) :: good
    integer :: peer
    call require(selected_depth>=0.and.selected_depth<3,'ordered MOL callback level')
    call require(stage==ordered_mol_calls(selected_depth+1).and.interval==dt,'ordered MOL callback stage and dt')
    if(selected_depth>0) call require(all(ordered_mol_calls(:selected_depth)==2),'coarse MOL completed before fine')
    ordered_mol_calls(selected_depth+1)=ordered_mol_calls(selected_depth+1)+1
    do peer=1,3
      if(owners(peer)/=rank) cycle
      rate%fields(peer)%values=0
      rate%fields(peer)%values(5,:,:,:)=merge(2e6_dp,4e6_dp,stage==0)
    enddo
    good=.not.(auxiliary_failure.and.rank==nranks-1.and.selected_depth==2.and.stage==1)
  end subroutine ordered_mol_auxiliary
  subroutine check_sdc()
    type(amr_tree_state_3d) :: ir,serial_ir,retained_ir,heat,serial_heat,retained_heat
    real(dp) :: external(12),serial_external(12),force(3),visc
    type(amr_tree_state_3d) :: start,finish,owned_start,owned_end,src0,src1,owned_src0,owned_src1,op,np,sop,snp
    type(amr_tree_flux_3d), allocatable :: flux(:),serial_flux(:)
    integer(int64) :: counts(3),serial_counts(3)
    integer :: peer,dims(3),variant
    logical :: flags(4)
    owners=mod([0,1,0],nranks);dt=1e-10_dp;faces='periodic';force=[20000._dp,-30000._dp,40000._dp]
    allocate(local%fields(3),temperature%fields(3),ir%fields(3),serial_ir%fields(3))
    do peer=1,3
      dims=tree%nodes(peer)%shape
      allocate(serial_ir%fields(peer)%values(8,dims(1),dims(2),dims(3)))
      serial_ir%fields(peer)%values=0
      if(owners(peer)/=rank) cycle
      local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
      ir%fields(peer)=serial_ir%fields(peer)
    end do
    heat=temperature;serial_heat=reference
    do peer=1,3
      serial_heat%fields(peer)%values=0
      if(owners(peer)==rank) heat%fields(peer)%values=0
    enddo
    oracle=global;oracle_t=reference
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
      ok,reaction_source=serial_ir,evaluate_reaction=reactor_adapter,heat_release=serial_heat)
    call require(ok.and.reactor_calls>0,'serial coupled SDC callback')
    reactor_calls=0
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok,reaction_source=ir, &
      evaluate_reaction=reactor_adapter,heat_release=heat)
    call require(ok,'MPI coupled SDC iteration')
    call require((reactor_calls>0).eqv.any(owners==rank),'callback only on owning ranks')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'exact MPI SDC boundary/source inventory')
    call require(largest==serial_largest.and.theta==serial_theta,'exact MPI SDC diagnostics')
    do peer=1,3
      call require(allocated(ir%fields(peer)%values).eqv.(owners(peer)==rank),'owner-local SDC estimate')
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'exact MPI SDC reaction estimates')
      call require(all(heat%fields(peer)%values==serial_heat%fields(peer)%values),'regular heat exact')
    end do
    retained_heat=heat
    saved=local;saved_t=temperature;retained_ir=ir;visc=.1_dp
    if(rank==nranks-1) visc=-1
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,difmag=visc,heat_release=heat)
    call require(.not.ok,'idle-rank invalid hydro after transport rejects whole SDC step')
    call unchanged()
    call require(all(external==serial_external),'source inventory rollback')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'reaction estimate rollback')
      call require(all(heat%fields(peer)%values==retained_heat%fields(peer)%values),'regular heat rollback')
    end do
    allocate(start%fields(3),owned_start%fields(3),owned_end%fields(3),owned_src0%fields(3),owned_src1%fields(3))
    src0=global;src1=global
    do peer=1,3
      dims=tree%nodes(peer)%shape
      allocate(start%fields(peer)%values(13,dims(1),dims(2),dims(3)))
      start%fields(peer)%values(:12,:,:,:)=oracle%fields(peer)%values
      start%fields(peer)%values(13,:,:,:)=oracle_t%fields(peer)%values(1,:,:,:)
      src0%fields(peer)%values=0;src1%fields(peer)%values=0
      src1%fields(peer)%values(2,:,:,:)=10000
    end do
    finish=start
    do peer=1,2
      finish%fields(peer)%values(:12,:,:,:)=1.001_dp*start%fields(peer)%values(:12,:,:,:)
      if(owners(peer)/=rank) cycle
      owned_start%fields(peer)=start%fields(peer);owned_end%fields(peer)=finish%fields(peer)
      owned_src0%fields(peer)=src0%fields(peer);owned_src1%fields(peer)=src1%fields(peer)
    end do
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
      ok,reaction_source=serial_ir,target_depth=2,coarse_start=start,coarse_end=finish, &
      fractions_start=[.25_dp,.5_dp],fractions_end=[.5_dp,1._dp],source_start=src0,source_end=src1, &
      integrated_flux=serial_flux,old_predictor=sop,new_predictor=snp,heat_release=serial_heat)
    call require(ok,'serial coupled level SDC')
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,target_depth=2,coarse_start=owned_start,coarse_end=owned_end, &
      fractions_start=[.25_dp,.5_dp],fractions_end=[.5_dp,1._dp],source_start=owned_src0,source_end=owned_src1, &
      integrated_flux=flux,old_predictor=op,new_predictor=np,heat_release=heat)
    call require(ok,'MPI coupled level SDC')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'level boundary/source exact')
    call require(largest==serial_largest.and.theta==serial_theta,'level diagnostics exact')
    do peer=1,3
      call require(allocated(flux(peer)%x).eqv.(peer==3.and.owners(peer)==rank),'level flux ownership')
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'level IR exact')
      call require(all(heat%fields(peer)%values==serial_heat%fields(peer)%values),'regular heat exact')
      call require(all(op%fields(peer)%values==sop%fields(peer)%values),'old predictor exact')
      call require(all(np%fields(peer)%values==snp%fields(peer)%values),'new predictor exact')
      if(peer<3) then
        call require(all(heat%fields(peer)%values==retained_heat%fields(peer)%values),'inactive regular heat')
        call require(all(local%fields(peer)%values==saved%fields(peer)%values),'inactive state unchanged')
        call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'inactive T unchanged')
      else
        call require(all(flux(peer)%x==serial_flux(peer)%x).and.all(flux(peer)%y==serial_flux(peer)%y).and. &
          all(flux(peer)%z==serial_flux(peer)%z),'level flux exact')
      end if
    end do
    call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
      ok,reaction_source=serial_ir,level_steps=serial_counts,heat_release=serial_heat)
    call require(ok,'serial recursive SDC')
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      2,.true.,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,level_steps=counts,evaluate_reaction=reactor_adapter,heat_release=heat)
    call require(ok.and.all(counts==serial_counts),'MPI recursive SDC and counts')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'recursive inventories exact')
    call require(largest==serial_largest.and.theta==serial_theta,'recursive diagnostics exact')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'recursive reaction estimates exact')
      call require(all(heat%fields(peer)%values==serial_heat%fields(peer)%values),'regular heat exact')
    end do
    retained_heat=heat
    saved=local;saved_t=temperature;retained_ir=ir;visc=.1_dp
    if(rank==nranks-1) visc=-1
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      2,.true.,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,level_steps=counts,difmag=visc,heat_release=heat)
    call require(.not.ok.and.all(counts==serial_counts),'recursive failure retains counts')
    call unchanged()
    call require(all(external==serial_external),'recursive source ledger rollback')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'recursive IR rollback')
      call require(all(heat%fields(peer)%values==retained_heat%fields(peer)%values),'regular heat rollback')
    end do
    ! Fail after successful cells, and reject callback presence disagreement on an idle rank.
    do variant=1,3
      reactor_calls=0;reactor_failure=variant<=2
      if(variant==1) then
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          2,.true.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta, &
          MPI_COMM_WORLD,ok,reaction_source=ir,evaluate_reaction=reactor_adapter)
      else if(variant==2.or.rank/=nranks-1) then
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing, &
          dt,faces,2,.true.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta, &
          MPI_COMM_WORLD,ok,reaction_source=ir,level_steps=counts,evaluate_reaction=reactor_adapter)
      else
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing, &
          dt,faces,2,.true.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta, &
          MPI_COMM_WORLD,ok,reaction_source=ir,level_steps=counts)
      end if
      call require(.not.ok,'reactor failure or rank mismatch rejected')
      call unchanged()
      call require(all(counts==serial_counts).and.all(external==serial_external),'reactor rollback ledgers')
      do peer=1,3
        if(owners(peer)/=rank) cycle
        call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'reactor history rollback')
      end do
      if(variant==3) call require(reactor_calls==0,'presence mismatch rejected before reaction')
    end do
    reactor_failure=.false.
    do variant=1,2
      flags=variant==1
      oracle=global;oracle_t=reference
      do peer=1,3
        if(owners(peer)/=rank) cycle
        local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
      end do
      if(variant==1) then
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.false., &
          flags,force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
          ok,evaluate_diffusion=test_diffusion)
        call require(ok,'serial diffusion callback')
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
          2,.false.,flags,force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
          evaluate_diffusion=test_diffusion)
      else
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
          2,.false.,flags,force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
          ok,level_steps=serial_counts,evaluate_diffusion=test_diffusion)
        call require(ok,'serial recursive diffusion callback')
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing, &
          dt,faces,2,.false.,flags,force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
          level_steps=counts,evaluate_diffusion=test_diffusion)
        call require(all(counts==serial_counts),'diffusion callback recursive counts')
      end if
      call require(ok,'MPI diffusion callback')
      call compare_state()
      call require(all(budget==serial_budget).and.all(external==serial_external),'diffusion callback inventories exact')
      saved=local;saved_t=temperature;diffusion_failure=.true.;diffusion_calls=0
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing, &
        dt,faces,2,.false.,flags,force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
        evaluate_diffusion=test_diffusion)
      call require(.not.ok,'single-rank callback failure rejects collective step')
      call unchanged()
      call require(all(budget==serial_budget).and.all(external==serial_external),'callback failure ledger rollback')
      diffusion_failure=.false.
    end do
  end subroutine check_sdc
  subroutine test_diffusion(state,temps,selected_depth,fractions,coarse_start,coarse_end,rates,flux,good, &
      internal_energy)
      type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d), intent(in) :: state,temps
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    integer, intent(in) :: selected_depth
    real(dp), optional, intent(in) :: fractions(:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: good
    integer :: peer,cells(3),face,ancestor
    real(dp) :: dx
    rates=state
    if(allocated(flux)) deallocate(flux)
    allocate(flux(size(tree%nodes)))
    do peer=1,size(tree%nodes)
      if(.not.allocated(state%fields(peer)%values)) cycle
      rates%fields(peer)%values=0
      if(selected_depth>=0.and.tree%nodes(peer)%level/=selected_depth) cycle
      cells=tree%nodes(peer)%shape;dx=spacing(1);ancestor=peer
      do while(ancestor>1)
        dx=dx/real(tree%nodes(ancestor)%patch%refinement_ratio,dp)
        ancestor=tree%nodes(ancestor)%parent
      end do
      allocate(flux(peer)%x(12,cells(1)+1,cells(2),cells(3)), &
        flux(peer)%y(12,cells(1),cells(2)+1,cells(3)),flux(peer)%z(12,cells(1),cells(2),cells(3)+1))
      flux(peer)%x=0;flux(peer)%y=0;flux(peer)%z=0
      do face=2,cells(1)
        flux(peer)%x(2,face,:,:)=1000._dp*real((face-1)*(cells(1)+1-face),dp)
      end do
      do face=1,cells(1)
        rates%fields(peer)%values(:,face,:,:)=-(flux(peer)%x(:,face+1,:,:)-flux(peer)%x(:,face,:,:))/dx
      end do
    end do
    diffusion_calls=diffusion_calls+1
    good=.not.(diffusion_failure.and.rank==0.and.diffusion_calls==2)
  end subroutine test_diffusion
  subroutine compare_state()
    integer :: peer
    do peer=1,3
      call require(allocated(local%fields(peer)%values).eqv.(owners(peer)==rank),'owned state only')
      call require(allocated(temperature%fields(peer)%values).eqv.(owners(peer)==rank),'owned temperature only')
      if (owners(peer)/=rank) cycle
      if(any(local%fields(peer)%values/=oracle%fields(peer)%values)) &
        print *, 'state mismatch:',rank,peer,maxloc(abs(local%fields(peer)%values-oracle%fields(peer)%values)), &
          maxval(abs(local%fields(peer)%values-oracle%fields(peer)%values))
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
    call require(all(budget==serial_budget).and.largest==serial_largest.and.theta==serial_theta,'diagnostic rollback')
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
