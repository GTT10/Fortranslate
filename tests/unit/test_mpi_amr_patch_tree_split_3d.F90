program test_mpi_amr_patch_tree_split_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d
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
  real(dp) :: coefficient
  logical :: diffusion_failure=.false.
  integer :: diffusion_calls=0,reactor_calls=0
  logical :: reactor_failure=.false.
  logical :: ppm
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
  subroutine check_sdc()
    type(amr_tree_state_3d) :: ir,serial_ir,retained_ir
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
    oracle=global;oracle_t=reference
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
      ok,reaction_source=serial_ir,evaluate_reaction=reactor_adapter)
    call require(ok.and.reactor_calls>0,'serial coupled SDC callback')
    reactor_calls=0
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok,reaction_source=ir, &
      evaluate_reaction=reactor_adapter)
    call require(ok,'MPI coupled SDC iteration')
    call require((reactor_calls>0).eqv.any(owners==rank),'callback only on owning ranks')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'exact MPI SDC boundary/source inventory')
    call require(largest==serial_largest.and.theta==serial_theta,'exact MPI SDC diagnostics')
    do peer=1,3
      call require(allocated(ir%fields(peer)%values).eqv.(owners(peer)==rank),'owner-local SDC estimate')
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'exact MPI SDC reaction estimates')
    end do
    saved=local;saved_t=temperature;retained_ir=ir;visc=.1_dp
    if(rank==nranks-1) visc=-1
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,difmag=visc)
    call require(.not.ok,'idle-rank invalid hydro after transport rejects whole SDC step')
    call unchanged()
    call require(all(external==serial_external),'source inventory rollback')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'reaction estimate rollback')
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
      integrated_flux=serial_flux,old_predictor=sop,new_predictor=snp)
    call require(ok,'serial coupled level SDC')
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,target_depth=2,coarse_start=owned_start,coarse_end=owned_end, &
      fractions_start=[.25_dp,.5_dp],fractions_end=[.5_dp,1._dp],source_start=owned_src0,source_end=owned_src1, &
      integrated_flux=flux,old_predictor=op,new_predictor=np)
    call require(ok,'MPI coupled level SDC')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'level boundary/source exact')
    call require(largest==serial_largest.and.theta==serial_theta,'level diagnostics exact')
    do peer=1,3
      call require(allocated(flux(peer)%x).eqv.(peer==3.and.owners(peer)==rank),'level flux ownership')
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'level IR exact')
      call require(all(op%fields(peer)%values==sop%fields(peer)%values),'old predictor exact')
      call require(all(np%fields(peer)%values==snp%fields(peer)%values),'new predictor exact')
      if(peer<3) then
        call require(all(local%fields(peer)%values==saved%fields(peer)%values),'inactive state unchanged')
        call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'inactive T unchanged')
      else
        call require(all(flux(peer)%x==serial_flux(peer)%x).and.all(flux(peer)%y==serial_flux(peer)%y).and. &
          all(flux(peer)%z==serial_flux(peer)%z),'level flux exact')
      end if
    end do
    call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
      ok,reaction_source=serial_ir,level_steps=serial_counts)
    call require(ok,'serial recursive SDC')
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      2,.true.,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,level_steps=counts,evaluate_reaction=reactor_adapter)
    call require(ok.and.all(counts==serial_counts),'MPI recursive SDC and counts')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'recursive inventories exact')
    call require(largest==serial_largest.and.theta==serial_theta,'recursive diagnostics exact')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'recursive reaction estimates exact')
    end do
    saved=local;saved_t=temperature;retained_ir=ir;visc=.1_dp
    if(rank==nranks-1) visc=-1
    call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      2,.true.,[.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,level_steps=counts,difmag=visc)
    call require(.not.ok.and.all(counts==serial_counts),'recursive failure retains counts')
    call unchanged()
    call require(all(external==serial_external),'recursive source ledger rollback')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'recursive IR rollback')
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
  subroutine test_diffusion(state,temps,selected_depth,fractions,coarse_start,coarse_end,rates,flux,good)
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
