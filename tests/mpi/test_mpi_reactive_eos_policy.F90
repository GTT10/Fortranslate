program test_mpi_reactive_eos_policy
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, configure_nasa7_eos_policy
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: gas_transport_species, load_h2o2_full_transport
  use simulation_config_reactive_1d_mod, only: reactive_1d_config
  use reactive_1d_mod, only: initialize_reactive_1d
  use mpi_domain_1d_mod, only: mpi_domain_1d, initialize_mpi_domain_1d
  use mpi_reactive_eos_policy_mod, only: agree_mpi_reactive_eos_1d
  use mpi_reactive_transport_1d_mod
  use mpi_reactive_1d_mod
  use amr_patch_tree_1d_mod, only: amr_patch_level_plan_1d
  use amr_patch_tree_reactive_1d_mod, only: amr_patch_tree_reactive_solution_1d, initialize_patch_tree_reactive_1d
  use mpi_amr_patch_1d_mod
  use mpi_amr_sparse_patch_1d_mod
  implicit none
  type(nasa7_species), allocatable :: species(:), eos(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_1d_config) :: config
  type(mpi_domain_1d) :: domain
  type(amr_patch_level_plan_1d), allocatable :: plans(:)
  type(amr_patch_tree_reactive_solution_1d) :: tree, saved_tree
  type(mpi_amr_patch_distribution_1d) :: distribution, new_distribution
  type(mpi_amr_sparse_reactive_solution_1d) :: sparse, saved_sparse
  real(dp), allocatable :: global_state(:,:), global_temp(:), state(:,:), temperature(:), saved_state(:,:), saved_temp(:)
  real(dp) :: dx,time_step,wave,diffusivity,accepted
  integer :: ierr,rank,variant,operation,rejected,nvar,level,patch,tagged,transferred
  logical :: ok, all_pass, changed
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'load thermo')
  call load_h2o2_full_mechanism(reactions,ok)
  call require(ok,'load reactions')
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'load transport')
  config=reactive_1d_config()
  config%nx=8
  config%problem='uniform_reactor'
  config%initial_temperature=1000
  config%hotspot_temperature_rise=0
  config%chemistry_enabled=.false.
  config%transport_enabled=.false.
  config%chemistry_model='full_h2o2'
  call initialize_reactive_1d(species,config,global_state,global_temp,dx,ok)
  call require(ok,'initialize regular')
  call initialize_mpi_domain_1d(domain,8,MPI_COMM_WORLD,ok)
  call require(ok,'initialize domain')
  nvar=size(global_state,1)
  allocate(state(nvar,0:domain%local_cells+1),temperature(0:domain%local_cells+1))
  state=spread(global_state(:,1),2,domain%local_cells+2)
  temperature=1000
  saved_state=state;saved_temp=temperature
  config%amr_enabled=.true.
  config%amr_minimum_patch_cells=4
  allocate(plans(0))
  call initialize_patch_tree_reactive_1d(species,config,plans,tree,ok)
  call require(ok,'initialize tree')
  call initialize_mpi_amr_patch_distribution_1d(tree%hierarchy,MPI_COMM_WORLD,distribution,ok)
  call require(ok,'initialize distribution')
  call scatter_owned_patch_tree_reactive_1d(distribution,tree,sparse,ok)
  call require(ok,'scatter sparse')
  saved_tree=tree;saved_sparse=sparse
  eos=species
  call configure_nasa7_eos_policy(eos,.true.,ok,200._dp,3500._dp)
  call require(ok,'configure active')
  call mpi_reactive_timestep(domain,eos,transport,state,temperature,dx,.2_dp,.2_dp, &
    .false.,.false.,.false.,.false.,time_step,wave,diffusivity,ok)
  call require(ok.and.time_step>0,'matching active regular control')
  call sparse_patch_tree_reactive_timestep_1d(eos,config,distribution,sparse,time_step,ok)
  call require(ok.and.time_step>0,'matching active sparse control')
  do variant=1,3
    eos=species
    call configure_nasa7_eos_policy(eos,.true.,ok,200._dp,3500._dp)
    call require(ok,'configure variant')
    if(rank==0) then
      select case(variant)
      case(1)
        call configure_nasa7_eos_policy(eos,.false.,ok)
      case(2)
        call configure_nasa7_eos_policy(eos,.true.,ok,200._dp,4000._dp)
      case(3)
        eos(1)%pinned_eos_extrapolation=.false.
      end select
    end if
    call agree_mpi_reactive_eos_1d(eos,MPI_COMM_WORLD,ok)
    call require(.not.ok,'collective policy refusal')
    do operation=1,8
      select case(operation)
      case(1)
        call mpi_reactive_timestep(domain,eos,transport,state,temperature,dx,.2_dp,.2_dp, &
          .false.,.false.,.false.,.false.,time_step,wave,diffusivity,ok)
      case(2)
        call advance_mpi_reactive_hydro(domain,eos,state,temperature,dx,1e-10_dp,'rusanov',ok)
      case(3)
        call advance_mpi_reactive_chemistry(domain,eos,reactions,state,temperature,1e-10_dp,1e-6_dp,1e-12_dp,ok)
      case(4)
        call advance_mpi_reactive_strang(domain,eos,reactions,transport,state,temperature,dx,1e-10_dp, &
          .true.,1e-6_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,'rusanov',ok)
      case(5)
        call advance_mpi_reactive_strang_adaptive(domain,eos,reactions,transport,state,temperature,dx, &
          1e-10_dp,1e-14_dp,.true.,1e-6_dp,1e-12_dp,.true.,.true.,.true.,.true.,.true.,'rusanov', &
          accepted,rejected,ok)
      case(6)
        call mpi_reactive_transport_timestep(domain,eos,transport,state,temperature,dx,.2_dp, &
          .true.,.true.,.true.,time_step,diffusivity,ok)
      case(7)
        call advance_mpi_reactive_transport(domain,eos,transport,state,temperature,dx,1e-10_dp, &
          .true.,.true.,.true.,.true.,ok)
      case(8)
        call advance_mpi_reactive_transport_adaptive(domain,eos,transport,state,temperature,dx, &
          1e-10_dp,1e-14_dp,.true.,.true.,.true.,.true.,accepted,rejected,ok)
      end select
      call require(.not.ok.and.all(state==saved_state).and.all(temperature==saved_temp),'regular refusal and rollback')
    end do
    do operation=1,7
      select case(operation)
      case(1)
        call sparse_patch_tree_reactive_timestep_1d(eos,config,distribution,sparse,time_step,ok)
      case(2)
        call advance_sparse_patch_tree_chemistry_1d(eos,reactions,config,1e-10_dp,distribution,sparse,ok)
      case(3)
        call advance_sparse_patch_tree_hydro_1d(eos,config,1e-10_dp,distribution,sparse,ok)
      case(4)
        call advance_sparse_patch_tree_transport_1d(eos,transport,config,1e-10_dp,distribution,sparse,ok)
      case(5)
        call advance_sparse_patch_tree_reactive_1d(eos,reactions,config,1e-10_dp,distribution,sparse,ok,transport)
      case(6)
        call regrid_sparse_patch_tree_reactive_1d(eos,config,plans,distribution,sparse,new_distribution, &
          changed,transferred,ok)
      case(7)
        call regrid_tagged_sparse_patch_tree_reactive_1d(eos,config,distribution,sparse,new_distribution, &
          changed,tagged,transferred,ok)
      end select
      call require(.not.ok.and.sparse%time==saved_sparse%time.and.sparse%steps==saved_sparse%steps, &
        'sparse refusal and clocks')
      do level=1,size(sparse%levels)
        do patch=1,size(sparse%levels(level)%patches)
          if(.not.sparse%levels(level)%is_local(patch)) cycle
          call require_local(all(sparse%levels(level)%patches(patch)%state== &
            saved_sparse%levels(level)%patches(patch)%state).and. &
            all(sparse%levels(level)%patches(patch)%temperature== &
            saved_sparse%levels(level)%patches(patch)%temperature),'sparse state rollback')
        end do
      end do
      call require(all(sparse%level_advances==saved_sparse%level_advances).and. &
        all(sparse%transport_level_advances==saved_sparse%transport_level_advances),'sparse counter rollback')
    end do
    do operation=1,4
      select case(operation)
      case(1)
        call advance_owned_patch_tree_chemistry_1d(eos,reactions,config,1e-10_dp,distribution,tree,ok)
      case(2)
        call advance_owned_patch_tree_hydro_1d(eos,config,1e-10_dp,distribution,tree,ok)
      case(3)
        call advance_owned_patch_tree_transport_1d(eos,transport,config,1e-10_dp,distribution,tree,ok)
      case(4)
        call advance_owned_patch_tree_reactive_1d(eos,reactions,config,1e-10_dp,distribution,tree,ok,transport)
      end select
      call require(.not.ok.and.tree%time==saved_tree%time.and.tree%steps==saved_tree%steps,'replicated refusal and clocks')
      do level=1,tree%level_count()
        do patch=1,size(tree%levels(level)%patches)
          call require(all(tree%levels(level)%patches(patch)%state==saved_tree%levels(level)%patches(patch)%state).and. &
            all(tree%levels(level)%patches(patch)%temperature==saved_tree%levels(level)%patches(patch)%temperature), &
            'replicated state rollback')
        end do
      end do
    end do
  end do
  if(rank==0) print *, 'test_mpi_reactive_eos_policy: PASS; 57 API refusals and 3 consensus refusals; 2 active controls'
  call MPI_Finalize(ierr)
contains
  subroutine require_local(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    if(condition) return
    print *,rank,trim(label)
    call MPI_Abort(MPI_COMM_WORLD,9,ierr)
  end subroutine require_local
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(len=*),intent(in)::label
    call MPI_Allreduce(condition,all_pass,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    call require_local(ierr==MPI_SUCCESS.and.all_pass,label)
  end subroutine require
end program test_mpi_reactive_eos_policy
