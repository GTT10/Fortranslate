program test_les_mpi_amr_tree
  use mpi_f08
  use les_mpi_amr_tree_mod, only: advance_mpi_amr_tree_les_3d,compute_mpi_amr_tree_les_timestep_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_3d,compute_amr_patch_tree_les_timestep_3d
  use les_smagorinsky_3d_mod, only: advance_les_ssprk2_3d,compute_les_rhs_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,saved,saved_t
  type(amr_patch_3d) :: patches(5)
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),serial(12,4,4,4),st(4,4,4),q(12),sound
  real(dp) :: spacing(3),dt,budget(12),serial_budget(12),largest,initial(12),final(12),wall_t(6),wall_v(3,6),outlet(6)
  real(dp) :: chosen,diffusivity,reference_diff,reference_rhs(12,4,4,4),half_dt,half_diff
  real(dp) :: inlet(12,6),inlet_t(6)
  integer :: i,j,k,mode,geometry,node,m
  character(len=32) :: faces(6),model
  type(amr_tree_state_3d) :: local,local_t,before,before_t
  integer, allocatable :: owners(:)
  integer :: rank,nranks,ierr,owner_mode
  real(dp) :: mpi_budget(12),mpi_largest,trial_ci
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  spacing=[.001_dp,.0015_dp,.002_dp];dt=1e-9_dp;wall_t=1100
  do k=1,4
    do j=1,4
      do i=1,4
        q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=101325*(1+.02_dp*k)
        q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
        q(6)=.025_dp+.003_dp*i;q(9)=.22_dp;q(12)=1-q(6)-q(9)
        call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
        call require(ok,'initial state')
      end do
    end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,2,3,1,4,2,patches(1),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,2,3,1,4,2,patches(2),ok)
  call initialize_amr_patch_3d(4,4,8,3,4,1,2,1,2,2,patches(3),ok)
  call initialize_amr_patch_3d(4,4,8,1,2,1,2,1,2,2,patches(4),ok)
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(5),ok)
  inlet=spread(root(:,1,1,1),2,6);inlet_t=rt(1,1,1,1)
  do geometry=1,2
    if (geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3,4],patches,tree,ok)
    end if
    call require(ok,'topology')
    do m=1,4
      select case(m)
      case(1)
        model='smagorinsky'
      case(2)
        model='wale'
      case(3)
        model='vreman'
      case(4)
        model='dynamic'
      end select
      do mode=1,4
        if(mode==4.and.m/=4) cycle
        faces='periodic';wall_v=0;outlet=0
        if (mode>=2) faces='outflow'
        if (mode==3) then
          faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
          outlet(2)=102325
        end if
        if(m==4.and.mode==3) then
          faces='outflow';faces(1)='inflow';outlet=0;wall_v=0
        end if
        if(mode==4) then
          faces(1:2)='slip_wall';faces(3:4)='no_slip_wall';faces(5:6)='isothermal_wall'
          wall_v(1,3:6)=2;outlet=0
        end if
        call prolong_amr_patch_tree_3d(tree,root,state,ok)
        call require(ok,'state initialization')
        call prolong_amr_patch_tree_3d(tree,rt,temperature,ok)
        call require(ok,'temperature initialization')
        call composite_integrals_amr_patch_tree_3d(tree,state,spacing,initial,ok)
        call require(ok,'initial inventory')
        saved=state;saved_t=temperature
        call advance_amr_patch_tree_les_3d(tree,species,state,temperature,spacing,dt,faces, &
          .17_dp,.7_dp,.7_dp,budget,largest,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
        call require(ok.and.largest>=0,'AMR LES step')
        if(m==4.and.mode>1) call require(maxval(abs(budget))>1e-25_dp,'dynamic open boundary transfers inventory')
        call require(any(state%fields(1)%values/=saved%fields(1)%values),'AMR LES evolves')
        call composite_integrals_amr_patch_tree_3d(tree,state,spacing,final,ok)
        call require(ok,'final inventory')
        call require(maxval(abs(final-initial-budget)/max(1e-12_dp,abs(initial)))<2e-12_dp,'composite conservation')
        do owner_mode=1,2
          if (allocated(owners)) deallocate(owners)
          allocate(owners(size(tree%nodes)))
          owners=0
          if (owner_mode==2) then
            do node=1,size(tree%nodes)
              owners(node)=mod(node-1,nranks)
            end do
          end if
          if (allocated(local%fields)) deallocate(local%fields,local_t%fields)
          allocate(local%fields(size(tree%nodes)),local_t%fields(size(tree%nodes)))
          do node=1,size(tree%nodes)
            if (owners(node)/=rank) cycle
            local%fields(node)%values=saved%fields(node)%values
            local_t%fields(node)%values=saved_t%fields(node)%values
          end do
          call compute_amr_patch_tree_les_timestep_3d(tree,species,saved,saved_t,spacing,.05_dp,faces, &
            .17_dp,.7_dp,.7_dp,chosen,diffusivity,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
            coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
          call require(ok,'serial timestep oracle')
          call compute_mpi_amr_tree_les_timestep_3d(tree,owners,species,local,local_t,spacing,.05_dp,faces, &
            .17_dp,.7_dp,.7_dp,half_dt,half_diff,MPI_COMM_WORLD,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
            coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
          call require(ok,'MPI timestep')
          call require(abs(half_dt-chosen)<=2e-12_dp*chosen,'MPI timestep parity')
          call require(abs(half_diff-diffusivity)<=2e-12_dp*diffusivity,'MPI diffusivity parity')
          half_dt=-7;half_diff=-9;trial_ci=.05_dp
          if (rank==nranks-1) trial_ci=.1_dp
          call compute_mpi_amr_tree_les_timestep_3d(tree,owners,species,local,local_t,spacing,trial_ci,faces, &
            .17_dp,.7_dp,.7_dp,half_dt,half_diff,MPI_COMM_WORLD,ok,model=model)
          call require(.not.ok.and.half_dt==-7.and.half_diff==-9,'CFL disagreement preserves outputs')
          call advance_mpi_amr_tree_les_3d(tree,owners,species,local,local_t,spacing,dt,faces, &
            .17_dp,.7_dp,.7_dp,mpi_budget,mpi_largest,MPI_COMM_WORLD,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
            coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
          call require(ok,'MPI AMR LES')
          do node=1,size(tree%nodes)
            if (owners(node)/=rank) then
              call require(.not.allocated(local%fields(node)%values),'nonowner remains empty')
              cycle
            end if
            call require(maxval(abs(local%fields(node)%values-state%fields(node)%values)/ &
              max(1._dp,abs(state%fields(node)%values)))<2e-12_dp,'MPI state parity')
            call require(maxval(abs(local_t%fields(node)%values-temperature%fields(node)%values))<1e-8_dp, &
              'MPI temperature parity')
          end do
          call require(maxval(abs(mpi_budget-budget))<1e-14_dp,'MPI boundary parity')
          call require(abs(mpi_largest-largest)<1e-12_dp*max(1._dp,largest),'MPI correction parity')
          before=local;before_t=local_t;mpi_budget=-17;mpi_largest=-19
          if(m==4) then
            call advance_mpi_amr_tree_les_3d(tree,owners,species,local,local_t,spacing,dt,faces, &
              .17_dp,.7_dp,.7_dp,mpi_budget,mpi_largest,MPI_COMM_WORLD,ok,model=model,test_ratio=2+2*rank)
            call require(.not.ok,'rank dynamic filter ratio disagreement rejected')
            call check_rollback()
          end if
          trial_ci=.1_dp
          if (rank==nranks-1) trial_ci=.2_dp
          call advance_mpi_amr_tree_les_3d(tree,owners,species,local,local_t,spacing,dt,faces, &
            .17_dp,.7_dp,.7_dp,mpi_budget,mpi_largest,MPI_COMM_WORLD,ok,model=model,ci=trial_ci)
          call require(.not.ok,'rank CI disagreement rejected')
          call check_rollback()
          call advance_mpi_amr_tree_les_3d(tree,owners,species,local,local_t,spacing,1._dp,faces, &
            .17_dp,.7_dp,.7_dp,mpi_budget,mpi_largest,MPI_COMM_WORLD,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
            coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
          call require(.not.ok,'MPI interval rejected')
          call check_rollback()
        end do
        saved=state;saved_t=temperature;budget=-17;largest=-19
        call advance_amr_patch_tree_les_3d(tree,species,state,temperature,spacing,1._dp,faces, &
          .17_dp,.7_dp,.7_dp,budget,largest,ok,wall_temperature=wall_t,wall_velocity=wall_v, &
          coarse_limiter='mc',outlet_pressure=outlet,model=model,ci=.1_dp,inflow_state=inlet,inflow_temperature=inlet_t)
        call require(.not.ok,'excessive interval rejected')
        do node=1,size(tree%nodes)
          call require(all(state%fields(node)%values==saved%fields(node)%values),'state rollback')
          call require(all(temperature%fields(node)%values==saved_t%fields(node)%values),'temperature rollback')
        end do
        call require(all(budget==-17).and.largest==-19,'diagnostics rollback')
      end do
    end do
  end do
  if (rank==0) print *, 'PASS: distributed AMR LES serial parity and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_rollback()
    integer :: p
    do p=1,size(tree%nodes)
      if (owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==before%fields(p)%values),'MPI state rollback')
      call require(all(local_t%fields(p)%values==before_t%fields(p)%values),'MPI temperature rollback')
    end do
    call require(all(mpi_budget==-17).and.mpi_largest==-19,'MPI diagnostic rollback')
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (condition) return
    print *, 'FAIL: ',label
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
end program
