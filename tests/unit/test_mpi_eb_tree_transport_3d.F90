program test_mpi_eb_tree_transport_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d, &
    extract_eb_geometry_3d,build_level_set_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_transport_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_timestep_3d_mod, only: mpi_eb_tree_timestep_3d
  use eb_reactive_patch_tree_3d_mod, only: reactive_eb_tree_timestep_3d
  use mpi_eb_tree_transport_3d_mod, only: advance_mpi_eb_tree_transport_3d,compute_sdc_mpi_eb_tree_transport_3d
  use eb_reactive_patch_tree_3d_mod, only: compute_sdc_eb_tree_transport_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d,fill_sdc_mpi_eb_tree_source_halo_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use eb_reactive_patch_tree_3d_mod, only: compute_sdc_eb_tree_transport_halos_3d
  use mpi_eb_tree_transport_3d_mod, only: compute_sdc_mpi_eb_tree_transport_halos_3d
  use eb_reactive_patch_tree_3d_mod, only: build_eb_tree_extended_geometry_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,compute_state_redist_level_stage_3d
  use mpi_eb_tree_redistribution_3d_mod, only: compute_state_redist_mpi_level_stage_3d
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d
  use mpi_eb_tree_hydro_3d_mod, only: compute_sdc_mpi_eb_tree_hydro_3d, &
    compute_godunov_state_redist_mpi_eb_tree_hydro_3d
  use eb_reactive_patch_tree_3d_mod, only: compute_sdc_eb_tree_hydro_3d, &
    compute_godunov_state_redist_eb_tree_hydro_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(rotating_frame_3d), allocatable :: frame
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:)
  real(dp) :: residual,theta,speed,position(3),potential
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,rank,nranks,status
  character(len=16) :: test_mode
  character(len=16), allocatable :: labels(:)
  real(dp), allocatable :: wall_t(:)
  integer :: rejected_stage=0
  real(dp) :: plane_cut
  logical :: ok
  call get_command_argument(1,test_mode)
  plane_cut=.31_dp
  if(test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or. &
    test_mode=='godunov_halo'.or.test_mode=='active_ctu') plane_cut=.39_dp
  if(test_mode=='isothermal') then
    allocate(labels(6),wall_t(6));labels='isothermal_wall';wall_t=2000
  end if
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  if(test_mode=='physical_boundary') then
    call check_physical_boundary_stage()
    call MPI_Finalize(status)
    stop
  end if
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  do axis=1,3
  do count=4,5
  if(allocated(labels).and.(axis/=3.or.count/=5)) cycle
  if((test_mode=='rotating_dt'.or.(test_mode=='active_ctu'.or.test_mode=='godunov_halo'.or.test_mode=='spatial_halo'.or. &
    test_mode=='constant_halo'.or.test_mode=='constant_sdc'.or. &
    test_mode=='prepared_sdc'.or. &
    (test_mode=='multi_sdc'.or.test_mode=='inert_iters'))).and.(axis/=3.or.count/=5)) cycle
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if(test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or. &
    test_mode=='godunov_halo'.or.test_mode=='active_ctu') then
    ! Put the coarse cut donor immediately below the refined regular cells.
    ! The original boxes covered both members of its vertical neighborhood.
    call initialize_amr_patch_3d(4,4,4,1,2,1,2,3,4,2,boxes(1),ok)
    call require(ok,'coarse-interface child left')
    call initialize_amr_patch_3d(4,4,4,3,4,1,2,3,4,2,boxes(2),ok)
    call require(ok,'coarse-interface child right')
  endif
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),plane_cut,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
    tree%nodes(node)%state(5,:,:,:)=1.1_dp*tree%nodes(node)%state(5,:,:,:)
    do k=1,4
    do j=1,4
    do i=1,4
      if (tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
      call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k),t(i,j,k),primitive, &
        tree%nodes(node)%temperature(i,j,k),speed,ok)
      call require(ok,'temperature perturbation EOS')
    end do
    end do
    end do
  end do
  if(test_mode=='rotating_dt') then
    allocate(frame)
    frame%axis=1;frame%angular_speed=100._dp;frame%axis_point=[.5_dp,.5_dp,.5_dp]
    do node=1,count
    do k=1,4
    do j=1,4
    do i=1,4
      if(tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
      associate(g=>tree%nodes(node)%geometry)
      position=[g%x_lower,g%y_lower,g%z_lower]+(real([i,j,k],dp)-.5_dp)*[g%dx,g%dy,g%dz]
      end associate
      call rotating_potential_3d(frame,position,potential,ok);call require(ok,'rotating cell potential')
      tree%nodes(node)%state(5,i,j,k)=tree%nodes(node)%state(5,i,j,k)- &
        tree%nodes(node)%state(1,i,j,k)*potential
    end do
    end do
    end do
    end do
  end if
  call average_down_reactive_eb_tree_3d(species,tree,ok,frame)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  if(test_mode=='active_ctu') then
    call check_active_ctu(tree)
    deallocate(u,t)
    cycle
  endif
  if((test_mode=='godunov_halo'.or.test_mode=='spatial_halo'.or.test_mode=='constant_halo'.or.test_mode=='constant_sdc'.or. &
    test_mode=='prepared_sdc'.or. &
    (test_mode=='multi_sdc'.or.test_mode=='inert_iters'))) then
    call check_transport_spatial_halos()
    deallocate(u,t)
    cycle
  end if
  if(test_mode=='rotating_dt') then
    call check_rotating_dt()
    deallocate(u,t,frame)
    cycle
  end if
  if(axis==3.and.count==5) call check_sdc_transport()
  boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=labels,wall_temperature=wall_t)
  call require(ok,'tree transport step')
  if(allocated(labels)) then
    call require(boundary(1)==0.and.all(boundary(6:)==0).and.boundary(5)>0,'heated impermeable walls')
  else
    call require(all(boundary==0._dp),'closed physical diffusion faces')
  end if
  call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'transport conservation')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep transport update')
  call check_mpi_transport()
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_transport_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_physical_boundary_stage()
    type(amr_patch_tree_3d) :: mesh
    type(amr_patch_3d) :: empty(0)
    type(eb_geometry_3d) :: wide(1),physical
    type(amr_tree_state_3d) :: old,rhs,serial_rate,local,lrhs,rate,kept
    type(state_redist_register_3d), allocatable :: serial_regs(:),regs(:),kept_regs(:)
    real(dp) :: phi(0:16,0:16,0:16),x,y,integral(2),expected(2),dt
    integer :: i,j,k,c,order,owner,ownership(1)
    logical :: good
    character(len=32) :: faces(6)
    call require(nranks>=2,'physical boundary MPI rank count')
    call initialize_amr_patch_tree_3d([4,4,4],[1],empty,mesh,good)
    call require(good,'physical boundary root mesh')
    do k=0,16
    do j=0,16
    do i=0,16
      phi(i,j,k)=sum(-1.5_dp+.25_dp*real([i,j,k],dp))-2.3_dp
    end do
    end do
    end do
    call build_level_set_eb_geometry_3d(phi,[-1.5_dp,-1.5_dp,-1.5_dp],[2.5_dp,2.5_dp,2.5_dp],wide(1),good)
    call require(good,'physical boundary oblique geometry')
    call extract_eb_geometry_3d(wide(1),[6,6,6],[10,10,10],physical,good)
    call require(good,'physical boundary valid geometry')
    allocate(old%fields(1),rhs%fields(1),local%fields(1),lrhs%fields(1))
    allocate(old%fields(1)%values(2,16,16,16),rhs%fields(1)%values(2,16,16,16))
    do k=1,16
    do j=1,16
    do i=1,16
      x=-1.5_dp+.25_dp*(i-.5_dp+wide(1)%cell_centroid_x(i,j,k))
      y=-1.5_dp+.25_dp*(j-.5_dp+wide(1)%cell_centroid_y(i,j,k))
      old%fields(1)%values(:,i,j,k)=[2._dp+.03_dp*x+.02_dp*y,3._dp+.01_dp*x-.01_dp*y]
      rhs%fields(1)%values(:,i,j,k)=[.1_dp+.01_dp*x,-.03_dp+.01_dp*y]
    end do
    end do
    end do
    dt=.01_dp;faces='outflow'
    do c=1,2
      expected(c)=sum(rhs%fields(1)%values(c,7:10,7:10,7:10)*physical%volume_fraction)
    end do
    do order=0,4,2
      call compute_state_redist_level_stage_3d(mesh,wide,old,rhs,2,0,dt,serial_rate,serial_regs,good, &
        boundary_faces=faces,max_order=order,register_scale=.5_dp)
      call require(good,'physical boundary serial stage')
      do c=1,2
        integral(c)=sum(serial_rate%fields(1)%values(c,:,:,:)*physical%volume_fraction)
      end do
      call require(all(abs(integral-expected)<=3e-12_dp*max(1._dp,abs(expected))), &
        'physical boundary stage keeps the physical source integral')
      do owner=0,nranks-1
        ownership=owner
        if(rank==owner) then
          local=old;lrhs=rhs
        end if
        call compute_state_redist_mpi_level_stage_3d(mesh,ownership,wide,local,lrhs,2,0,dt,rate,regs, &
          MPI_COMM_WORLD,good,boundary_faces=faces,max_order=order,register_scale=.5_dp)
        call require(good,'physical boundary MPI stage')
        if(rank==owner) then
          call require(all(rate%fields(1)%values==serial_rate%fields(1)%values),'physical boundary exact serial/MPI stage')
          call require(all(local%fields(1)%values==old%fields(1)%values).and. &
            all(lrhs%fields(1)%values==rhs%fields(1)%values),'physical boundary immutable stage inputs')
        else
          call require(.not.allocated(rate%fields(1)%values),'physical boundary remote source remains absent')
        end if
        if(allocated(local%fields(1)%values)) deallocate(local%fields(1)%values)
        if(allocated(lrhs%fields(1)%values)) deallocate(lrhs%fields(1)%values)
      end do
    end do
    ! An isolated corner has less than the target physical fluid volume.
    phi=phi-.5_dp
    call build_level_set_eb_geometry_3d(phi,[-1.5_dp,-1.5_dp,-1.5_dp],[2.5_dp,2.5_dp,2.5_dp],wide(1),good)
    call require(good,'physical boundary insufficient-volume geometry')
    ownership=nranks-1;kept=rate;kept_regs=regs
    if(rank==ownership(1)) then
      local=old;lrhs=rhs
    end if
    call compute_state_redist_mpi_level_stage_3d(mesh,ownership,wide,local,lrhs,2,0,dt,rate,regs, &
      MPI_COMM_WORLD,good,boundary_faces=faces,max_order=2,register_scale=.5_dp)
    call require(.not.good,'physical boundary insufficient volume collective refusal')
    if(rank==ownership(1)) then
      call require(all(rate%fields(1)%values==kept%fields(1)%values),'physical boundary failed stage source rollback')
      call require(regs(1)%scale==kept_regs(1)%scale,'physical boundary failed register scale rollback')
      call require(allocated(regs(1)%valid_cells).eqv.allocated(kept_regs(1)%valid_cells), &
        'physical boundary register valid_cells allocation rollback')
      if(allocated(regs(1)%valid_cells)) then
        call require(all(shape(regs(1)%valid_cells)==shape(kept_regs(1)%valid_cells)), &
          'physical boundary register valid_cells shape rollback')
        call require(all(regs(1)%valid_cells.eqv.kept_regs(1)%valid_cells), &
          'physical boundary register valid_cells value rollback')
      end if
      call require(allocated(regs(1)%coarse_flags).eqv.allocated(kept_regs(1)%coarse_flags), &
        'physical boundary register coarse_flags allocation rollback')
      if(allocated(regs(1)%coarse_flags)) then
        call require(all(shape(regs(1)%coarse_flags)==shape(kept_regs(1)%coarse_flags)), &
          'physical boundary register coarse_flags shape rollback')
        call require(all(regs(1)%coarse_flags==kept_regs(1)%coarse_flags), &
          'physical boundary register coarse_flags value rollback')
      end if
      call require(allocated(regs(1)%fine_ghost).eqv.allocated(kept_regs(1)%fine_ghost), &
        'physical boundary register fine_ghost allocation rollback')
      if(allocated(regs(1)%fine_ghost)) then
        call require(all(shape(regs(1)%fine_ghost)==shape(kept_regs(1)%fine_ghost)), &
          'physical boundary register fine_ghost shape rollback')
        call require(all(regs(1)%fine_ghost.eqv.kept_regs(1)%fine_ghost), &
          'physical boundary register fine_ghost value rollback')
      end if
      call require(allocated(regs(1)%coarse_delta).eqv.allocated(kept_regs(1)%coarse_delta), &
        'physical boundary register coarse_delta allocation rollback')
      if(allocated(regs(1)%coarse_delta)) then
        call require(all(shape(regs(1)%coarse_delta)==shape(kept_regs(1)%coarse_delta)), &
          'physical boundary register coarse_delta shape rollback')
        call require(all(regs(1)%coarse_delta==kept_regs(1)%coarse_delta), &
          'physical boundary register coarse_delta value rollback')
      end if
      call require(allocated(regs(1)%fine_delta).eqv.allocated(kept_regs(1)%fine_delta), &
        'physical boundary register fine_delta allocation rollback')
      if(allocated(regs(1)%fine_delta)) then
        call require(all(shape(regs(1)%fine_delta)==shape(kept_regs(1)%fine_delta)), &
          'physical boundary register fine_delta shape rollback')
        call require(all(regs(1)%fine_delta==kept_regs(1)%fine_delta), &
          'physical boundary register fine_delta value rollback')
      end if
    end if
    if(rank==0) print *, 'PASS: physical StateRedist stages orders0/2/4, owner layouts, source conservation and rollback'
  end subroutine check_physical_boundary_stage
  ! Require an actually limited CTU face before checking fresh ML registers.
  ! The tiny-dt spatial-halo tests exercise the same path with inactive limiting.
  subroutine check_active_ctu(fixture)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: geo(:),wide(:),redist(:)
    type(amr_tree_state_3d) :: predictor,raw,corrected,halo,old,expected,local,lt,lp,actual,raw_halo,raw_old
    type(state_redist_register_3d), allocatable :: expected_regs(:),regs(:),raw_regs(:)
    type(eb_integrated_flux_3d), allocatable :: raw_flux(:),corrected_flux(:),flux(:),baseline(:)
    real(dp) :: interval,change,scale,budget,sink,volume
    integer :: p,owners(count),x,y,z
    logical :: good
    character(len=32) :: faces(6)
    faces='outflow'
    allocate(geo(count),redist(count),predictor%fields(count))
    do p=1,count
      geo(p)=fixture%nodes(p)%geometry
      allocate(predictor%fields(p)%values(nv,4,4,4),source=0._dp)
    enddo
    call build_axis_plane_eb_geometry_3d(96,96,96,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp, &
      -2.5_dp,3.5_dp,'z',plane_cut,padded,good)
    call require(good,'active CTU padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,geo,padded,wide,good,halo_width=10)
    call require(good,'active CTU extended geometry')
    do p=1,count
      call extract_eb_geometry_3d(wide(p),[5,5,5],[20,20,20],redist(p),good)
      call require(good,'active CTU register geometry')
    enddo
    interval=1.e-6_dp
    ! Deplete one constituent through an explicit density/species source.
    ! Use an interior cut donor so physical outflow extrapolation cannot
    ! apply a valid-cell allowance to a different exterior PCM divergence.
    ! The source leaves a positive known PCM forward budget. This isolates
    ! correction activity from a CFL violation or invalid low-order state.
    call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,wide,interval,raw,raw_flux,good, &
      boundary_faces=faces,target_depth=0,baseline_flux=baseline)
    call require(good,'active CTU source-free baseline')
    associate(g=>geo(1),f=>baseline(1))
    do z=1,4
    do y=1,4
    do x=1,4
      if(any([x,y,z]/=[2,2,2])) cycle
      if(g%volume_fraction(x,y,z)<=0) cycle
      volume=g%dx*g%dy*g%dz*g%volume_fraction(x,y,z)
      budget=fixture%nodes(1)%state(6,x,y,z)+( &
        g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(6,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(6,x,y,z))+ &
        g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(6,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(6,x,y,z))+ &
        g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(6,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(6,x,y,z)))/volume
      call require(budget>=0,'active CTU funded PCM constituent')
      sink=-(1-1.e-8_dp)*budget/interval
      predictor%fields(1)%values(6,x,y,z)=sink
      predictor%fields(1)%values(1,x,y,z)=sink
    enddo
    enddo
    enddo
    end associate
    call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,wide,interval,raw,raw_flux,good, &
      boundary_faces=faces,target_depth=0,spatial_rate_halos=halo,old_state_halos=old)
    call require(good,'active CTU depleted-source raw flux')
    raw_halo=halo;raw_old=old
    call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,wide,interval,corrected,corrected_flux,good, &
      boundary_faces=faces,target_depth=0,spatial_rate_halos=halo,old_state_halos=old,correct_forward_flux=.true.)
    call require(good,'active CTU corrected serial source')
    change=maxval(abs(corrected_flux(1)%x-raw_flux(1)%x))+ &
      maxval(abs(corrected_flux(1)%y-raw_flux(1)%y))+maxval(abs(corrected_flux(1)%z-raw_flux(1)%z))
    if(rank==0) print *, 'active CTU dt, accepted face change: ',interval,change
    call require(change>0,'active CTU fixture must limit an accepted face')
    expected=corrected
    call compute_state_redist_level_stage_3d(topology,redist,old,halo,nv,0,interval,expected,expected_regs,good, &
      boundary_faces=faces,max_order=2,register_scale=1._dp)
    call require(good,'active CTU independent fresh registers')
    call compute_state_redist_level_stage_3d(topology,redist,raw_old,raw_halo,nv,0,interval,raw,raw_regs,good, &
      boundary_faces=faces,max_order=2,register_scale=1._dp)
    call require(good,'active CTU independent uncorrected registers')
    change=maxval(abs(expected_regs(1)%coarse_delta-raw_regs(1)%coarse_delta))
    if(rank==0) print *, 'active CTU coarse register change: ',change
    call require(change>0,'active CTU correction changes the coarse register')
    owners=[(modulo(count-p,nranks),p=1,count)]
    allocate(local%fields(count),lt%fields(count),lp%fields(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      local%fields(p)%values=fixture%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=fixture%nodes(p)%temperature
      lp%fields(p)%values=predictor%fields(p)%values
    enddo
    call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,lp, &
      wide,interval,actual,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=0,registers=regs, &
      max_order=2,correct_forward_flux=.true.)
    call require(good,'active CTU MPI correction and redistribution')
    do p=1,count
      if(owners(p)/=rank) cycle
      scale=max(1._dp,maxval(abs(expected%fields(p)%values)))
      call require(maxval(abs(actual%fields(p)%values-expected%fields(p)%values))/scale<1.e-12_dp, &
        'active CTU corrected redistributed source')
      if(topology%nodes(p)%level/=0) cycle
      call require(all(flux(p)%x==corrected_flux(p)%x).and.all(flux(p)%y==corrected_flux(p)%y).and. &
        all(flux(p)%z==corrected_flux(p)%z),'active CTU corrected accepted faces')
      call require(regs(p)%scale==1,'active CTU factor1 register scale')
      if(allocated(regs(p)%coarse_delta)) call require(all(regs(p)%coarse_delta==expected_regs(p)%coarse_delta), &
        'active CTU fresh coarse registers')
      if(allocated(regs(p)%fine_delta)) call require(all(regs(p)%fine_delta==expected_regs(p)%fine_delta), &
        'active CTU fresh fine registers')
    enddo
  end subroutine check_active_ctu
  subroutine check_godunov_halos_mpi(fixture,geo)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: geo(:)
    type(state_redist_register_3d), allocatable :: serial_regs(:),regs(:)
    type(eb_geometry_3d) :: padded
    type(eb_geometry_3d), allocatable :: wide(:),transport_geometry(:)
    type(amr_tree_state_3d) :: predictor,local_predictor,local,lt,rate,halo,old,control,serial_halo,serial_old
    type(amr_tree_state_3d) :: kept,kept_halo,kept_old,provided,source_metrics,kept_provided
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:),baseline(:),serial_baseline(:),kept_flux(:)
    integer :: owners(count),p,layout,depth
    real(dp) :: scale,limiter
    logical :: good
    character(len=32) :: faces(6)
    faces='outflow'
    call build_axis_plane_eb_geometry_3d(96,96,96,-2.5_dp,3.5_dp,-2.5_dp,3.5_dp, &
      -2.5_dp,3.5_dp,'z',plane_cut,padded,good)
    call require(good,'MPI Godunov padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,geo,padded,wide,good,halo_width=10)
    call require(good,'MPI Godunov wide geometry')
    allocate(transport_geometry(count))
    do p=1,count
      call extract_eb_geometry_3d(wide(p),[5,5,5],[20,20,20],transport_geometry(p),good)
      call require(good,'MPI non-MOL diffusion geometry crop')
    enddo
    allocate(predictor%fields(count))
    do p=1,count
      allocate(predictor%fields(p)%values(nv,4,4,4),source=0._dp)
      predictor%fields(p)%values(5,:,:,:)=1000
    enddo
    do layout=0,1
      owners=0
      if(layout==1) owners=[(modulo(count-p,nranks),p=1,count)]
      if(allocated(local%fields)) deallocate(local%fields,lt%fields,local_predictor%fields)
      allocate(local%fields(count),lt%fields(count),local_predictor%fields(count))
      if(allocated(source_metrics%fields)) deallocate(source_metrics%fields)
      allocate(source_metrics%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=fixture%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=fixture%nodes(p)%temperature
        local_predictor%fields(p)%values=predictor%fields(p)%values
        allocate(source_metrics%fields(p)%values(1,18,18,18))
        source_metrics%fields(p)%values(1,:,:,:)=wide(p)%volume_fraction(4:21,4:21,4:21)
      enddo
      call fill_sdc_mpi_eb_tree_source_halo_3d(topology,owners,geo,local_predictor,nv,7,provided, &
        MPI_COMM_WORLD,good,boundary_faces=faces,target_volume_fraction=source_metrics)
      call require(good,'MPI Godunov cached grow7 source preparation')
      do depth=0,1
        call compute_sdc_eb_tree_hydro_3d(species,fixture,predictor,wide,1.e-10_dp,control,serial_flux,good, &
          boundary_faces=faces,target_depth=depth,baseline_flux=serial_baseline, &
          spatial_rate_halos=serial_halo,old_state_halos=serial_old)
        call require(good,'MPI Godunov serial spatial control')
        call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
          wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=depth, &
          baseline_flux=baseline,spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
        call require(good,'MPI Godunov collective spatial evaluation')
        do p=1,count
          if(owners(p)/=rank) then
            call require(.not.allocated(rate%fields(p)%values).and..not.allocated(halo%fields(p)%values).and. &
              .not.allocated(old%fields(p)%values),'MPI Godunov remote fields absent')
            cycle
          endif
          scale=max(1._dp,maxval(abs(control%fields(p)%values)))
          call require(maxval(abs(rate%fields(p)%values-control%fields(p)%values))/scale<1.e-12_dp, &
            'MPI Godunov valid RHS serial control')
          if(topology%nodes(p)%level/=depth) then
            call require(.not.allocated(halo%fields(p)%values).and..not.allocated(old%fields(p)%values), &
              'MPI Godunov inactive halo absent')
            cycle
          endif
          call require(all(halo%fields(p)%values(:,7:10,7:10,7:10)==rate%fields(p)%values), &
            'MPI Godunov exact valid RHS crop')
          scale=max(1._dp,maxval(abs(serial_halo%fields(p)%values)))
          call require(maxval(abs(halo%fields(p)%values-serial_halo%fields(p)%values))/scale<1.e-12_dp, &
            'MPI Godunov full spatial halo control')
          call require(all(old%fields(p)%values==serial_old%fields(p)%values),'MPI Godunov exact old halo control')
          call require(maxval(abs(flux(p)%x-serial_flux(p)%x))+maxval(abs(flux(p)%y-serial_flux(p)%y))+ &
            maxval(abs(flux(p)%z-serial_flux(p)%z))<1.e-13_dp,'MPI Godunov accepted faces control')
          call require(all(baseline(p)%x==serial_baseline(p)%x).and.all(baseline(p)%y==serial_baseline(p)%y).and. &
            all(baseline(p)%z==serial_baseline(p)%z),'MPI Godunov baseline faces control')
        enddo
        call compute_godunov_state_redist_eb_tree_hydro_3d(species,fixture,predictor,wide,1.e-10_dp, &
          control,serial_flux,good,boundary_faces=faces,target_depth=depth,registers=serial_regs,max_order=2, &
          correct_forward_flux=.true.)
        call require(good,'MPI Godunov serial ML composition')
        if(depth==0) call require(maxval(abs(serial_regs(1)%coarse_delta))>0,'MPI Godunov nonzero coarse ML ledger')
        call compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
          wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=depth, &
          registers=regs,max_order=2,correct_forward_flux=.true.,predictor_halo=provided)
        call require(good,'MPI Godunov ML composition')
        do p=1,count
          if(owners(p)/=rank) cycle
          scale=max(1._dp,maxval(abs(control%fields(p)%values)))
          call require(maxval(abs(rate%fields(p)%values-control%fields(p)%values))/scale<1.e-12_dp, &
            'MPI Godunov redistributed RHS serial control')
          if(topology%nodes(p)%level/=depth) cycle
          call require(regs(p)%scale==1,'MPI Godunov factor1 ledger')
          if(allocated(regs(p)%coarse_delta)) call require(all(regs(p)%coarse_delta==serial_regs(p)%coarse_delta), &
            'MPI Godunov exact coarse ML ledger')
          if(allocated(regs(p)%fine_delta)) call require(all(regs(p)%fine_delta==serial_regs(p)%fine_delta), &
            'MPI Godunov exact fine ML ledger')
        enddo
        call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,transport_geometry,1.e-10_dp, &
          [.true.,.true.,.true.,.true.],control,serial_flux,serial_halo,serial_old,limiter,good, &
          boundary_faces=faces,target_depth=depth,registers=serial_regs,max_order=2)
        call require(good,'MPI non-MOL diffusion serial ML control')
        call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,transport_geometry, &
          local,lt,1.e-10_dp,[.true.,.true.,.true.,.true.],rate,flux,halo,old,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=faces,target_depth=depth,registers=regs,max_order=2)
        call require(good,'MPI non-MOL diffusion ML export')
        do p=1,count
          if(owners(p)/=rank) cycle
          scale=max(1._dp,maxval(abs(control%fields(p)%values)))
          call require(maxval(abs(rate%fields(p)%values-control%fields(p)%values))/scale<1.e-12_dp, &
            'MPI non-MOL diffusion redistributed RHS control')
          if(topology%nodes(p)%level/=depth) cycle
          call require(regs(p)%scale==1,'MPI non-MOL diffusion factor1 ledger')
          if(allocated(regs(p)%coarse_delta)) call require(all(regs(p)%coarse_delta==serial_regs(p)%coarse_delta), &
            'MPI non-MOL diffusion exact coarse ledger')
          if(allocated(regs(p)%fine_delta)) call require(all(regs(p)%fine_delta==serial_regs(p)%fine_delta), &
            'MPI non-MOL diffusion exact fine ledger')
        enddo
      enddo
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1, &
        spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
      call require(good,'MPI Godunov canonical peer predictor control')
      kept=rate;kept_halo=halo;kept_old=old;kept_flux=flux
      ! Node2's first x+ ghost is a real cell of its same-level sibling.
      ! A cached population source must not replace that cell's actual rate.
      if(owners(2)==rank) provided%fields(2)%values(5,12,9,9)=provided%fields(2)%values(5,12,9,9)+1.e12_dp
      kept_provided=provided
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1, &
        spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
      call require(good,'MPI Godunov peer source priority')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(rate%fields(p)%values==kept%fields(p)%values),'MPI Godunov authoritative peer RHS')
        call require(all(provided%fields(p)%values==kept_provided%fields(p)%values), &
          'MPI Godunov cached predictor remains immutable')
        if(.not.allocated(halo%fields(p)%values)) cycle
        call require(all(halo%fields(p)%values==kept_halo%fields(p)%values),'MPI Godunov authoritative peer halo')
        call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
          all(flux(p)%z==kept_flux(p)%z),'MPI Godunov authoritative shared faces')
      enddo
      faces(1:2)='periodic'
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=0, &
        spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
      call require(good,'MPI Godunov periodic self predictor control')
      kept=rate;kept_halo=halo;kept_flux=flux
      if(owners(1)==rank) provided%fields(1)%values(5,7,9,10)=provided%fields(1)%values(5,7,9,10)+1.e12_dp
      kept_provided=provided
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=0, &
        spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
      call require(good,'MPI Godunov periodic self source priority')
      if(owners(1)==rank) then
        call require(all(rate%fields(1)%values==kept%fields(1)%values).and. &
          all(halo%fields(1)%values==kept_halo%fields(1)%values),'MPI Godunov exact periodic self source')
        call require(all(flux(1)%x==kept_flux(1)%x).and.all(flux(1)%y==kept_flux(1)%y).and. &
          all(flux(1)%z==kept_flux(1)%z),'MPI Godunov periodic self faces')
        call require(all(provided%fields(1)%values==kept_provided%fields(1)%values), &
          'MPI Godunov periodic cache remains immutable')
        provided%fields(1)%values(5,7,9,10)=provided%fields(1)%values(5,7,9,10)-1.e12_dp
      endif
      faces='outflow'
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1, &
        spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
      call require(good,'MPI Godunov coarse ghost predictor control')
      kept=rate
      ! This y+ ghost lies inside the domain but outside every level1 patch.
      ! Its ancestor-supplied source is where cached virtual/ghost populations
      ! may change the CTU predictor rather than replace an actual peer rate.
      if(owners(2)==rank) provided%fields(2)%values(5,9,12,9)=provided%fields(2)%values(5,9,12,9)+1.e12_dp
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1, &
        spatial_rate_halos=halo,old_state_halos=old,predictor_halo=provided)
      call require(good,'MPI Godunov coarse ghost source retained')
      if(owners(2)==rank) call require(any(rate%fields(2)%values/=kept%fields(2)%values), &
        'MPI Godunov preserves custom source without a same-level recipient')
      kept=rate;kept_halo=halo;kept_old=old;kept_flux=flux
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
        wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1, &
        spatial_rate_halos=halo,old_state_halos=old,correct_forward_flux=rank==0)
      call require(.not.good,'MPI Godunov correction mode mismatch rejected')
      if(rank==0) then
        call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
          wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1,old_state_halos=old)
      else
        call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,geo,local,lt,local_predictor, &
          wide,1.e-10_dp,rate,flux,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=1, &
          spatial_rate_halos=halo,old_state_halos=old)
      endif
      call require(.not.good,'MPI Godunov output presence mismatch rejected')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(rate%fields(p)%values==kept%fields(p)%values),'MPI Godunov RHS collective rollback')
        if(.not.allocated(halo%fields(p)%values)) cycle
        call require(all(halo%fields(p)%values==kept_halo%fields(p)%values).and. &
          all(old%fields(p)%values==kept_old%fields(p)%values),'MPI Godunov halo collective rollback')
        call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
          all(flux(p)%z==kept_flux(p)%z),'MPI Godunov flux collective rollback')
      enddo
    enddo
    call check_nonmol_sdc_mpi(fixture,geo,wide)
  end subroutine check_godunov_halos_mpi
  subroutine check_nonmol_sdc_mpi(fixture,geo,wide)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: geo(:),wide(:)
    type(rotating_frame_3d) :: spin
    type(eb_reactive_tree_3d) :: control
    type(elementary_reaction) :: reactions(0)
    type(amr_tree_state_3d) :: local,lt,kept,kept_t
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:),kept_flux(:)
    type(state_redist_register_3d), allocatable :: regs(:),serial_regs(:),kept_regs(:)
    real(dp) :: b(nv),w(nv),ch(nv),sc(nv),sb(nv),sw(nv),sch(nv),ssc(nv),lim,slim,aux(nv)
    integer :: owners(count),p,layout,depth
    logical :: good
    character(len=32) :: faces(6)
    faces='outflow'
    do layout=0,1
      owners=0
      if(layout==1) owners=[(modulo(count-p,nranks),p=1,count)]
      do depth=0,1
        if(allocated(local%fields)) deallocate(local%fields,lt%fields)
        allocate(local%fields(count),lt%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          local%fields(p)%values=fixture%nodes(p)%state
          allocate(lt%fields(p)%values(1,4,4,4))
          lt%fields(p)%values(1,:,:,:)=fixture%nodes(p)%temperature
        enddo
        control=fixture
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,2,.false., &
          [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,control,wide,faces, &
          sb,sw,sch,ssc,slim,good,target_depth=depth,integrated_flux=serial_flux, &
          redistribution='state_redist',mol_mode=.false.,state_redist_max_order=2,state_redist_registers=serial_regs)
        call require(good,'MPI non-MOL SDC serial control')
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
          1.e-10_dp,2,.false.,[.true.,.true.,.true. , .true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
          faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=depth,integrated_flux=flux, &
          redistribution='state_redist',mol_mode=.false.,state_redist_max_order=2,state_redist_registers=regs)
        call require(good,'MPI non-MOL two-iteration SDC')
        call require(maxval(abs(b-sb))+maxval(abs(w-sw))+maxval(abs(ch-sch))+maxval(abs(sc-ssc))<1.e-12_dp, &
          'MPI non-MOL SDC composite diagnostic control')
        do p=1,count
          if(owners(p)/=rank) then
            call require(.not.allocated(local%fields(p)%values).and..not.allocated(regs(p)%coarse_delta).and. &
              .not.allocated(regs(p)%fine_delta),'MPI non-MOL remote state/register absent')
            cycle
          endif
          call require(maxval(abs(local%fields(p)%values-control%nodes(p)%state)/ &
            max(1._dp,abs(control%nodes(p)%state)))<1.e-12_dp,'MPI non-MOL SDC gas control')
          call require(maxval(abs(lt%fields(p)%values(1,:,:,:)-control%nodes(p)%temperature))<1.e-9_dp, &
            'MPI non-MOL SDC temperature control')
          if(topology%nodes(p)%level/=depth) cycle
          call require(maxval(abs(flux(p)%x-serial_flux(p)%x))+maxval(abs(flux(p)%y-serial_flux(p)%y))+ &
            maxval(abs(flux(p)%z-serial_flux(p)%z))<1.e-13_dp,'MPI non-MOL accepted stage faces control')
          if(allocated(regs(p)%coarse_delta)) call require(all(regs(p)%coarse_delta==serial_regs(p)%coarse_delta), &
            'MPI non-MOL exact all-stage coarse ledger')
          if(allocated(regs(p)%fine_delta)) call require(all(regs(p)%fine_delta==serial_regs(p)%fine_delta), &
            'MPI non-MOL exact accepted fine ledger')
        enddo
      enddo
      kept=local;kept_t=lt;kept_flux=flux;kept_regs=regs
      rejected_stage=0
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
        1.e-10_dp,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
        faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=1,integrated_flux=flux, &
        redistribution='state_redist',mol_mode=.false.,state_redist_max_order=2,state_redist_registers=regs, &
        auxiliary_source=reject_mpi_endpoint,auxiliary_change=aux)
      call require(.not.good.and.rejected_stage==1,'MPI non-MOL first endpoint collective rejection')
      spin%axis=3;spin%angular_speed=100;spin%axis_point=[.5_dp,.5_dp,.5_dp]
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
        1.e-10_dp,2,.false.,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
        faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=1,integrated_flux=flux, &
        redistribution='state_redist',mol_mode=.false.,state_redist_max_order=2,state_redist_registers=regs,frame=spin)
      call require(.not.good,'MPI selected-depth validation retains non-MOL rotating-frame rejection')

      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
          all(lt%fields(p)%values==kept_t%fields(p)%values),'MPI non-MOL gas/temperature rollback')
        if(allocated(regs(p)%coarse_delta)) call require(all(regs(p)%coarse_delta==kept_regs(p)%coarse_delta), &
          'MPI non-MOL coarse ledger rollback')
        if(allocated(regs(p)%fine_delta)) call require(all(regs(p)%fine_delta==kept_regs(p)%fine_delta), &
          'MPI non-MOL fine ledger rollback')
        if(.not.allocated(flux(p)%x)) cycle
        call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
          all(flux(p)%z==kept_flux(p)%z),'MPI non-MOL flux rollback')
      enddo
    enddo
  end subroutine check_nonmol_sdc_mpi
  subroutine check_prepared_sdc_mpi(fixture,geo,wide,faces,spin)
    real(dp), allocatable :: coefficients(:)
    type(eb_reactive_tree_3d), intent(in) :: fixture
    type(eb_geometry_3d), intent(in) :: geo(:),wide(:)
    character(len=*), intent(in) :: faces(6)
    type(rotating_frame_3d), optional, intent(in) :: spin
    type(elementary_reaction) :: reactions(1)
    type(elementary_reaction), allocatable :: loaded(:)
    type(eb_reactive_tree_3d) :: control
    type(amr_tree_state_3d) :: local,lt,kept,kept_t
    type(eb_integrated_flux_3d), allocatable :: flux(:),serial_flux(:),kept_flux(:)
    type(state_redist_register_3d), allocatable :: regs(:),serial_regs(:),kept_regs(:)
    real(dp), allocatable :: heat
    real(dp) :: b(nv),w(nv),ch(nv),sc(nv),sb(nv),sw(nv),sch(nv),ssc(nv),lim,slim,aux(nv)
    integer :: owners(count),p,ownership,depth,body_mode
    integer :: iterations
    logical :: good,chemistry
    if(test_mode=='constant_sdc') coefficients=[1.e-5_dp,.02_dp,3.e-5_dp,2.e-5_dp]
    ! Keep a valid stoichiometric reaction while disabling its rate so the
    ! reacting iteration lifecycle is isolated from chemical evolution.
    call load_h2o2_elementary_mechanism(loaded,good)
    call require(good,'load zero-rate iteration mechanism')
    reactions(1)=loaded(1)
    reactions(1)%forward_rate%pre_exponential=0
    reactions(1)%reversible=.false.
    chemistry=test_mode=='multi_sdc'
    iterations=1
    if((test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or.test_mode=='godunov_halo')) iterations=2
    do body_mode=0,1
    if(body_mode==1) then
      allocate(heat);heat=700
    endif
    do ownership=0,1
      owners=0
      if(ownership==1) owners=[(modulo(count-p,nranks),p=1,count)]
      do depth=0,2
        if(allocated(local%fields)) deallocate(local%fields,lt%fields)
        allocate(local%fields(count),lt%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          local%fields(p)%values=fixture%nodes(p)%state
          allocate(lt%fields(p)%values(1,4,4,4))
          lt%fields(p)%values(1,:,:,:)=fixture%nodes(p)%temperature
        enddo
        control=fixture
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-10_dp,iterations,chemistry, &
          [.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp,control,wide,faces, &
          sb,sw,sch,ssc,slim,good,target_depth=depth,integrated_flux=serial_flux, &
          redistribution='state_redist',mol_mode=.true.,frame=spin,state_redist_max_order=2, &
          state_redist_registers=serial_regs,embedded_wall_temperature=heat,embedded_noslip=body_mode==1, &
          constant_coefficients=coefficients)
        call require(good,'MPI prepared SDC serial control')
        if((test_mode=='multi_sdc'.or.test_mode=='inert_iters').and.depth==0) &
          call require(maxval(abs(serial_regs(1)%coarse_delta))>0, &
          'MPI nonzero coarse register control')
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
          1.e-10_dp,iterations,chemistry,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
          faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=depth,integrated_flux=flux, &
          redistribution='state_redist',mol_mode=.true.,frame=spin,state_redist_max_order=2, &
          state_redist_registers=regs,embedded_wall_temperature=heat,embedded_noslip=body_mode==1, &
          constant_coefficients=coefficients)
        call require(good,'MPI prepared SDC selected-level step')
        call require(maxval(abs(b-sb))+maxval(abs(w-sw))+maxval(abs(ch-sch))+maxval(abs(sc-ssc))<1.e-12_dp, &
          'MPI prepared SDC composite diagnostics control')
        call require(abs(lim-slim)<1.e-14_dp,'MPI prepared SDC limiter control')
        do p=1,count
          if(owners(p)/=rank) then
            call require(.not.allocated(local%fields(p)%values).and..not.allocated(regs(p)%coarse_delta).and. &
              .not.allocated(regs(p)%fine_delta),'MPI prepared SDC no remote state or register')
            cycle
          endif
          call require(maxval(abs(local%fields(p)%values-control%nodes(p)%state)/ &
            max(1._dp,abs(control%nodes(p)%state)))<1.e-12_dp,'MPI prepared SDC gas serial control')
          call require(maxval(abs(lt%fields(p)%values(1,:,:,:)-control%nodes(p)%temperature))<1.e-9_dp, &
            'MPI prepared SDC temperature serial control')
          if(topology%nodes(p)%level/=depth) cycle
          call require(maxval(abs(flux(p)%x-serial_flux(p)%x))+maxval(abs(flux(p)%y-serial_flux(p)%y))+ &
            maxval(abs(flux(p)%z-serial_flux(p)%z))<1.e-13_dp,'MPI prepared SDC stage flux control')
          if(allocated(regs(p)%coarse_delta)) then
            call require(maxval(abs(regs(p)%coarse_delta-serial_regs(p)%coarse_delta))<1.e-13_dp, &
              'MPI prepared SDC coarse register control')
          endif
          if(allocated(regs(p)%fine_delta)) then
            call require(maxval(abs(regs(p)%fine_delta-serial_regs(p)%fine_delta))<1.e-13_dp, &
              'MPI prepared SDC fine register control')
          endif
        enddo
      enddo
    enddo
    if(allocated(heat)) deallocate(heat)
    enddo
    kept=local;kept_t=lt;kept_regs=regs;kept_flux=flux
    b=17;w=19;ch=23;sc=29;aux=31;lim=37
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
      1.e-10_dp,iterations,chemistry,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
      faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=2,integrated_flux=flux, &
      redistribution='state_redist',mol_mode=.true.,frame=spin,state_redist_max_order=2,state_redist_registers=regs, &
      auxiliary_source=reject_mpi_endpoint,auxiliary_change=aux, &
      constant_coefficients=coefficients)
    call require(.not.good,'MPI prepared SDC late owner endpoint rejection')
    call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(sc==29).and.all(aux==31).and.lim==37, &
      'MPI prepared SDC diagnostics rollback')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
        all(lt%fields(p)%values==kept_t%fields(p)%values),'MPI prepared SDC gas and EOS rollback')
      if(allocated(regs(p)%coarse_delta)) then
        call require(all(regs(p)%coarse_delta==kept_regs(p)%coarse_delta),'MPI prepared SDC coarse register rollback')
      endif
      if(allocated(regs(p)%fine_delta)) then
        call require(all(regs(p)%fine_delta==kept_regs(p)%fine_delta),'MPI prepared SDC fine register rollback')
      endif
      if(allocated(flux(p)%x)) then
        call require(all(flux(p)%x==kept_flux(p)%x).and.all(flux(p)%y==kept_flux(p)%y).and. &
          all(flux(p)%z==kept_flux(p)%z),'MPI prepared SDC flux rollback')
      endif
    enddo
    if(rank==0) then
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
        1.e-10_dp,iterations,chemistry,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
        faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=2,integrated_flux=flux, &
        redistribution='state_redist',mol_mode=.true.,state_redist_registers=regs, &
        constant_coefficients=coefficients)
    else
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geo,local,lt,wide, &
        1.e-10_dp,iterations,chemistry,[.true.,.true.,.true.,.true.],[0._dp,0._dp,0._dp],1.e-6_dp,1.e-12_dp, &
        faces,b,w,ch,sc,lim,MPI_COMM_WORLD,good,target_depth=2,integrated_flux=flux, &
        redistribution='state_redist',mol_mode=.true., &
        constant_coefficients=coefficients)
    endif
    call require(.not.good,'MPI prepared SDC collective register-output presence')
  end subroutine check_prepared_sdc_mpi
  subroutine reject_mpi_endpoint(gas,temp,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: gas,temp
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: p
    rejected_stage=max(rejected_stage,stage)
    good=stage==0.or.rank/=0
    if(.not.good) return
    do p=1,size(gas%fields)
      if(allocated(rate%fields(p)%values)) rate%fields(p)%values=0
    enddo
  end subroutine reject_mpi_endpoint
  subroutine check_transport_spatial_halos()
    real(dp), allocatable :: coefficients(:)
    type(amr_tree_state_3d), allocatable :: serial_held,local_held
    type(eb_reactive_tree_3d) :: fixture
    type(eb_geometry_3d), allocatable :: geo(:),wide(:)
    type(eb_geometry_3d) :: padded
    type(amr_tree_state_3d) :: local,lt,rates,halos,old,serial_rates,serial_halos,serial_old, &
      kept_rates,kept_halos,kept_old,redist
    type(eb_integrated_flux_3d), allocatable :: fluxes(:),serial_fluxes(:),kept_fluxes(:)
    type(state_redist_register_3d), allocatable :: regs(:)
    type(rotating_frame_3d), allocatable :: spin
    real(dp), allocatable :: phi
    real(dp) :: q(reactive_nprim(size(species))),pos(3),temp_cell,speed,limiter,serial_limiter, &
      serial_dt,mpi_dt,exchange(nv),serial_exchange(nv),frequency,serial_frequency,kept_exchange(nv),kept_frequency,kept_limiter
    integer, allocatable :: selected
    integer :: owners(count),p,x,y,z,c,ownership,rotation_id,depth,policy_id,refusal_id
    character(len=32) :: faces(6)
    logical :: good
    if(test_mode=='constant_halo'.or.test_mode=='constant_sdc') then
      coefficients=[1.e-5_dp,.02_dp,3.e-5_dp,2.e-5_dp]
      transport%diameter=0
    endif
    fixture=saved
    allocate(geo(count))
    do p=1,count
      geo(p)=fixture%nodes(p)%geometry
    enddo
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',plane_cut,padded,good)
    call require(good,'MPI transport halo padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,geo,padded,wide,good)
    call require(good,'MPI transport halo prepared geometry')
    do rotation_id=0,1
      if(rotation_id==1) then
        allocate(spin,phi)
        spin%axis=3;spin%angular_speed=100;spin%axis_point=[.5_dp,.5_dp,.5_dp]
      endif
      fixture=saved
      do p=1,count
        do z=1,4
        do y=1,4
        do x=1,4
          if(geo(p)%volume_fraction(x,y,z)<=0) cycle
          pos=[geo(p)%x_lower+(x-.5_dp)*geo(p)%dx,geo(p)%y_lower+(y-.5_dp)*geo(p)%dy, &
            geo(p)%z_lower+(z-.5_dp)*geo(p)%dz]
          q=0;q(1)=1+.01_dp*pos(1);q(2:4)=.1_dp*pos;q(5)=250000+1000*pos(2)
          if((test_mode=='multi_sdc'.or.test_mode=='inert_iters'.or.test_mode=='godunov_halo')) then
            q(1)=1+.01_dp*pos(1)**2+.02_dp*pos(3)**2
            q(5)=250000+1000*pos(2)**2+3000*pos(3)**2
          endif
          do c=1,size(species)
            if(species(c)%name=='N2') q(c+5)=.95_dp
            if(species(c)%name=='H2') q(c+5)=.05_dp
          enddo
          if(allocated(spin)) then
            call rotating_potential_3d(spin,pos,phi,good)
            call require(good,'MPI transport halo initial potential')
          endif
          call reactive_primitive_to_conserved(species,q,fixture%nodes(p)%state(:,x,y,z),temp_cell,speed,good)
          call require(good,'MPI transport halo smooth state')
          fixture%nodes(p)%temperature(x,y,z)=temp_cell
        enddo
        enddo
        enddo
      enddo
      call average_down_reactive_eb_tree_3d(species,fixture,good)
      call require(good,'MPI transport halo restriction')
      if(allocated(spin)) then
        do p=1,count
        do z=1,4
        do y=1,4
        do x=1,4
          if(geo(p)%volume_fraction(x,y,z)<=0) cycle
          pos=[geo(p)%x_lower+(x-.5_dp)*geo(p)%dx,geo(p)%y_lower+(y-.5_dp)*geo(p)%dy, &
            geo(p)%z_lower+(z-.5_dp)*geo(p)%dz]
          call rotating_potential_3d(spin,pos,phi,good)
          call require(good,'MPI transport halo post-restriction potential')
          fixture%nodes(p)%state(5,x,y,z)=fixture%nodes(p)%state(5,x,y,z)-fixture%nodes(p)%state(1,x,y,z)*phi
        enddo
        enddo
        enddo
        enddo
      endif
      faces='outflow'
      if(rotation_id==1) faces(1:4)='periodic'
      if(test_mode=='godunov_halo') then
        call check_godunov_halos_mpi(fixture,geo)
        return
      endif
      if(test_mode=='prepared_sdc'.or.test_mode=='constant_sdc'.or. &
    (test_mode=='multi_sdc'.or.test_mode=='inert_iters')) then
        call check_prepared_sdc_mpi(fixture,geo,wide,faces,spin)
        if(allocated(spin)) deallocate(spin,phi)
        cycle
      endif
      do ownership=0,2
        owners=0
        if(ownership==1) owners=[(modulo(p-1,nranks),p=1,count)]
        if(ownership==2) owners=[(modulo(count-p,nranks),p=1,count)]
        allocate(local%fields(count),lt%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          local%fields(p)%values=fixture%nodes(p)%state
          allocate(lt%fields(p)%values(1,4,4,4))
          lt%fields(p)%values(1,:,:,:)=fixture%nodes(p)%temperature
        enddo
        if(allocated(coefficients)) then
          call reactive_eb_tree_timestep_3d(species,transport,fixture,.2_dp,.2_dp,.true., &
            .true.,.true.,.true.,serial_dt,good,boundary_faces=faces,frame=spin, &
            constant_coefficients=coefficients)
          call require(good,'constant serial EB timestep')
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geo,local,lt,.2_dp,.2_dp,.true., &
            .true.,.true.,.true.,mpi_dt,MPI_COMM_WORLD,good,boundary_faces=faces,frame=spin, &
            constant_coefficients=coefficients)
          call require(good.and.abs(mpi_dt-serial_dt)<1.e-13_dp*serial_dt,'constant MPI EB timestep equality')
        endif
        ! Exercise both prepared policies on a nonrotating held-energy state.
        ! Policy zero retains the pre-existing default-path regression.
        do policy_id=0,3
        if(rotation_id/=0.and.policy_id>0) cycle
        if(policy_id>0) then
          allocate(serial_held,local_held)
          allocate(serial_held%fields(count),local_held%fields(count))
          do p=1,count
            allocate(serial_held%fields(p)%values(1,4,4,4),source=0._dp)
            do z=1,4
            do y=1,4
            do x=1,4
              if(geo(p)%volume_fraction(x,y,z)<=0) cycle
              serial_held%fields(p)%values(1,x,y,z)=fixture%nodes(p)%state(5,x,y,z)- &
                sum(fixture%nodes(p)%state(2:4,x,y,z)**2)/(2*fixture%nodes(p)%state(1,x,y,z))
            enddo
            enddo
            enddo
            if(policy_id==3) serial_held%fields(p)%values=0.99_dp*serial_held%fields(p)%values
            if(owners(p)==rank) local_held%fields(p)=serial_held%fields(p)
          enddo
        endif
        do depth=-1,2
          if(allocated(selected)) deallocate(selected)
          if(depth>=0) then
            allocate(selected);selected=depth
          endif
          call compute_sdc_eb_tree_transport_halos_3d(species,transport,fixture,wide,1.e-9_dp, &
            [.true.,.true.,.true.,.true.],serial_rates,serial_fluxes,serial_halos,serial_old,serial_limiter,good, &
            boundary_faces=faces,target_depth=selected,frame=spin,embedded_wall_temperature=700._dp, &
            embedded_noslip=.true.,embedded_exchange=serial_exchange,embedded_frequency=serial_frequency, &
            constant_coefficients=coefficients,internal_energy=serial_held, &
            coarse_linear=policy_id>0,separate_momentum=policy_id==2,primitive_from_total_energy=policy_id==3)
          call require(good,'MPI transport halo serial control')
          call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
            1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
            boundary_faces=faces,target_depth=selected,frame=spin,embedded_wall_temperature=700._dp, &
            embedded_noslip=.true.,embedded_exchange=exchange,embedded_frequency=frequency, &
            constant_coefficients=coefficients,internal_energy=local_held, &
            coarse_linear=policy_id>0,separate_momentum=policy_id==2,primitive_from_total_energy=policy_id==3)
          call require(good,'MPI transport halo collective computation')
          if(policy_id==3.and.depth==-1.and.nranks>1) then
          call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
            1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
            boundary_faces=faces,target_depth=selected,frame=spin,embedded_wall_temperature=700._dp, &
            embedded_noslip=.true.,embedded_exchange=exchange,embedded_frequency=frequency, &
            constant_coefficients=coefficients,internal_energy=local_held, &
            coarse_linear=policy_id>0,separate_momentum=policy_id==2,primitive_from_total_energy=rank==0)
            call require(.not.good,'MPI transport primitive policy disagreement rollback')
          endif
          call require(abs(limiter-serial_limiter)<1.e-14_dp,'MPI transport halo limiter control')
          call require(abs(frequency-serial_frequency)<1.e-12_dp*max(1._dp,abs(serial_frequency)), &
            'MPI transport halo global wall frequency')
          call require(maxval(abs(exchange-serial_exchange)/max(1._dp,abs(serial_exchange)))<1.e-12_dp, &
            'MPI transport halo composite wall energy and impulse')
          do p=1,count
            if(owners(p)/=rank) then
              call require(.not.allocated(halos%fields(p)%values).and..not.allocated(old%fields(p)%values), &
                'MPI transport halo no remote state/source allocation')
              cycle
            endif
            call require(maxval(abs(rates%fields(p)%values-serial_rates%fields(p)%values)/ &
              max(1._dp,abs(serial_rates%fields(p)%values)))<1.e-12_dp,'MPI transport halo valid serial equality')
            if(depth>=0.and.topology%nodes(p)%level/=depth) then
              call require(.not.allocated(halos%fields(p)%values),'MPI transport halo inactive level absent')
              cycle
            endif
            call require(maxval(abs(fluxes(p)%x-serial_fluxes(p)%x)/max(1._dp,abs(serial_fluxes(p)%x)))<1.e-12_dp.and. &
              maxval(abs(fluxes(p)%y-serial_fluxes(p)%y)/max(1._dp,abs(serial_fluxes(p)%y)))<1.e-12_dp.and. &
              maxval(abs(fluxes(p)%z-serial_fluxes(p)%z)/max(1._dp,abs(serial_fluxes(p)%z)))<1.e-12_dp, &
              'MPI transport halo integrated flux serial equality')
            call require(all(lbound(old%fields(p)%values)==1),'MPI transport halo state indices')
            call require(maxval(abs(old%fields(p)%values-serial_old%fields(p)%values)/ &
              max(1._dp,abs(serial_old%fields(p)%values)))<1.e-12_dp,'MPI transport halo old serial equality')
            call require(maxval(abs(halos%fields(p)%values-serial_halos%fields(p)%values)/ &
              max(1._dp,abs(serial_halos%fields(p)%values)))<1.e-12_dp,'MPI transport halo entire source serial equality')
            call require(all(halos%fields(p)%values(:,7:10,7:10,7:10)==rates%fields(p)%values), &
              'MPI transport halo valid source after peer limiting')
          enddo
          if(depth>=0) then
            redist=rates
            call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,old,halos,nv,depth,1.e-9_dp, &
              redist,regs,MPI_COMM_WORLD,good,boundary_faces=faces,max_order=2,register_scale=.5_dp)
            call require(good,'MPI transport halo feeds StateRedist level')
          endif
        enddo
        if(policy_id==2.and.nranks>1) then
          do refusal_id=1,2
            call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
              1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
              boundary_faces=faces,internal_energy=local_held,coarse_linear=refusal_id==1, &
              separate_momentum=rank==0.or.refusal_id==2,constant_coefficients=coefficients)
            call require(.not.good,'MPI transport rejects mismatched or unsupported momentum policy')
            call require(abs(limiter-serial_limiter)<1.e-14_dp,'MPI transport policy refusal preserves limiter')
            do p=1,count
              if(owners(p)/=rank) cycle
              call require(maxval(abs(rates%fields(p)%values-serial_rates%fields(p)%values)/ &
                max(1._dp,abs(serial_rates%fields(p)%values)))<1.e-12_dp, &
                'MPI transport policy refusal preserves rates')
              if(allocated(old%fields(p)%values)) then
                call require(maxval(abs(old%fields(p)%values-serial_old%fields(p)%values)/ &
                  max(1._dp,abs(serial_old%fields(p)%values)))<1.e-12_dp, &
                  'MPI transport policy refusal preserves old halo')
                call require(maxval(abs(halos%fields(p)%values-serial_halos%fields(p)%values)/ &
                  max(1._dp,abs(serial_halos%fields(p)%values)))<1.e-12_dp, &
                  'MPI transport policy refusal preserves source halo')
              endif
            enddo
          enddo
        endif
        if(allocated(serial_held)) deallocate(serial_held,local_held)
        enddo
        kept_rates=rates;kept_halos=halos;kept_old=old;kept_fluxes=fluxes
        kept_exchange=exchange;kept_frequency=frequency;kept_limiter=limiter
        if(allocated(coefficients)) then
          if(rank==nranks-1) coefficients(3)=2*coefficients(3)
          call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
            1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
            boundary_faces=faces,frame=spin,embedded_wall_temperature=700._dp,embedded_noslip=.true., &
            embedded_exchange=exchange,embedded_frequency=frequency, &
            constant_coefficients=coefficients)
          call require(.not.good,'constant MPI EB controls mismatch')
        call require(all(exchange==kept_exchange).and.frequency==kept_frequency.and.limiter==kept_limiter, &
          'MPI transport halo wall and limiter rollback')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(rates%fields(p)%values==kept_rates%fields(p)%values),'MPI transport halo rates rollback')
          call require(allocated(fluxes(p)%x).eqv.allocated(kept_fluxes(p)%x),'MPI transport halo flux allocation rollback')
          if(allocated(fluxes(p)%x)) then
            call require(all(fluxes(p)%x==kept_fluxes(p)%x).and.all(fluxes(p)%y==kept_fluxes(p)%y).and. &
              all(fluxes(p)%z==kept_fluxes(p)%z),'MPI transport halo flux rollback')
          endif
          if(allocated(old%fields(p)%values)) then
            call require(all(old%fields(p)%values==kept_old%fields(p)%values),'MPI transport halo snapshot rollback')
            call require(all(halos%fields(p)%values==kept_halos%fields(p)%values),'MPI transport halo RHS rollback')
          endif
        enddo
          coefficients(3)=3.e-5_dp
        endif
        if(owners(count)==rank) local%fields(count)%values(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
        call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
          1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=faces,frame=spin,embedded_wall_temperature=700._dp,embedded_noslip=.true., &
          embedded_exchange=exchange,embedded_frequency=frequency, &
          constant_coefficients=coefficients)
        call require(.not.good,'MPI transport halo remote EOS failure')
        call require(all(exchange==kept_exchange).and.frequency==kept_frequency.and.limiter==kept_limiter, &
          'MPI transport halo wall and limiter rollback')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(rates%fields(p)%values==kept_rates%fields(p)%values),'MPI transport halo rates rollback')
          call require(allocated(fluxes(p)%x).eqv.allocated(kept_fluxes(p)%x),'MPI transport halo flux allocation rollback')
          if(allocated(fluxes(p)%x)) then
            call require(all(fluxes(p)%x==kept_fluxes(p)%x).and.all(fluxes(p)%y==kept_fluxes(p)%y).and. &
              all(fluxes(p)%z==kept_fluxes(p)%z),'MPI transport halo flux rollback')
          endif
          if(allocated(old%fields(p)%values)) then
            call require(all(old%fields(p)%values==kept_old%fields(p)%values),'MPI transport halo snapshot rollback')
            call require(all(halos%fields(p)%values==kept_halos%fields(p)%values),'MPI transport halo RHS rollback')
          endif
        enddo
        if(owners(count)==rank) local%fields(count)%values=fixture%nodes(count)%state
        if(nranks>1) then
          if(rank==0) then
            call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
              1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
              embedded_wall_temperature=700._dp,embedded_exchange=exchange,embedded_frequency=frequency, &
              constant_coefficients=coefficients)
          else
            call compute_sdc_mpi_eb_tree_transport_halos_3d(species,transport,topology,owners,geo,wide,local,lt, &
              1.e-9_dp,[.true.,.true.,.true.,.true.],rates,fluxes,halos,old,limiter,MPI_COMM_WORLD,good, &
              constant_coefficients=coefficients)
          endif
          call require(.not.good,'MPI transport halo mismatched wall controls')
        endif
        deallocate(local%fields,lt%fields)
      enddo
      if(allocated(spin)) deallocate(spin,phi)
    enddo
  end subroutine check_transport_spatial_halos
  subroutine check_rotating_dt()
    type(amr_tree_state_3d) :: local,lt
    type(eb_geometry_3d) :: geometry(count)
    type(rotating_frame_3d) :: mismatched
    real(dp) :: expected_dt,actual_dt,mach(6),tin(6),wall_temperature(6)
    character(len=16) :: faces(6)
    integer :: owners(count),p,mode,recursive
    logical :: good,enabled
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(local%fields(count),lt%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      if(owners(p)/=rank) cycle
      local%fields(p)%values=saved%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
    end do
    if(nranks>1) then
      mismatched=frame
      if(rank==1) mismatched%angular_speed=mismatched%angular_speed+1._dp
      actual_dt=17._dp
      call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,.35_dp,.true., &
        .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,frame=mismatched)
      call require(.not.good.and.actual_dt==17._dp,'rank-mismatched rotating timestep rollback')
    end if
    mach=-1._dp;mach(1)=.2_dp;tin=1000._dp
    faces='outflow';faces(1)='isothermal_wall';wall_temperature=1500._dp
    do recursive=0,1
    do mode=1,3
      enabled=mode/=1
      select case(mode)
      case(1)
        call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,.35_dp,enabled, &
          .true.,.true.,.true.,expected_dt,good,subcycling=recursive==1,frame=frame)
        call require(good,'serial rotating hydro timestep')
        call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,.35_dp,enabled, &
          .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive==1,frame=frame)
      case(2)
        call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,.35_dp,enabled, &
          .true.,.true.,.true.,expected_dt,good,subcycling=recursive==1,frame=frame, &
          inflow_mach=mach,inflow_temperature=tin)
        call require(good,'serial rotating Mach timestep')
        call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,.35_dp,enabled, &
          .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive==1,frame=frame, &
          inflow_mach=mach,inflow_temperature=tin)
      case(3)
        call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,.35_dp,enabled, &
          .true.,.true.,.true.,expected_dt,good,subcycling=recursive==1,frame=frame, &
          boundary_faces=faces,wall_temperature=wall_temperature)
        call require(good,'serial rotating wall timestep')
        call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,.35_dp,enabled, &
          .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive==1,frame=frame, &
          boundary_faces=faces,wall_temperature=wall_temperature)
      end select
      call require(good.and.actual_dt==expected_dt,'serial/MPI rotating timestep')
    end do
    end do
  end subroutine check_rotating_dt
  subroutine check_sdc_transport()
    type(eb_reactive_tree_3d) :: stage_tree,end_tree
    type(amr_tree_state_3d) :: local,lt,rates,reference,prior,ghost,gt,end_local,level_rates,level_reference
    type(eb_integrated_flux_3d), allocatable :: flux(:),expected(:),old_flux(:),level_flux(:),level_expected(:)
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: limiter,expected_theta,dt,level_theta,level_expected_theta
    integer :: owners(count),p,mode,failure
    logical :: good
    stage_tree=saved
    ! Provisional parents deliberately differ from the restricted fine state.
    stage_tree%nodes(1)%state=1.01_dp*stage_tree%nodes(1)%state
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      stage_tree,reference,expected,expected_theta,good,boundary_faces=labels,wall_temperature=wall_t)
    call require(good,'serial SDC stage')
    do p=1,count
      geometry(p)=stage_tree%nodes(p)%geometry
    end do
    do mode=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(mode==2) owners=0
      allocate(local%fields(count),lt%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=stage_tree%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=stage_tree%nodes(p)%temperature
      end do
      call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1.e-8_dp, &
        .true.,.true.,.true.,.true.,local,lt,rates,flux,limiter,MPI_COMM_WORLD,good, &
        boundary_faces=labels,wall_temperature=wall_t)
      call require(good.and.limiter==expected_theta,'MPI SDC limiter')
      do p=1,count
        call require(allocated(rates%fields(p)%values).eqv.(owners(p)==rank),'owner-local SDC rates')
        call require(allocated(flux(p)%x).eqv.(owners(p)==rank),'owner-local SDC flux')
        if(owners(p)/=rank) cycle
        call require(all(rates%fields(p)%values==reference%fields(p)%values),'exact serial SDC rate')
        call require(all(flux(p)%x==expected(p)%x).and.all(flux(p)%y==expected(p)%y).and. &
          all(flux(p)%z==expected(p)%z),'exact serial SDC flux')
        call require(all(local%fields(p)%values==stage_tree%nodes(p)%state),'SDC preserves input U')
        call require(all(lt%fields(p)%values(1,:,:,:)==stage_tree%nodes(p)%temperature),'SDC preserves input T')
      end do
      end_tree=stage_tree;end_local=local
      do p=1,count
        end_tree%nodes(p)%state=1.01_dp*stage_tree%nodes(p)%state
        if(owners(p)==rank) end_local%fields(p)%values=end_tree%nodes(p)%state
      end do
      call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
        stage_tree,level_reference,level_expected,level_expected_theta,good,boundary_faces=labels,wall_temperature=wall_t, &
        target_depth=2,ancestor_start=stage_tree,ancestor_end=end_tree,ancestor_fraction=[.5_dp,.5_dp])
      call require(good,'serial selected transport stage')
      call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1.e-8_dp, &
        .true.,.true.,.true.,.true.,local,lt,level_rates,level_flux,level_theta,MPI_COMM_WORLD,good, &
        boundary_faces=labels,wall_temperature=wall_t,target_depth=2,ancestor_start=local,ancestor_end=end_local, &
        temperature_start=lt,temperature_end=lt,ancestor_fraction=[.5_dp,.5_dp])
      call require(good.and.level_theta==level_expected_theta,'selected MPI transport limiter')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(level_rates%fields(p)%values==level_reference%fields(p)%values),'selected MPI transport rate')
        if(topology%nodes(p)%level/=2) then
          call require(.not.allocated(level_flux(p)%x),'inactive transport flux absent')
          cycle
        end if
        call require(all(level_flux(p)%x==level_expected(p)%x).and.all(level_flux(p)%y==level_expected(p)%y).and. &
          all(level_flux(p)%z==level_expected(p)%z),'selected MPI transport flux')
      end do
      prior=rates;old_flux=flux
      do failure=1,2
        dt=1.e-8_dp;limiter=19
        if(failure==1.and.rank==owners(count)) local%fields(count)%values(1,4,4,4)=-1
        if(failure==2.and.rank==nranks-1) dt=-dt
        call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,dt, &
          .true.,.true.,.true.,.true.,local,lt,rates,flux,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=labels,wall_temperature=wall_t)
        call require(.not.good.and.limiter==19,'collective SDC rejection')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(rates%fields(p)%values==prior%fields(p)%values),'SDC rates rollback')
          call require(all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
            all(flux(p)%z==old_flux(p)%z),'SDC flux rollback')
        end do
        if(rank==owners(count)) local%fields(count)%values=stage_tree%nodes(count)%state
      end do
      call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gt,MPI_COMM_WORLD, &
        good,restrict_snapshot=rank==0)
      call require(.not.good.and..not.allocated(ghost%fields).and..not.allocated(gt%fields), &
        'rank-dependent restriction rejected before halo exchange')
      deallocate(local%fields,lt%fields)
    end do
  end subroutine check_sdc_transport
  subroutine check_mpi_transport()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: mach(6),tin(6),pin(6),old_dt
    real(dp) :: b(nv),limiter,inventory(nv),error,dt,expected_dt,actual_dt,hcfl,tcfl
    integer :: owners(count),p,mode,failure,estimate_mode,recursive_mode
    logical :: good,visc
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
    end do
    do mode=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if (mode==2) owners=nranks-1
      allocate(local%fields(count),lt%fields(count))
      do p=1,count
        if (owners(p)/=rank) cycle
        local%fields(p)%values=saved%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
      end do
      initial=local;it=lt;b=17._dp;limiter=19._dp
      if(allocated(wall_t)) then
        do recursive_mode=0,2
          call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,old_dt,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2)
          call require(good,'thermal baseline timestep')
          call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,expected_dt,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2, &
            boundary_faces=labels,wall_temperature=wall_t)
          call require(good.and.expected_dt<=old_dt,'serial thermal bound')
          if(recursive_mode==1) call require(expected_dt<old_dt,'thermal bound is active')
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2, &
            boundary_faces=labels,wall_temperature=wall_t)
          call require(good.and.actual_dt==expected_dt,'exact MPI thermal timestep')
          if(rank==0) print *, 'thermal dt baseline, serial, MPI: ',old_dt,expected_dt,actual_dt
          actual_dt=17
          if(rank==0) wall_t(1)=2001
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,1.e-12_dp,.true., &
            .false.,.true.,.false.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2, &
            boundary_faces=labels,wall_temperature=wall_t)
          call require(.not.good.and.actual_dt==17,'rank-dependent wall temperature preserves timestep')
          wall_t=2000
        end do
      end if

      do recursive_mode=0,2
      do estimate_mode=1,3
        hcfl=.5_dp;tcfl=.35_dp
        if (estimate_mode==3) tcfl=1.e-10_dp
        call reactive_eb_tree_timestep_3d(species,transport,saved,hcfl,tcfl,estimate_mode/=1, &
          .true.,.true.,.true.,expected_dt,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2)
        call require(good,'serial estimate')
        call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,estimate_mode/=1, &
          .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2)
        call require(good.and.actual_dt==expected_dt,'exact MPI hydro/diffusion bound')
        if (axis==3.and.count==5) then
          old_dt=expected_dt;mach=10;tin=2000;pin=101325
          call reactive_eb_tree_timestep_3d(species,transport,saved,hcfl,tcfl,estimate_mode/=1, &
            .true.,.true.,.true.,expected_dt,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2, &
            inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
          call require(good.and.expected_dt<old_dt,'hot fast inlet tightens hydro and diffusion bound')
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,estimate_mode/=1, &
            .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2, &
            inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
          call require(good.and.actual_dt==expected_dt,'exact MPI inlet timestep')
          actual_dt=17
          if (rank==0) pin(1)=101326
          call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,estimate_mode/=1, &
            .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=recursive_mode>=1, &
            coupled_subcycling=recursive_mode==2, &
            inflow_mach=mach,inflow_temperature=tin,inflow_pressure=pin)
          call require(.not.good.and.actual_dt==17,'rank-dependent inlet preserves timestep')
        end if
      end do
      end do
      actual_dt=17._dp
      if (rank==0) hcfl=.25_dp
      call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,hcfl,tcfl,.true., &
        .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good)
      call require(.not.good.and.actual_dt==17._dp,'inconsistent CFL preserves estimate')
      call mpi_eb_tree_timestep_3d(species,transport,topology,owners,geometry,local,lt,.5_dp,tcfl,.true., &
        .true.,.true.,.true.,actual_dt,MPI_COMM_WORLD,good,subcycling=rank==0)
      call require(.not.good.and.actual_dt==17._dp,'inconsistent subcycling preserves estimate')
      do p=1,count
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'estimate preserves U')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'estimate preserves T')
      end do
      call advance_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1.e-8_dp, &
        .true.,.true.,.true.,.true.,local,lt,b,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=labels,wall_temperature=wall_t)
      call require(good,'MPI transport')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local U')
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial transport U')
        call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial transport T')
      end do
      call require(all(b==boundary),'exact boundary ledger')
      call require(limiter==theta,'exact serial minimum limiter')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
      call require(good,'MPI inventory')
      error=maxval(abs(inventory-before-b)/max(1._dp,abs(before)))
      call require(error<3.e-14_dp,'MPI transport conservation')
      if (rank==0) print *, 'MPI axis, nodes, mode, residual: ',axes(axis),count,mode,error
      do failure=1,5
        if(failure==5.and..not.allocated(wall_t)) cycle
        local=initial;lt=it;b=17._dp;limiter=19._dp;dt=1.e-8_dp
        visc=.true.
        if (failure==1) dt=-1.e-8_dp
        if (failure==2.and.rank==0) dt=2.e-8_dp
        if (failure==3.and.rank==owners(count)) local%fields(count)%values(1,4,4,4)=-1._dp
        if (failure==4.and.rank==0) visc=.false.
        if(failure==5.and.rank==0) wall_t(1)=2001
        call advance_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,dt, &
          visc,.true.,.true.,.true.,local,lt,b,limiter,MPI_COMM_WORLD,good, &
          boundary_faces=labels,wall_temperature=wall_t)
        if(allocated(wall_t)) wall_t=2000
        call require(.not.good,'collective transport rejection')
        if (failure==3.and.rank==owners(count)) then
          call require(local%fields(count)%values(1,4,4,4)==-1._dp,'bad density preserved')
          local%fields(count)%values(1,4,4,4)=initial%fields(count)%values(1,4,4,4)
        end if
        do p=1,count
          if (owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values),'U rollback')
          call require(all(lt%fields(p)%values==it%fields(p)%values),'T rollback')
        end do
        call require(all(b==17._dp).and.limiter==19._dp,'ledger rollback')
      end do
      deallocate(local%fields,lt%fields)
    end do
  end subroutine check_mpi_transport
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_transport_3d
