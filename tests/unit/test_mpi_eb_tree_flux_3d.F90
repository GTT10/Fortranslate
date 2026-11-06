program test_mpi_eb_tree_flux_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use mixture_thermo_mod, only: mixture_mass_properties,mixture_temperature_bounds
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d, &
    build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,synchronize_reactive_eb_tree_flux_3d, &
    assemble_reactive_eb_tree_reflux_3d,repair_reactive_eb_tree_interface_budget_3d, &
    repair_reactive_eb_tree_held_interface_budget_3d, &
    build_eb_tree_extended_geometry_3d,synchronize_reactive_eb_tree_register_flux_3d, &
    synchronize_reactive_eb_tree_scalar_register_flux_3d,synchronize_reactive_eb_tree_scalar_flux_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,integrate_mpi_eb_tree_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d,assemble_mpi_eb_tree_reflux_3d, &
    repair_mpi_eb_tree_interface_budget_3d,repair_mpi_eb_tree_held_interface_budget_3d, &
    synchronize_mpi_eb_tree_register_flux_3d, &
    synchronize_mpi_eb_tree_scalar_register_flux_3d,synchronize_mpi_eb_tree_scalar_flux_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d,initialize_state_redist_tree_register_3d, &
    assemble_state_redist_interface_rereflux_3d,assemble_state_redist_interface_source_3d, &
    assemble_state_redist_patch_rereflux_3d,state_redist_patch_rereflux_budget_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  type(rotating_frame_3d) :: frame
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  type(eb_integrated_flux_3d), allocatable :: flux(:),raw(:)
  real(dp), allocatable :: before(:),after(:),boundary(:),basis(:)
  real(dp) :: volume,residual,potential,position(3)
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,rank,nranks,status
  character(len=16), allocatable :: faces(:)
  character(len=16) :: test_mode
  logical :: ok,rotating
  logical :: native_query_mode,reject_idle_query=.false.
  integer :: native_query_calls=0
  call get_command_argument(1,test_mode)
  rotating=test_mode=='rotating'
  native_query_mode=test_mode=='native_query'
  frame%axis=2;frame%angular_speed=100._dp;frame%axis_point=[.25_dp,.4_dp,.5_dp]
  if(test_mode=='periodic') allocate(faces(6))
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  if(test_mode=='held_thermfull'.or.test_mode=='gas_thermfull'.or. &
    test_mode=='gas_fine_eos'.or.test_mode=='held_fine_eos'.or. &
    test_mode=='gas_fine_cov'.or.test_mode=='held_fine_cov'.or. &
    test_mode=='gas_raw_eos'.or.test_mode=='held_raw_eos') then
    call load_h2o2_full_thermo(species,ok)
  else
    call load_h2o2_elementary_thermo(species,ok)
  endif
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  if(native_query_mode) then
    call check_held_ancestor_repair(uniform_inventory=.true.,far_endpoint=.true.)
    if(rank==0) print *, 'PASS: collective adjacent native packet query, conserved gas/UEINT and idle refusal rollback'
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_raw_eos'.or.test_mode=='held_raw_eos') then
    do i=1,5
      do j=0,1
        do k=0,1
          call check_held_ancestor_repair(gas_only=test_mode=='gas_raw_eos',thermal_case=i, &
            raw_endpoint=.true.,empty_register=k==1,far_endpoint=.true.,geometry_case=j)
          call check_held_ancestor_repair(gas_only=test_mode=='gas_raw_eos',thermal_case=i, &
            raw_endpoint=.true.,empty_register=k==1,fine_endpoint=.true.,geometry_case=j)
          call check_held_ancestor_repair(gas_only=test_mode=='gas_raw_eos',thermal_case=i, &
            raw_endpoint=.true.,empty_register=k==1,fine_covered=.true.,geometry_case=j)
        enddo
      enddo
    enddo
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_fine_cov'.or.test_mode=='held_fine_cov') then
    do i=1,5
      call check_held_ancestor_repair(gas_only=test_mode=='gas_fine_cov',thermal_case=i,fine_covered=.true.)
      call check_held_ancestor_repair(gas_only=test_mode=='gas_fine_cov',thermal_case=i, &
        fine_covered=.true.,geometry_case=1)
    enddo
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_fine_eos'.or.test_mode=='held_fine_eos') then
    do i=1,5
      call check_held_ancestor_repair(gas_only=test_mode=='gas_fine_eos',thermal_case=i,fine_endpoint=.true.)
      call check_held_ancestor_repair(gas_only=test_mode=='gas_fine_eos',thermal_case=i, &
        fine_endpoint=.true.,geometry_case=1)
    enddo
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_cut') then
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2,geometry_case=2)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='held_cut') then
    call check_held_ancestor_repair(covered_depth=2,geometry_case=2)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_far') then
    call check_held_ancestor_repair(gas_only=.true.,uniform_inventory=.true.,far_endpoint=.true.)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_thermal'.or.test_mode=='gas_thermfull') then
    call check_held_ancestor_repair(gas_only=.true.,far_endpoint=.true.,thermal_case=1)
    call check_held_ancestor_repair(gas_only=.true.,far_endpoint=.true.,thermal_case=2)
    call check_held_ancestor_repair(gas_only=.true.,far_endpoint=.true.,thermal_case=3)
    call check_held_ancestor_repair(gas_only=.true.,far_endpoint=.true.,thermal_case=4)
    call check_held_ancestor_repair(gas_only=.true.,far_endpoint=.true.,thermal_case=5)
    call check_held_ancestor_repair(gas_only=.true.,far_endpoint=.true.,thermal_case=6,geometry_case=3)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='held_thermal'.or.test_mode=='held_thermfull') then
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=1)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=2)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=3)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=4)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=5)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='held_far') then
    call check_held_ancestor_repair(uniform_inventory=.true.,far_endpoint=.true.)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_uniform') then
    call check_held_ancestor_repair(gas_only=.true.,uniform_inventory=.true.)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='held_uniform') then
    call check_held_ancestor_repair(uniform_inventory=.true.)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='gas_ancestor') then
    call check_held_ancestor_repair(gas_only=.true.,uniform_inventory=.true.,far_endpoint=.true.)
    call check_held_ancestor_repair(gas_only=.true.,uniform_inventory=.true.)
    call check_held_ancestor_repair(gas_only=.true.)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=1)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=1,geometry_case=1)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2,geometry_case=1)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=1,geometry_case=2)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2,geometry_case=2)
    call check_ancestor_register_close()
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='held_ancestor') then
    call check_held_ancestor_repair()
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='held_budget') then
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=1)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=2)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=3)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=4)
    call check_held_ancestor_repair(far_endpoint=.true.,thermal_case=5)
    call check_held_ancestor_repair(uniform_inventory=.true.,far_endpoint=.true.)
    call check_held_ancestor_repair(uniform_inventory=.true.)
    call check_held_interface_repair()
    call check_held_ancestor_repair()
    call check_held_ancestor_repair(covered_depth=1)
    call check_held_ancestor_repair(covered_depth=2)
    call check_held_ancestor_repair(covered_depth=1,geometry_case=1)
    call check_held_ancestor_repair(covered_depth=2,geometry_case=1)
    call check_held_ancestor_repair(covered_depth=1,geometry_case=2)
    call check_held_ancestor_repair(covered_depth=2,geometry_case=2)
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='register_close') then
    call check_register_close()
    call check_held_ancestor_repair(gas_only=.true.)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=1)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=1,geometry_case=1)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2,geometry_case=1)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=1,geometry_case=2)
    call check_held_ancestor_repair(gas_only=.true.,covered_depth=2,geometry_case=2)
    call check_ancestor_register_close()
    if(rank==0) print *, 'PASS: MPI interface register close and EOS rollback'
    call MPI_Finalize(status)
    stop
  endif
  allocate(before(nv),after(nv),boundary(nv),basis(nv))
  if(.not.rotating.and.test_mode/="periodic") then
    call check_interface_repair()
    call check_cross_patch_repair()
    call check_peer_coarse_donor()
    call check_refined_interface_repair()
    call check_refined_corridor_repair()
  endif
  do axis=1,3
  if(allocated(faces)) then
    faces='periodic';faces(2*axis-1:2*axis)='outflow'
  end if
  do count=4,5
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4),flux(count))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  if(rotating) then
    do k=1,4
    do j=1,4
    do i=1,4
      if(coarse%volume_fraction(i,j,k)<=0) cycle
      position=[coarse%x_lower,coarse%y_lower,coarse%z_lower]+ &
        (real([i,j,k],dp)-.5_dp)*[coarse%dx,coarse%dy,coarse%dz]
      call rotating_potential_3d(frame,position,potential,ok)
      call require(ok,'root potential')
      u(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
    end do
    end do
    end do
    call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok,frame)
  else
    call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  end if
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  if(rotating) then
    call average_down_reactive_eb_tree_3d(species,tree,ok,frame)
  else
    call average_down_reactive_eb_tree_3d(species,tree,ok)
  end if
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  ! Manufactured time-integrated transfers, deliberately unequal between nodes.
  ! A single physical conserved vector keeps species closure and positive EOS.
  basis=u(:,4,4,4)
  do node=1,count
    allocate(flux(node)%x(nv,0:4,4,4),flux(node)%y(nv,4,0:4,4),flux(node)%z(nv,4,4,0:4))
    do k=1,4
    do j=1,4
    do i=0,4
      flux(node)%x(:,i,j,k)=basis*1.e-6_dp*real(node*(1+i+2*j+3*k),dp)
      flux(node)%y(:,j,i,k)=basis*2.e-6_dp*real(node*(2+2*i+j+k),dp)
      flux(node)%z(:,j,k,i)=basis*3.e-6_dp*real(node*(3+i+j+2*k),dp)
    end do
    end do
    end do
    associate(g=>tree%nodes(node)%geometry,s=>tree%nodes(node)%state,f=>flux(node))
    do k=1,4
    do j=1,4
    do i=1,4
      if (g%volume_fraction(i,j,k)<=0._dp) cycle
      volume=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
      s(:,i,j,k)=s(:,i,j,k)+( &
        (g%x_face_fraction(i-1,j,k)*f%x(:,i-1,j,k)-g%x_face_fraction(i,j,k)*f%x(:,i,j,k))*g%dy*g%dz+ &
        (g%y_face_fraction(i,j-1,k)*f%y(:,i,j-1,k)-g%y_face_fraction(i,j,k)*f%y(:,i,j,k))*g%dx*g%dz+ &
        (g%z_face_fraction(i,j,k-1)*f%z(:,i,j,k-1)-g%z_face_fraction(i,j,k)*f%z(:,i,j,k))*g%dx*g%dy)/volume
    end do
    end do
    end do
    end associate
  end do
  saved=tree;raw=flux
  if(rotating) then
    call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok,boundary_faces=faces,frame=frame)
  else
    call synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok,boundary_faces=faces)
  end if
  call require(ok,'flux synchronization')
  call require(any(flux(1)%x/=raw(1)%x),'nontrivial coarse replacement')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  boundary=0._dp
  associate(g=>tree%nodes(1)%geometry,f=>flux(1))
  do j=1,4
  do i=1,4
    boundary=boundary+(g%x_face_fraction(0,i,j)*f%x(:,0,i,j)-g%x_face_fraction(4,i,j)*f%x(:,4,i,j))*g%dy*g%dz &
      +(g%y_face_fraction(i,0,j)*f%y(:,i,0,j)-g%y_face_fraction(i,4,j)*f%y(:,i,4,j))*g%dx*g%dz &
      +(g%z_face_fraction(i,j,0)*f%z(:,i,j,0)-g%z_face_fraction(i,j,4)*f%z(:,i,j,4))*g%dx*g%dy
  end do
  end do
  end associate
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'composite conservation including physical boundaries')
  call check_mpi_sync()
  call check_mpi_rereflux()
  deallocate(u,t,flux)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_flux_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_held_ancestor_repair(gas_only,covered_depth,geometry_case,uniform_inventory,far_endpoint, &
      thermal_case,fine_endpoint,fine_covered,raw_endpoint,empty_register)
    logical, optional, intent(in) :: gas_only,uniform_inventory,far_endpoint,fine_endpoint,fine_covered
    logical, optional, intent(in) :: raw_endpoint,empty_register
    integer, optional, intent(in) :: covered_depth,geometry_case,thermal_case
    type(amr_patch_3d) :: patches(2)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: finest_grid,middle_grid,metrics(3)
    type(amr_tree_state_3d) :: gas,temp,held,gas_delta,held_delta,paired_delta,packed,increment
    type(amr_tree_state_3d) :: kept_gas,kept_temp,kept_held,kept_delta,kept_h_delta
    type(eb_integrated_flux_3d), allocatable :: transfers(:),scalars(:),packed_faces(:),effective(:)
    type(eb_integrated_flux_3d), allocatable :: kept_transfers(:),kept_scalars(:)
    real(dp) :: y(size(species)),vector(nv),mw,rg,cp,cv,gamma,h,e700,e1000,entropy
    real(dp) :: initial(nv+1),expected(nv+1),actual(nv+1),amount,scalar_debt
    real(dp) :: wet(nv),donor_delta(nv),pure_held,donor_h_delta,donor_volume,target_volume
    real(dp) :: trace(nv),trace_held,trace_water,lower_t,upper_t,lower_e,upper_e,energy_debt,background_e
    real(dp), allocatable :: nodal_phi(:,:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,sound
    integer :: ownership(3),layout,mode,p,d(3),i,j,k,water,nitrogen,first_level
    integer :: donor_depth,donor_root(3),donor_node,donor_leaf(3),first(3),last(3),max_mode,geom_case
    logical :: good,gas_mode,accepted_mode,uniform_mode,far_mode,fine_mode,fine_covered_mode,raw_mode,empty_mode
    integer :: target_cell(3),target_node,repair_depth,tx,ty,tz,thermal_mode
    raw_mode=.false.;empty_mode=.false.
    if(present(raw_endpoint)) raw_mode=raw_endpoint
    if(present(empty_register)) empty_mode=empty_register
    target_node=1;fine_mode=.false.
    if(present(fine_endpoint)) fine_mode=fine_endpoint
    fine_covered_mode=.false.
    if(present(fine_covered)) fine_covered_mode=fine_covered
    fine_mode=fine_mode.or.fine_covered_mode
    thermal_mode=0
    if(present(thermal_case)) thermal_mode=thermal_case
    far_mode=.false.
    if(present(far_endpoint)) far_mode=far_endpoint
    target_cell=[3,2,2];repair_depth=1
    if(far_mode) then
      target_cell=[4,3,2];repair_depth=0
    endif
    if(thermal_mode==6) target_cell=[4,2,2]
    if(fine_mode) then
      target_node=3;target_cell=[3,6,4];repair_depth=1
    endif
    if(fine_covered_mode) then
      target_node=2;target_cell=[3,3,2];repair_depth=0
    endif
    tx=target_cell(1);ty=target_cell(2);tz=target_cell(3)
    uniform_mode=.false.
    if(present(uniform_inventory)) uniform_mode=uniform_inventory
    gas_mode=.false.
    if(present(gas_only)) gas_mode=gas_only
    geom_case=0
    if(present(geometry_case)) geom_case=geometry_case
    call require(geom_case>=0.and.geom_case<=3,'ancestor geometry case')
    donor_depth=0
    if(present(covered_depth)) donor_depth=covered_depth
    call require(donor_depth>=0.and.donor_depth<=2,'ancestor donor depth')
    call require(nranks>=2,'ancestor budget owners')
    if(.not.gas_mode) call require(nranks>=3,'held ancestor budget owner and idle ranks')
    water=0;nitrogen=0
    do p=1,size(species)
      if(species(p)%name=='H2O') water=p
      if(species(p)%name=='N2') nitrogen=p
    enddo
    call require(water>0.and.nitrogen>0,'ancestor named species')
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),good)
    call require(good,'ancestor middle patch')
    call initialize_amr_patch_3d(4,8,8,1,2,1,8,1,8,2,patches(2),good)
    call require(good,'ancestor fine patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,mesh,good)
    call require(good,'ancestor three-level tree')
    if(geom_case==0) then
      call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        'x',-.1_dp,finest_grid,good)
    else if(geom_case==1) then
      call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        'y',.37_dp,finest_grid,good)
    else if(geom_case==3) then
      call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
        'y',.49875_dp,finest_grid,good)
    else
      allocate(nodal_phi(0:16,0:16,0:16))
      do k=0,16
      do j=0,16
      do i=0,16
        nodal_phi(i,j,k)=real(j,dp)/16-.05_dp*real(i,dp)/16-.33_dp
      enddo
      enddo
      enddo
      call build_level_set_eb_geometry_3d(nodal_phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],finest_grid,good)
    endif
    call require(good,'ancestor finest geometry')
    call coarsen_eb_geometry_3d(finest_grid,2,middle_grid,good)
    call require(good,'ancestor middle geometry')
    call coarsen_eb_geometry_3d(middle_grid,2,metrics(1),good)
    call require(good,'ancestor root geometry')
    call extract_eb_geometry_3d(middle_grid,[1,1,1],[4,8,8],metrics(2),good)
    call require(good,'ancestor middle extract')
    call extract_eb_geometry_3d(finest_grid,[1,1,1],[4,16,16],metrics(3),good)
    call require(good,'ancestor fine extract')
    do layout=0,1
    ownership=[0,nranks-1,nranks-1]
    if(fine_covered_mode) ownership(3)=0
    if(layout==1) ownership=0
    max_mode=2
    if(donor_depth>0) max_mode=4
    if(uniform_mode.or.thermal_mode/=0) max_mode=0
    if(fine_covered_mode) max_mode=2
    do mode=0,max_mode
      allocate(gas%fields(3),temp%fields(3),held%fields(3),gas_delta%fields(3),held_delta%fields(3), &
        paired_delta%fields(3),transfers(3),scalars(3))
      y=0;y(nitrogen)=1
      call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
      call require(good,'ancestor initial held EOS')
      pure_held=e700
      call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e1000,entropy,good)
      call require(good,'ancestor initial bulk EOS')
      vector=0;vector(1)=1;vector(5)=e1000;vector(6:)=y
      if(thermal_mode==6) then
        vector(2)=-10000;vector(5)=e1000+.5_dp*vector(2)**2
      endif
      if(thermal_mode/=0) then
        call mixture_temperature_bounds(species,lower_t,upper_t,good)
        call require(good,'thermal donor NASA temperature bounds')
        call mixture_mass_properties(species,y,lower_t,mw,rg,cp,cv,gamma,h,lower_e,entropy,good)
        call require(good,'thermal donor cold EOS bound')
        call mixture_mass_properties(species,y,upper_t,mw,rg,cp,cv,gamma,h,upper_e,entropy,good)
        call require(good,'thermal donor hot EOS bound')
      endif
      do p=1,3
        if(ownership(p)/=rank) cycle
        d=mesh%nodes(p)%shape
        allocate(gas%fields(p)%values(nv,d(1),d(2),d(3)),temp%fields(p)%values(1,d(1),d(2),d(3)), &
          held%fields(p)%values(1,d(1),d(2),d(3)),gas_delta%fields(p)%values(nv,d(1),d(2),d(3)), &
          held_delta%fields(p)%values(1,d(1),d(2),d(3)),paired_delta%fields(p)%values(nv+1,d(1),d(2),d(3)))
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          gas%fields(p)%values(:,i,j,k)=vector
        enddo
        enddo
        enddo
        temp%fields(p)%values=700;held%fields(p)%values=e700
        gas_delta%fields(p)%values=0;held_delta%fields(p)%values=0
        allocate(transfers(p)%x(nv,0:d(1),d(2),d(3)),transfers(p)%y(nv,d(1),0:d(2),d(3)), &
          transfers(p)%z(nv,d(1),d(2),0:d(3)),scalars(p)%x(1,0:d(1),d(2),d(3)), &
          scalars(p)%y(1,d(1),0:d(2),d(3)),scalars(p)%z(1,d(1),d(2),0:d(3)))
        transfers(p)%x=0;transfers(p)%y=0;transfers(p)%z=0
        scalars(p)%x=0;scalars(p)%y=0;scalars(p)%z=0
      enddo
      y=0;y(water)=.9_dp;y(nitrogen)=.1_dp
      call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
      call require(good.and.e700<0,'ancestor donor signed held energy')
      call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e1000,entropy,good)
      call require(good,'ancestor donor bulk energy')
      amount=1e-8_dp
      if(mode==1) amount=2._dp
      scalar_debt=-3e-6_dp
      if(gas_mode) scalar_debt=0
      wet=0;wet(1)=1;wet(5)=e1000;wet(6:)=y
      if(uniform_mode) then
        wet=vector;e700=pure_held
        if(.not.gas_mode) scalar_debt=-1.01_dp*pure_held
      endif
      energy_debt=0
      if(thermal_mode/=0) then
        amount=0;scalar_debt=1.1_dp*(lower_e-pure_held)
        if(thermal_mode==3) scalar_debt=-1e30_dp
        if(thermal_mode==5) scalar_debt=2.1_dp*(lower_e-pure_held)
        if(thermal_mode==4) scalar_debt=upper_e+.1_dp*(upper_e-pure_held)-pure_held
        if(thermal_mode/=2) then
          wet=vector;e700=pure_held
        endif
        if(gas_mode) then
          scalar_debt=0;energy_debt=1.1_dp*(lower_e-vector(5))
          if(thermal_mode==3) energy_debt=-1e30_dp
          if(thermal_mode==4) energy_debt=1.1_dp*(upper_e-vector(5))
          if(thermal_mode==5) energy_debt=2.1_dp*(lower_e-vector(5))
          if(thermal_mode==6) energy_debt=0
        endif
      endif
      donor_root=[3,3,2];donor_node=1;donor_leaf=donor_root
      if(donor_depth>0) then
        donor_root=[3-donor_depth,2,2];donor_node=2
        first=(donor_root-1)*2+1;last=first+1
        if(donor_depth==2) then
          donor_node=3;first=(first-1)*2+1;last=first+3
        endif
        donor_leaf=[last(1),first(2),first(3)]
        do while(metrics(donor_node)%volume_fraction(donor_leaf(1),donor_leaf(2),donor_leaf(3))<=0)
          donor_leaf(2)=donor_leaf(2)+1
          call require(donor_leaf(2)<=last(2),'ancestor donor has a physical leaf')
        enddo
      endif
      trace=vector;trace_held=pure_held;trace_water=0
      if(geom_case==2.and.donor_depth==2) then
        trace_water=1e-6_dp
        trace=vector+(trace_water/.9_dp)*(wet-vector)
        trace_held=pure_held+(trace_water/.9_dp)*(e700-pure_held)
        if(mode/=1) amount=1e-5_dp
        amount=amount+trace_water
        if(rank==ownership(1)) then
          gas%fields(1)%values(:,2,2,2)=trace;held%fields(1)%values(1,2,2,2)=trace_held
          gas%fields(1)%values(:,3,2,2)=trace;held%fields(1)%values(1,3,2,2)=trace_held
        endif
        if(rank==ownership(2)) then
          do k=3,4
          do j=3,4
          do i=3,4
            gas%fields(2)%values(:,i,j,k)=trace;held%fields(2)%values(1,i,j,k)=trace_held
          enddo
          enddo
          enddo
        endif
      endif
      if(fine_covered_mode) then
        background_e=lower_e
        if(thermal_mode==4) background_e=upper_e
        do p=1,3
          if(ownership(p)/=rank) cycle
          if(gas_mode) gas%fields(p)%values(5,:,:,:)=background_e
          if(.not.gas_mode) held%fields(p)%values=background_e
        enddo
        if(rank==ownership(target_node)) then
          gas%fields(target_node)%values(:,tx,ty,tz)=vector
          held%fields(target_node)%values(1,tx,ty,tz)=pure_held
          if(mode/=1) then
            gas%fields(2)%values(:,2,ty,tz)=wet
            held%fields(2)%values(1,2,ty,tz)=e700
            gas%fields(2)%values(:,1,ty,tz)=2*wet
            held%fields(2)%values(1,1,ty,tz)=2*e700
          endif
        endif
        if(mode/=1.and.rank==ownership(3)) then
          do k=2*tz-1,2*tz
          do j=2*ty-1,2*ty
          do i=1,4
            gas%fields(3)%values(:,i,j,k)=wet
            held%fields(3)%values(1,i,j,k)=e700
            if(i<=2) then
              gas%fields(3)%values(:,i,j,k)=2*wet
              held%fields(3)%values(1,i,j,k)=2*e700
            endif
          enddo
          enddo
          enddo
        endif
      endif
      if(rank==ownership(1)) then
        gas%fields(1)%values(:,donor_root(1),donor_root(2),donor_root(3))=wet
        held%fields(1)%values(1,donor_root(1),donor_root(2),donor_root(3))=e700
      endif
      if(rank==ownership(target_node)) then
        if(fine_mode.and..not.fine_covered_mode) then
          gas%fields(target_node)%values(:,tx-1,ty,tz)=wet
          held%fields(target_node)%values(1,tx-1,ty,tz)=e700
        endif
        gas_delta%fields(target_node)%values(5+water,tx,ty,tz)=-amount
        gas_delta%fields(target_node)%values(5+nitrogen,tx,ty,tz)=amount
        if(uniform_mode) gas_delta%fields(target_node)%values(:,tx,ty,tz)=-1.01_dp*vector
        held_delta%fields(target_node)%values(1,tx,ty,tz)=scalar_debt
        if(gas_mode.and.thermal_mode/=0) gas_delta%fields(target_node)%values(5,tx,ty,tz)=energy_debt
        if(thermal_mode==6) then
          gas%fields(target_node)%values(2,tx,ty,tz)=10000
          gas_delta%fields(target_node)%values(2,tx,ty,tz)=10000
        endif
      endif
      if(donor_depth>0) then
        if(rank==ownership(2)) then
          first=(donor_root-1)*2+1;last=first+1
          do k=first(3),last(3)
          do j=first(2),last(2)
          do i=first(1),last(1)
            gas%fields(2)%values(:,i,j,k)=wet;held%fields(2)%values(1,i,j,k)=e700
          enddo
          enddo
          enddo
        endif
        if(donor_depth==2.and.rank==ownership(3)) then
          first=((donor_root-1)*2)*2+1;last=first+3
          do k=first(3),last(3)
          do j=first(2),last(2)
          do i=first(1),last(1)
            gas%fields(3)%values(:,i,j,k)=wet;held%fields(3)%values(1,i,j,k)=e700
          enddo
          enddo
          enddo
        endif
      endif
      if(fine_mode) then
        donor_node=target_node;donor_leaf=[tx-1,ty,tz]
      endif
      donor_delta=0;donor_h_delta=0
      if(mode==3.or.(fine_mode.and.thermal_mode==2.and.mode/=1)) then
        donor_delta=.01_dp*(vector-wet)
        if(.not.gas_mode) donor_h_delta=.01_dp*(pure_held-e700)
      endif
      if(mode==4) then
        donor_delta(5+water)=-.95_dp;donor_delta(5+nitrogen)=.95_dp
      endif
      if(fine_covered_mode.and.mode==2.and.rank==ownership(3)) then
        donor_leaf=[4,2*ty,2*tz-1]
        if(gas_mode) then
          gas%fields(3)%values(5,donor_leaf(1),donor_leaf(2),donor_leaf(3))=-1e30_dp
        else
          held%fields(3)%values(1,donor_leaf(1),donor_leaf(2),donor_leaf(3))=-1e30_dp
        endif
      endif
      if(fine_covered_mode) donor_leaf=[tx-1,ty,tz]
      if(rank==ownership(donor_node)) then
        if(mode==2.and..not.fine_covered_mode) then
          if(gas_mode) then
            gas%fields(donor_node)%values(5,donor_leaf(1),donor_leaf(2),donor_leaf(3))=-1e30_dp
          else
            held%fields(donor_node)%values(1,donor_leaf(1),donor_leaf(2),donor_leaf(3))=-1e30_dp
          endif
        endif
        gas_delta%fields(donor_node)%values(:,donor_leaf(1),donor_leaf(2),donor_leaf(3))= &
          gas_delta%fields(donor_node)%values(:,donor_leaf(1),donor_leaf(2),donor_leaf(3))+donor_delta
        held_delta%fields(donor_node)%values(1,donor_leaf(1),donor_leaf(2),donor_leaf(3))= &
          held_delta%fields(donor_node)%values(1,donor_leaf(1),donor_leaf(2),donor_leaf(3))+donor_h_delta
      endif
      do p=1,3
        if(ownership(p)/=rank) cycle
        paired_delta%fields(p)%values(:nv,:,:,:)=gas_delta%fields(p)%values
        paired_delta%fields(p)%values(nv+1:nv+1,:,:,:)=held_delta%fields(p)%values
      enddo
      if(raw_mode) then
        ! Materialize the actual routed source into the raw updated family,
        ! then remove the register. Repair must fund the raw EOS deficit with
        ! either an absent or explicitly zero ledger, without recharging it.
        call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
        call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,packed,temp,packed_faces,increment,effective, &
          MPI_COMM_WORLD,good,repair_depth,2,rereflux_delta=paired_delta,rereflux_first_level=0,passive_components=1)
        call require(good,'raw thermal fixture materializes routed source once')
        do p=1,3
          if(ownership(p)/=rank) cycle
          gas%fields(p)%values=gas%fields(p)%values+increment%fields(p)%values(:nv,:,:,:)
          held%fields(p)%values=held%fields(p)%values+increment%fields(p)%values(nv+1:nv+1,:,:,:)
          gas_delta%fields(p)%values=0;held_delta%fields(p)%values=0;paired_delta%fields(p)%values=0
        enddo
      endif
      if(thermal_mode/=0.and.rank==ownership(target_node)) then
        if(gas_mode) then
          call recover_reactive_eb_cell_3d(species,metrics(target_node),target_cell, &
            gas%fields(target_node)%values(:,tx,ty,tz)+gas_delta%fields(target_node)%values(:,tx,ty,tz), &
            temp%fields(target_node)%values(1,tx,ty,tz),primitive,recovered,sound,good)
        else
          call recover_reactive_eb_cell_3d(species,metrics(target_node),target_cell,gas%fields(target_node)%values(:,tx,ty,tz), &
            temp%fields(target_node)%values(1,tx,ty,tz),primitive,recovered,sound,good, &
            internal_energy_density=held%fields(target_node)%values(1,tx,ty,tz)+ &
              held_delta%fields(target_node)%values(1,tx,ty,tz))
        endif
        call require(.not.good,'thermal-only registered endpoint is independently outside actual EOS bounds')
      endif
      kept_gas=gas;kept_temp=temp;kept_held=held
      kept_delta=gas_delta;kept_h_delta=held_delta;kept_transfers=transfers;kept_scalars=scalars
      call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
      call integrate_mpi_eb_tree_3d(mesh,ownership,metrics,packed,initial,MPI_COMM_WORLD,good)
      call require(good,'ancestor initial inventory')
      expected=initial
      target_volume=metrics(target_node)%volume_fraction(tx,ty,tz)* &
        metrics(target_node)%dx*metrics(target_node)%dy*metrics(target_node)%dz
      if(.not.raw_mode) then
      expected(5+water)=expected(5+water)-amount*target_volume
      expected(5+nitrogen)=expected(5+nitrogen)+amount*target_volume
      if(uniform_mode) expected(:nv)=initial(:nv)-1.01_dp*vector*target_volume
      expected(nv+1)=expected(nv+1)+scalar_debt*target_volume
      expected(5)=expected(5)+energy_debt*target_volume
      if(thermal_mode==6) expected(2)=expected(2)+10000*target_volume
      endif
      donor_volume=metrics(donor_node)%volume_fraction(donor_leaf(1),donor_leaf(2),donor_leaf(3))* &
        metrics(donor_node)%dx*metrics(donor_node)%dy*metrics(donor_node)%dz
      if(.not.raw_mode) expected(:nv)=expected(:nv)+donor_volume*donor_delta
      if(.not.raw_mode) expected(nv+1)=expected(nv+1)+donor_volume*donor_h_delta
      if(gas_mode) then
        if(mode==0) then
          first_level=0
          if(rank==1) first_level=1
          call repair_mpi_eb_tree_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,transfers, &
            MPI_COMM_WORLD,good,repair_depth,rereflux_delta=gas_delta,rereflux_first_level=first_level)
          call require(.not.good,'gas ancestor mixed range refusal')
          if(rank==1) then
            call repair_mpi_eb_tree_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,transfers, &
              MPI_COMM_WORLD,good,repair_depth,rereflux_delta=gas_delta)
          else
            call repair_mpi_eb_tree_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,transfers, &
              MPI_COMM_WORLD,good,repair_depth,rereflux_delta=gas_delta,rereflux_first_level=0)
          endif
          call require(.not.good,'gas ancestor mixed presence refusal')
          call repair_mpi_eb_tree_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,transfers, &
            MPI_COMM_WORLD,good,repair_depth,rereflux_delta=gas_delta,rereflux_first_level=2)
          call require(.not.good,'gas ancestor invalid range refusal')
          do p=1,3
            if(ownership(p)/=rank) cycle
            call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
              all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
              all(transfers(p)%x==kept_transfers(p)%x).and.all(transfers(p)%y==kept_transfers(p)%y).and. &
              all(transfers(p)%z==kept_transfers(p)%z),'gas ancestor metadata refusal rollback')
          enddo
        endif
        if(raw_mode.and..not.empty_mode) then
          call compare_interface_repair(mesh,ownership,metrics,gas,temp,transfers,good,minimum_level=repair_depth)
        else
          call compare_interface_repair(mesh,ownership,metrics,gas,temp,transfers,good,minimum_level=repair_depth, &
            rereflux_delta=gas_delta,rereflux_first_level=0)
        endif
      else
        if(raw_mode.and..not.empty_mode) then
          call compare_held_interface_repair(mesh,ownership,metrics,gas,temp,held,transfers,scalars,good, &
            minimum_level=repair_depth)
        else
          call compare_held_interface_repair(mesh,ownership,metrics,gas,temp,held,transfers,scalars,good, &
            gas_delta,held_delta,minimum_level=repair_depth,register_first_level=0)
        endif
      endif
      accepted_mode=mode==0.or.(donor_depth>0.and.mode==3)
      if(thermal_mode==3) accepted_mode=.false.
      if(good.neqv.accepted_mode) print *, 'Ancestor mismatch: geom/depth/layout/mode=', &
        geom_case,donor_depth,layout,mode,' thermal=',thermal_mode,' admitted=',good
      call require(good.eqv.accepted_mode,'ancestor debt admission')
      do p=1,3
        if(ownership(p)/=rank) cycle
        call require(all(gas_delta%fields(p)%values==kept_delta%fields(p)%values).and. &
          all(held_delta%fields(p)%values==kept_h_delta%fields(p)%values),'ancestor debt is read-only')
        if(accepted_mode) cycle
        call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
          all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
          all(held%fields(p)%values==kept_held%fields(p)%values),'ancestor refusal state rollback')
        call require(all(transfers(p)%x==kept_transfers(p)%x).and.all(transfers(p)%y==kept_transfers(p)%y).and. &
          all(transfers(p)%z==kept_transfers(p)%z).and.all(scalars(p)%x==kept_scalars(p)%x).and. &
          all(scalars(p)%y==kept_scalars(p)%y).and.all(scalars(p)%z==kept_scalars(p)%z), &
          'ancestor refusal both face families rollback')
      enddo
      if(accepted_mode) then
        if(thermal_mode==6.and.rank==ownership(target_node)) then
          call require(gas%fields(target_node)%values(1,tx,ty,tz)>50.and.gas%fields(target_node)%values(2,tx,ty,tz)<0, &
            'counterstreaming kinetic debt reaches the second admissible EOS interval')
          call require(transfers(1)%x(1,tx-1,ty,tz)/=0.and.transfers(1)%y(1,tx,ty,tz)/=0, &
            'counterstreaming gas debt pools cut and full-volume donors')
        endif
        if(thermal_mode==5.and.rank==ownership(target_node).and.fine_mode) then
          if(gas_mode) then
            call require(transfers(target_node)%x(5,tx-1,ty,tz)/=0.and. &
              transfers(target_node)%x(5,tx-2,ty,tz)/=0,'fine thermal gas pools distinct providers')
          else
            call require(scalars(target_node)%x(1,tx-1,ty,tz)/=0.and. &
              scalars(target_node)%x(1,tx-2,ty,tz)/=0,'fine thermal held pools distinct providers')
          endif
        endif
        if(thermal_mode==5.and.rank==ownership(target_node).and..not.fine_mode) then
          if(gas_mode) then
            call require(transfers(target_node)%x(5,tx-1,ty,tz)/=0.and.transfers(target_node)%y(5,tx,ty-1,tz)/=0, &
              'gas thermal registered endpoint pools distinct physical donors without clipping')
          else
          call require(scalars(target_node)%x(1,tx-1,ty,tz)/=0.and.scalars(target_node)%y(1,tx,ty-1,tz)/=0, &
            'thermal registered endpoint pools distinct physical donors without clipping')
          endif
        endif
        if(thermal_mode==2.and.rank==ownership(target_node)) then
          if(gas_mode) then
            call require(gas%fields(target_node)%values(5+water,tx,ty,tz)>0.and.transfers(target_node)%x(5,tx-1,ty,tz)<0, &
              'signed water donor funds gas EOS by moving every component and negative total energy together')
          else
          call require(gas%fields(target_node)%values(5+water,tx,ty,tz)>0.and.scalars(target_node)%x(1,tx-1,ty,tz)<0, &
            'signed formation-energy donor funds cold EOS by moving water and negative held energy together')
          endif
        endif
        if(fine_covered_mode.and.rank==ownership(3)) &
          call require(any(transfers(3)%x/=0),'covered fine packet debits physical descendant faces')
        if(donor_depth>0.and.rank==ownership(donor_node)) &
          call require(any(transfers(donor_node)%x/=0),'ancestor donor packet reaches physical descendant faces')
        if(rank==ownership(target_node).and..not.gas_mode) &
          call require(any(scalars(target_node)%x/=0).or.any(scalars(target_node)%y/=0), &
            'ancestor complete packet carries signed held face')
        call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
        call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,packed,temp,packed_faces,increment,effective, &
          MPI_COMM_WORLD,good,repair_depth,2,rereflux_delta=paired_delta,rereflux_first_level=0,passive_components=1)
        call require(good,'ancestor strict closure preview')
        do p=1,3
          if(ownership(p)/=rank) cycle
          packed%fields(p)%values=packed%fields(p)%values+increment%fields(p)%values
          call require(all(packed%fields(p)%values(6:nv,:,:,:)>=0),'ancestor post-debt species nonnegative')
          d=mesh%nodes(p)%shape
          do k=1,d(3)
          do j=1,d(2)
          do i=1,d(1)
            if(metrics(p)%volume_fraction(i,j,k)<=0) cycle
            if(gas_mode) then
              call recover_reactive_eb_cell_3d(species,metrics(p),[i,j,k],packed%fields(p)%values(:nv,i,j,k), &
                temp%fields(p)%values(1,i,j,k),primitive,recovered,sound,good)
            else
              call recover_reactive_eb_cell_3d(species,metrics(p),[i,j,k],packed%fields(p)%values(:nv,i,j,k), &
                temp%fields(p)%values(1,i,j,k),primitive,recovered,sound,good, &
                internal_energy_density=packed%fields(p)%values(nv+1,i,j,k))
            endif
            if(.not.good.and.uniform_mode) then
              print *, 'Uniform failed cell: node/index=',p,i,j,k
              print *, 'Uniform gas/H=',packed%fields(p)%values(:,i,j,k)
            endif
            call require(good,'ancestor post-debt actual held EOS')
          enddo
          enddo
          enddo
        enddo
        call integrate_mpi_eb_tree_3d(mesh,ownership,metrics,packed,actual,MPI_COMM_WORLD,good)
        call require(good.and.maxval(abs(actual-expected)/max(1._dp,abs(expected)))<128*epsilon(1._dp), &
          'ancestor gas and signed held inventory charge debt once')
      endif
      deallocate(gas%fields,temp%fields,held%fields,gas_delta%fields,held_delta%fields,paired_delta%fields, &
        transfers,scalars)
    enddo
    enddo
    if(rank==0) then
      if(gas_mode) then
        print *, 'PASS: serial/MPI gas ancestor debt, range contracts, idle owners and atomic refusal; depth/geom=', &
          donor_depth,geom_case
      else
        print *, 'PASS: serial/MPI ancestor debt, signed held packet, idle owners and atomic refusal; depth/geom=', &
          donor_depth,geom_case
      endif
    endif
  end subroutine check_held_ancestor_repair

  subroutine check_held_interface_repair()
    type(amr_patch_3d) :: box(1)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: fine_grid,metrics(2)
    type(amr_tree_state_3d) :: gas,temp,held,saved_gas,saved_temp,saved_held,packed,increment,expected_state
    type(amr_tree_state_3d) :: gas_delta,held_delta,paired_delta
    type(eb_integrated_flux_3d), allocatable :: transfers(:),scalars(:),kept_transfers(:),kept_scalars(:)
    type(eb_integrated_flux_3d), allocatable :: packed_faces(:),effective(:)
    real(dp) :: y(size(species)),vector(nv),mw,rg,cp,cv,gamma,h,e700,e1000,entropy
    real(dp) :: expected(nv+1),actual(nv+1),primitive(reactive_nprim(size(species))),tg,sound,scalar_change
    integer :: ownership(2),p,d(3),i,j,k,water,nitrogen
    logical :: good
    call require(nranks>=3,'paired budget needs owner and idle ranks')
    water=0;nitrogen=0
    do p=1,size(species)
      if(species(p)%name=='H2O') water=p
      if(species(p)%name=='N2') nitrogen=p
    enddo
    call require(water>0.and.nitrogen>0,'paired budget named species')
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,box(1),good)
    call require(good,'paired budget child')
    call initialize_amr_patch_tree_3d([4,4,4],[1],box,mesh,good)
    call require(good,'paired budget tree')
    call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,fine_grid,good)
    call require(good,'paired budget regular geometry')
    call coarsen_eb_geometry_3d(fine_grid,2,metrics(1),good)
    call require(good,'paired budget root geometry')
    call extract_eb_geometry_3d(fine_grid,[1,1,1],[4,8,8],metrics(2),good)
    call require(good,'paired budget child geometry')
    ownership=[0,1]
    allocate(gas%fields(2),temp%fields(2),held%fields(2),transfers(2),scalars(2))
    y=0;y(nitrogen)=1
    call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
    call require(good,'paired budget independent energy')
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e1000,entropy,good)
    call require(good,'paired budget bulk energy')
    vector=0;vector(1)=1;vector(5)=e1000;vector(6:)=y
    do p=1,2
      if(ownership(p)/=rank) cycle
      d=mesh%nodes(p)%shape
      allocate(gas%fields(p)%values(nv,d(1),d(2),d(3)),temp%fields(p)%values(1,d(1),d(2),d(3)), &
        held%fields(p)%values(1,d(1),d(2),d(3)))
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        gas%fields(p)%values(:,i,j,k)=vector
      enddo
      enddo
      enddo
      temp%fields(p)%values=700;held%fields(p)%values=e700
      allocate(transfers(p)%x(nv,0:d(1),d(2),d(3)),transfers(p)%y(nv,d(1),0:d(2),d(3)), &
        transfers(p)%z(nv,d(1),d(2),0:d(3)))
      allocate(scalars(p)%x(1,0:d(1),d(2),d(3)),scalars(p)%y(1,d(1),0:d(2),d(3)), &
        scalars(p)%z(1,d(1),d(2),0:d(3)))
      transfers(p)%x=0;transfers(p)%y=0;transfers(p)%z=0
      scalars(p)%x=0;scalars(p)%y=0;scalars(p)%z=0
    enddo
    ! One exposed coarse donor contains water, with its own signed formation
    ! energy. The interface withdraws water from an initially absent species.
    y=0;y(water)=.9_dp;y(nitrogen)=.1_dp
    call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
    call require(good,'paired water donor UEINT')
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e1000,entropy,good)
    call require(good,'paired water donor UEDEN')
    if(rank==0) then
      gas%fields(1)%values(5,3,3,2)=e1000;gas%fields(1)%values(6:,3,3,2)=y
      held%fields(1)%values(1,3,3,2)=e700
      call require(e700<0,'independent donor energy is signed')
    endif
    if(rank==1) then
      transfers(2)%x(5+water,4,3,3)=-1e-12_dp
      transfers(2)%x(5+nitrogen,4,3,3)=1e-12_dp
    endif
    saved_gas=gas;saved_temp=temp;saved_held=held;kept_transfers=transfers;kept_scalars=scalars
    call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
    call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,packed,temp,packed_faces,increment,effective, &
      MPI_COMM_WORLD,good,0,1,passive_components=1)
    call require(good,'paired flux-only preview')
    expected_state=packed
    do p=1,2
      if(ownership(p)/=rank) cycle
      expected_state%fields(p)%values=packed%fields(p)%values+increment%fields(p)%values
      if(p==1) call require(expected_state%fields(p)%values(5+water,3,2,2)<0,'actual coarse deficit represented')
    enddo
    call integrate_mpi_eb_tree_3d(mesh,ownership,metrics,expected_state,expected,MPI_COMM_WORLD,good)
    call require(good,'paired signed preview inventory')
    call compare_held_interface_repair(mesh,ownership,metrics,gas,temp,held,transfers,scalars,good)
    call require(good,'paired gas and independent energy repair')
    scalar_change=0
    do p=1,2
      if(ownership(p)/=rank) cycle
      scalar_change=max(scalar_change,maxval(abs(held%fields(p)%values-saved_held%fields(p)%values)))
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,scalar_change,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.scalar_change>0,'independent energy moves with donor packet')
    call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
    call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,packed,temp,packed_faces,increment,effective, &
      MPI_COMM_WORLD,good,0,1,passive_components=1)
    call require(good,'paired repaired preview')
    do p=1,2
      if(ownership(p)/=rank) cycle
      packed%fields(p)%values=packed%fields(p)%values+increment%fields(p)%values
      call require(all(packed%fields(p)%values(6:nv,:,:,:)>=0),'paired final species remain physical')
      d=mesh%nodes(p)%shape
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        call recover_reactive_eb_cell_3d(species,metrics(p),[i,j,k],packed%fields(p)%values(:nv,i,j,k), &
          temp%fields(p)%values(1,i,j,k),primitive,tg,sound,good, &
          internal_energy_density=packed%fields(p)%values(nv+1,i,j,k))
        call require(good.and.abs(tg-700._dp)<1e-6_dp,'paired final EOS uses actual independent energy')
      enddo
      enddo
      enddo
    enddo
    call integrate_mpi_eb_tree_3d(mesh,ownership,metrics,packed,actual,MPI_COMM_WORLD,good)
    call require(good.and.maxval(abs(actual-expected)/max(1._dp,abs(expected)))<128*epsilon(1._dp), &
      'paired repair conserves every gas and scalar inventory')
    ! A finite but thermodynamically impossible held donor cannot be hidden by
    ! its valid 1000K UEDEN. Both state and face families must remain private.
    gas=saved_gas;temp=saved_temp;held=saved_held;transfers=kept_transfers;scalars=kept_scalars
    if(rank==0) held%fields(1)%values(1,3,3,2)=1e20_dp
    saved_held=held
    call compare_held_interface_repair(mesh,ownership,metrics,gas,temp,held,transfers,scalars,good)
    call require(.not.good,'native held EOS refusal rejects paired donor')
    do p=1,2
      if(ownership(p)/=rank) cycle
      call require(all(gas%fields(p)%values==saved_gas%fields(p)%values).and. &
        all(temp%fields(p)%values==saved_temp%fields(p)%values).and. &
        all(held%fields(p)%values==saved_held%fields(p)%values),'paired EOS field rollback')
      call require(all(transfers(p)%x==kept_transfers(p)%x).and.all(transfers(p)%y==kept_transfers(p)%y).and. &
        all(transfers(p)%z==kept_transfers(p)%z).and.all(scalars(p)%x==kept_scalars(p)%x).and. &
        all(scalars(p)%y==kept_scalars(p)%y).and.all(scalars(p)%z==kept_scalars(p)%z),'paired EOS face rollback')
    enddo
    gas=saved_gas;temp=saved_temp;held=saved_held;transfers=kept_transfers;scalars=kept_scalars
    ! saved_held contains the refused scalar. Restore its actual 700K mixture,
    ! then add a separate conserved register impulse at the exposed target.
    if(rank==0) then
      y=0;y(water)=.9_dp;y(nitrogen)=.1_dp
      call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
      call require(good,'registered held donor energy')
      held%fields(1)%values(1,3,3,2)=e700
    else
      call require(.true.,'registered held donor energy')
    endif
    allocate(gas_delta%fields(2),held_delta%fields(2),paired_delta%fields(2))
    do p=1,2
      if(ownership(p)/=rank) cycle
      d=mesh%nodes(p)%shape
      allocate(gas_delta%fields(p)%values(nv,d(1),d(2),d(3)),source=0._dp)
      allocate(held_delta%fields(p)%values(1,d(1),d(2),d(3)),source=0._dp)
      allocate(paired_delta%fields(p)%values(nv+1,d(1),d(2),d(3)),source=0._dp)
    enddo
    if(rank==0) then
      gas_delta%fields(1)%values(5+water,3,2,2)=-5e-12_dp
      gas_delta%fields(1)%values(5+nitrogen,3,2,2)=5e-12_dp
      held_delta%fields(1)%values(1,3,2,2)=-1e-3_dp
    endif
    do p=1,2
      if(ownership(p)/=rank) cycle
      paired_delta%fields(p)%values(:nv,:,:,:)=gas_delta%fields(p)%values
      paired_delta%fields(p)%values(nv+1:nv+1,:,:,:)=held_delta%fields(p)%values
    enddo
    call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
    call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,packed,temp,packed_faces,increment,effective, &
      MPI_COMM_WORLD,good,0,1,rereflux_delta=paired_delta,rereflux_first_level=0,passive_components=1)
    call require(good,'paired register signed preview')
    expected_state=packed
    do p=1,2
      if(ownership(p)/=rank) cycle
      expected_state%fields(p)%values=packed%fields(p)%values+increment%fields(p)%values
    enddo
    call integrate_mpi_eb_tree_3d(mesh,ownership,metrics,expected_state,expected,MPI_COMM_WORLD,good)
    call require(good,'paired register expected composite')
    call compare_held_interface_repair(mesh,ownership,metrics,gas,temp,held,transfers,scalars,good,gas_delta,held_delta)
    call require(good,'paired register debt is conservatively funded')
    call pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
    call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,packed,temp,packed_faces,increment,effective, &
      MPI_COMM_WORLD,good,0,1,rereflux_delta=paired_delta,rereflux_first_level=0,passive_components=1)
    call require(good,'paired register final preview')
    do p=1,2
      if(ownership(p)/=rank) cycle
      packed%fields(p)%values=packed%fields(p)%values+increment%fields(p)%values
      call require(all(packed%fields(p)%values(6:nv,:,:,:)>=0),'registered final gas species nonnegative')
    enddo
    call integrate_mpi_eb_tree_3d(mesh,ownership,metrics,packed,actual,MPI_COMM_WORLD,good)
    call require(good.and.maxval(abs(actual-expected)/max(1._dp,abs(expected)))<128*epsilon(1._dp), &
      'paired registered closure charges gas and independent scalar impulse once')
    if(rank==0) print *, 'PASS: serial/MPI gas/UEINT donor, actual EOS, signed inventory, register debt and refusal'

  end subroutine check_held_interface_repair
  subroutine compare_held_interface_repair(mesh,ownership,metrics,gas,temp,held,transfers,scalars,accepted, &
      gas_delta,held_delta,minimum_level,register_first_level)
    integer, optional, intent(in) :: minimum_level,register_first_level
    type(amr_patch_tree_3d), intent(in) :: mesh
    integer, intent(in) :: ownership(:)
    type(eb_geometry_3d), intent(in) :: metrics(:)
    type(amr_tree_state_3d), intent(inout) :: gas,temp,held
    type(eb_integrated_flux_3d), intent(inout) :: transfers(:),scalars(:)
    type(amr_tree_state_3d), optional, intent(in) :: gas_delta,held_delta
    logical, intent(out) :: accepted
    type(eb_reactive_tree_3d) :: serial
    type(amr_tree_state_3d) :: serial_held,serial_delta,serial_h_delta
    type(amr_tree_state_3d) :: refused_gas,refused_temp,refused_held
    type(amr_tree_state_3d) :: original_gas,original_temp,original_held
    type(eb_integrated_flux_3d), allocatable :: serial_faces(:),serial_scalars(:)
    type(eb_integrated_flux_3d), allocatable :: refused_faces(:),refused_scalars(:)
    type(eb_integrated_flux_3d), allocatable :: original_faces(:),original_scalars(:)
    integer :: p,d(3),mpi_status,depth,first_depth,minimum_calls,maximum_calls
    logical :: serial_ok,refused_ok
    depth=0;first_depth=0
    if(present(minimum_level)) depth=minimum_level
    if(present(register_first_level)) first_depth=register_first_level
    call require(present(gas_delta).eqv.present(held_delta),'paired register fixture presence')
    if(present(gas_delta)) allocate(serial_delta%fields(size(ownership)),serial_h_delta%fields(size(ownership)))
    serial%topology=mesh
    allocate(serial%nodes(size(ownership)),serial_held%fields(size(ownership)), &
      serial_faces(size(ownership)),serial_scalars(size(ownership)))
    do p=1,size(ownership)
      d=mesh%nodes(p)%shape;serial%nodes(p)%geometry=metrics(p)
      allocate(serial%nodes(p)%state(nv,d(1),d(2),d(3)),serial%nodes(p)%temperature(d(1),d(2),d(3)))
      allocate(serial_held%fields(p)%values(1,d(1),d(2),d(3)))
      allocate(serial_faces(p)%x(nv,0:d(1),d(2),d(3)),serial_faces(p)%y(nv,d(1),0:d(2),d(3)), &
        serial_faces(p)%z(nv,d(1),d(2),0:d(3)))
      allocate(serial_scalars(p)%x(1,0:d(1),d(2),d(3)),serial_scalars(p)%y(1,d(1),0:d(2),d(3)), &
        serial_scalars(p)%z(1,d(1),d(2),0:d(3)))
      if(rank==ownership(p)) then
        serial%nodes(p)%state=gas%fields(p)%values
        serial%nodes(p)%temperature=temp%fields(p)%values(1,:,:,:)
        serial_held%fields(p)%values=held%fields(p)%values
        serial_faces(p)=transfers(p);serial_scalars(p)=scalars(p)
      endif
      call MPI_Bcast(serial%nodes(p)%state,size(serial%nodes(p)%state),MPI_DOUBLE_PRECISION, &
        ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial gas broadcast')
      call MPI_Bcast(serial%nodes(p)%temperature,size(serial%nodes(p)%temperature),MPI_DOUBLE_PRECISION, &
        ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial temperature broadcast')
      call MPI_Bcast(serial_held%fields(p)%values,size(serial_held%fields(p)%values),MPI_DOUBLE_PRECISION, &
        ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial independent energy broadcast')
      call MPI_Bcast(serial_faces(p)%x,size(serial_faces(p)%x),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial gas x broadcast')
      call MPI_Bcast(serial_faces(p)%y,size(serial_faces(p)%y),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial gas y broadcast')
      call MPI_Bcast(serial_faces(p)%z,size(serial_faces(p)%z),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial gas z broadcast')
      call MPI_Bcast(serial_scalars(p)%x,size(serial_scalars(p)%x),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial scalar x broadcast')
      call MPI_Bcast(serial_scalars(p)%y,size(serial_scalars(p)%y),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial scalar y broadcast')
      call MPI_Bcast(serial_scalars(p)%z,size(serial_scalars(p)%z),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'held serial scalar z broadcast')
      if(present(gas_delta)) then
        allocate(serial_delta%fields(p)%values(nv,d(1),d(2),d(3)),serial_h_delta%fields(p)%values(1,d(1),d(2),d(3)))
        if(rank==ownership(p)) then
          serial_delta%fields(p)%values=gas_delta%fields(p)%values
          serial_h_delta%fields(p)%values=held_delta%fields(p)%values
        endif
        call MPI_Bcast(serial_delta%fields(p)%values,size(serial_delta%fields(p)%values),MPI_DOUBLE_PRECISION, &
          ownership(p),MPI_COMM_WORLD,mpi_status)
        call require(mpi_status==MPI_SUCCESS,'held serial gas register broadcast')
        call MPI_Bcast(serial_h_delta%fields(p)%values,size(serial_h_delta%fields(p)%values),MPI_DOUBLE_PRECISION, &
          ownership(p),MPI_COMM_WORLD,mpi_status)
        call require(mpi_status==MPI_SUCCESS,'held serial scalar register broadcast')
      endif
    enddo
    if(present(gas_delta)) then
      if(native_query_mode) then
        refused_gas=gas;refused_temp=temp;refused_held=held;refused_faces=transfers;refused_scalars=scalars
        original_gas=gas;original_temp=temp;original_held=held;original_faces=transfers;original_scalars=scalars
        call repair_reactive_eb_tree_held_interface_budget_3d(species,serial,serial_held,serial_faces,serial_scalars, &
          serial_ok,depth,rereflux_delta=serial_delta,internal_energy_rereflux_delta=serial_h_delta, &
          rereflux_first_level=first_depth,native_projection=serial_native_query)
        native_query_calls=0
        call repair_mpi_eb_tree_held_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,held,transfers,scalars, &
          MPI_COMM_WORLD,accepted,depth,rereflux_delta=gas_delta,internal_energy_rereflux_delta=held_delta, &
          rereflux_first_level=first_depth,native_projection=mpi_native_query)
        call MPI_Allreduce(native_query_calls,minimum_calls,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_status)
        call require(mpi_status==MPI_SUCCESS,'native query minimum count')
        call MPI_Allreduce(native_query_calls,maximum_calls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_status)
        call require(mpi_status==MPI_SUCCESS.and.minimum_calls>0.and.minimum_calls==maximum_calls, &
          'owner and idle ranks enter every actual native packet query together')
      else
      call repair_reactive_eb_tree_held_interface_budget_3d(species,serial,serial_held,serial_faces,serial_scalars, &
        serial_ok,depth,rereflux_delta=serial_delta,internal_energy_rereflux_delta=serial_h_delta, &
        rereflux_first_level=first_depth)
      call repair_mpi_eb_tree_held_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,held,transfers,scalars, &
        MPI_COMM_WORLD,accepted,depth,rereflux_delta=gas_delta,internal_energy_rereflux_delta=held_delta, &
        rereflux_first_level=first_depth)
      endif
    else
    call repair_reactive_eb_tree_held_interface_budget_3d(species,serial,serial_held,serial_faces,serial_scalars,serial_ok,depth)
    call repair_mpi_eb_tree_held_interface_budget_3d(species,mesh,ownership,metrics,gas,temp,held,transfers,scalars, &
      MPI_COMM_WORLD,accepted,depth)
    endif
    call require(accepted.eqv.serial_ok,'serial/MPI paired donor admission agrees')
    do p=1,size(ownership)
      if(rank/=ownership(p)) cycle
      call require(all(gas%fields(p)%values==serial%nodes(p)%state), &
        'serial/MPI paired donor gas agrees exactly')
      call require(all(temp%fields(p)%values(1,:,:,:)==serial%nodes(p)%temperature), &
        'serial/MPI paired donor actual energy EOS agrees exactly')
      call require(all(held%fields(p)%values==serial_held%fields(p)%values), &
        'serial/MPI paired donor independent energy agrees exactly')
      call require(all(transfers(p)%x==serial_faces(p)%x).and.all(transfers(p)%y==serial_faces(p)%y).and. &
        all(transfers(p)%z==serial_faces(p)%z),'serial/MPI paired donor gas faces agree exactly')
      call require(all(scalars(p)%x==serial_scalars(p)%x).and.all(scalars(p)%y==serial_scalars(p)%y).and. &
        all(scalars(p)%z==serial_scalars(p)%z),'serial/MPI paired donor independent energy faces agree exactly')
    enddo
    if(native_query_mode) then
      call require(accepted,'native query fixture is admitted before refusal')
      reject_idle_query=.true.;native_query_calls=0
      call repair_mpi_eb_tree_held_interface_budget_3d(species,mesh,ownership,metrics,refused_gas,refused_temp,refused_held, &
        refused_faces,refused_scalars,MPI_COMM_WORLD,refused_ok,depth,rereflux_delta=gas_delta, &
        internal_energy_rereflux_delta=held_delta,rereflux_first_level=first_depth,native_projection=mpi_native_query)
      reject_idle_query=.false.
      call MPI_Allreduce(native_query_calls,minimum_calls,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS,'refused native query minimum count')
      call MPI_Allreduce(native_query_calls,maximum_calls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_status)
      call require(mpi_status==MPI_SUCCESS.and..not.refused_ok.and.minimum_calls>0.and.minimum_calls==maximum_calls, &
        'idle-rank native query refusal is collective with no packet publication')
      do p=1,size(ownership)
        if(rank/=ownership(p)) cycle
        call require(all(refused_gas%fields(p)%values==original_gas%fields(p)%values).and. &
          all(refused_temp%fields(p)%values==original_temp%fields(p)%values).and. &
          all(refused_held%fields(p)%values==original_held%fields(p)%values), &
          'collective idle refusal retains gas/T/UEINT exactly')
        call require(all(refused_faces(p)%x==original_faces(p)%x).and. &
          all(refused_faces(p)%y==original_faces(p)%y).and.all(refused_faces(p)%z==original_faces(p)%z).and. &
          all(refused_scalars(p)%x==original_scalars(p)%x).and.all(refused_scalars(p)%y==original_scalars(p)%y).and. &
          all(refused_scalars(p)%z==original_scalars(p)%z),'collective idle refusal retains both face families exactly')
      enddo
    endif
  end subroutine compare_held_interface_repair

  subroutine serial_native_query(node,cell,pre,require_post,accepted)
    integer, intent(in) :: node,cell(3)
    real(dp), intent(in) :: pre(:)
    logical, intent(in) :: require_post
    logical, intent(out) :: accepted
    accepted=node>0.and.all(cell>0).and.all(ieee_is_finite(pre))
  end subroutine serial_native_query

  subroutine mpi_native_query(node,cell,pre,require_post,accepted)
    integer, intent(in) :: node,cell(3)
    real(dp), intent(in) :: pre(:)
    logical, intent(in) :: require_post
    logical, intent(out) :: accepted
    logical :: local_ok
    integer :: ierr
    call serial_native_query(node,cell,pre,require_post,local_ok)
    local_ok=local_ok.and..not.(reject_idle_query.and.rank==1)
    native_query_calls=native_query_calls+1
    call MPI_Allreduce(local_ok,accepted,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) accepted=.false.
  end subroutine mpi_native_query

  subroutine pack_held_budget_fixture(mesh,ownership,gas,held,transfers,scalars,packed,packed_faces)
    type(amr_patch_tree_3d), intent(in) :: mesh
    integer, intent(in) :: ownership(:)
    type(amr_tree_state_3d), intent(in) :: gas,held
    type(eb_integrated_flux_3d), intent(in) :: transfers(:),scalars(:)
    type(amr_tree_state_3d), intent(out) :: packed
    type(eb_integrated_flux_3d), allocatable, intent(out) :: packed_faces(:)
    integer :: p,d(3)
    allocate(packed%fields(size(ownership)),packed_faces(size(ownership)))
    do p=1,size(ownership)
      if(ownership(p)/=rank) cycle
      d=mesh%nodes(p)%shape
      allocate(packed%fields(p)%values(nv+1,d(1),d(2),d(3)))
      packed%fields(p)%values(:nv,:,:,:)=gas%fields(p)%values
      packed%fields(p)%values(nv+1:nv+1,:,:,:)=held%fields(p)%values
      allocate(packed_faces(p)%x(nv+1,0:d(1),d(2),d(3)),packed_faces(p)%y(nv+1,d(1),0:d(2),d(3)), &
        packed_faces(p)%z(nv+1,d(1),d(2),0:d(3)))
      packed_faces(p)%x(:nv,:,:,:)=transfers(p)%x;packed_faces(p)%x(nv+1:nv+1,:,:,:)=scalars(p)%x
      packed_faces(p)%y(:nv,:,:,:)=transfers(p)%y;packed_faces(p)%y(nv+1:nv+1,:,:,:)=scalars(p)%y
      packed_faces(p)%z(:nv,:,:,:)=transfers(p)%z;packed_faces(p)%z(nv+1:nv+1,:,:,:)=scalars(p)%z
    enddo
  end subroutine pack_held_budget_fixture

  subroutine check_ancestor_register_close()
    type(amr_patch_3d) :: patches(2)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: fg,mg,cg,padded,physical(3)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(eb_reactive_tree_3d) :: seed,expected
    type(amr_tree_state_3d) :: correction,spill,gas,temp,kept_gas,kept_temp
    type(state_redist_register_3d) :: baseline(3),serial_regs(3),regs(3),kept_regs(3)
    type(eb_integrated_flux_3d) :: serial_flux(3),f(3),kept_flux(3)
    real(dp) :: y(size(species)),pure(nv),wet(nv),mw,rg,cp,cv,gamma,h,e,entropy
    real(dp) :: cs(nv,4,4,4),ct(4,4,4),before(nv),after(nv),impulse(nv)
    integer :: p,i,j,k,layout,water,nitrogen,ownership(3)
    logical :: good,has_spill
    water=0;nitrogen=0
    do p=1,size(species)
      if(species(p)%name=='H2O') water=p
      if(species(p)%name=='N2') nitrogen=p
    enddo
    call require(water>0.and.nitrogen>0,'cut ancestor named species')
    y=0;y(nitrogen)=1
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e,entropy,good)
    call require(good,'cut ancestor pure EOS')
    pure=0;pure(1)=1;pure(5)=e;pure(6:)=y
    y=0;y(water)=.9_dp;y(nitrogen)=.1_dp
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e,entropy,good)
    call require(good,'cut ancestor wet EOS')
    wet=0;wet(1)=1;wet(5)=e;wet(6:)=y
    do k=1,4
    do j=1,4
    do i=1,4
      cs(:,i,j,k)=pure
    enddo
    enddo
    enddo
    ct=1000
    call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',.37_dp,fg,good)
    call require(good,'cut ancestor finest geometry')
    call coarsen_eb_geometry_3d(fg,2,mg,good)
    call require(good,'cut ancestor middle geometry')
    call coarsen_eb_geometry_3d(mg,2,cg,good)
    call require(good,'cut ancestor coarse geometry')
    call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patches(1),good)
    call require(good,'cut ancestor middle patch')
    call initialize_amr_patch_3d(4,4,4,2,3,1,2,1,2,2,patches(2),good)
    call require(good,'cut ancestor fine patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,mesh,good)
    call require(good,'cut ancestor topology')
    call initialize_reactive_eb_tree_3d(species,mesh,fg,cs,ct,seed,good)
    call require(good,'cut ancestor seed')
    do p=1,3
      physical(p)=seed%nodes(p)%geometry
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical(p)%volume_fraction(i,j,k)==1) seed%nodes(p)%state(:,i,j,k)=wet
      enddo
      enddo
      enddo
      call initialize_state_redist_tree_register_3d(mesh,p,nv,6,baseline(p),good,scale=.5_dp)
      call require(good,'cut ancestor register masks')
      allocate(serial_flux(p)%x(nv,0:4,4,4),serial_flux(p)%y(nv,4,0:4,4),serial_flux(p)%z(nv,4,4,0:4))
      serial_flux(p)%x=0;serial_flux(p)%y=0;serial_flux(p)%z=0
    enddo
    call average_down_reactive_eb_tree_3d(species,seed,good)
    call require(good,'cut ancestor seed restriction')
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'x',.37_dp,padded,good)
    call require(good,'cut ancestor padded geometry')
    call build_eb_tree_extended_geometry_3d(mesh,physical,padded,wide,good)
    call require(good,'cut ancestor extended geometry')
    impulse=0;impulse(5+water)=-1e-8_dp;impulse(5+nitrogen)=1e-8_dp
    baseline(2)%coarse_delta(:,7,7,7)=impulse
    baseline(2)%coarse_delta(:,10,7,7)=-physical(2)%volume_fraction(1,1,1)*impulse
    ! Keep an unrelated older clock ledger live; depth-one closure must not consume it.
    baseline(1)%coarse_delta(:,9,7,8)=1e-7_dp*pure
    call assemble_state_redist_interface_rereflux_3d(mesh,wide,baseline,nv,1,correction,good,spill=spill, &
      spill_active=has_spill)
    call require(good.and.has_spill,'cut ancestor real register spill')
    call require(any(spill%fields(1)%values(5+water,:,:,:)<0),'cut ancestor register has negative ancestor water')
    call integrate_reactive_eb_tree_3d(seed,before,good)
    call require(good,'cut ancestor initial inventory')
    expected=seed;serial_regs=baseline
    call synchronize_reactive_eb_tree_register_flux_3d(species,expected,serial_flux,wide,serial_regs,1,good)
    call require(good,'cut ancestor serial registered deficit closure')
    call integrate_reactive_eb_tree_3d(expected,after,good)
    call require(good.and.maxval(abs(after-before)/max(1._dp,abs(before)))<128*epsilon(1._dp), &
      'cut ancestor opposite registration packets conserve every component')
    call require(any(serial_flux(1)%x/=0).or.any(serial_flux(1)%y/=0).or.any(serial_flux(1)%z/=0), &
      'cut ancestor registered debt funds an ancestor face')
    do layout=0,2
      ownership=[0,nranks-1,nranks-1]
      if(layout==1) ownership=nranks-1-ownership
      if(layout==2) ownership=0
      allocate(gas%fields(3),temp%fields(3))
      block
        type(state_redist_register_3d) :: empty
        do p=1,3
          regs(p)=empty
          if(ownership(p)/=rank) cycle
          regs(p)=baseline(p)
          gas%fields(p)%values=seed%nodes(p)%state
          allocate(temp%fields(p)%values(1,4,4,4))
          temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          allocate(f(p)%x(nv,0:4,4,4),f(p)%y(nv,4,0:4,4),f(p)%z(nv,4,4,0:4))
          f(p)%x=0;f(p)%y=0;f(p)%z=0
        enddo
      end block
      call synchronize_mpi_eb_tree_register_flux_3d(species,mesh,ownership,physical,wide,gas,temp,f,regs, &
        1,MPI_COMM_WORLD,good)
      call require(good,'cut ancestor MPI registered deficit closure')
      do p=1,3
        if(ownership(p)/=rank) cycle
        call require(all(gas%fields(p)%values==expected%nodes(p)%state).and. &
          all(temp%fields(p)%values(1,:,:,:)==expected%nodes(p)%temperature), &
          'cut ancestor exact registered gas and temperature')
        call require(all(f(p)%x==serial_flux(p)%x).and.all(f(p)%y==serial_flux(p)%y).and. &
          all(f(p)%z==serial_flux(p)%z),'cut ancestor exact registered faces')
        if(allocated(regs(p)%coarse_delta)) call require(all(regs(p)%coarse_delta==serial_regs(p)%coarse_delta), &
          'cut ancestor consumed and older coarse ledger')
        if(allocated(regs(p)%fine_delta)) call require(all(regs(p)%fine_delta==serial_regs(p)%fine_delta), &
          'cut ancestor consumed and older fine ledger')
        deallocate(f(p)%x,f(p)%y,f(p)%z)
      enddo
      deallocate(gas%fields,temp%fields)
    enddo
    if(rank==0) print *, 'PASS: cut-cell StateRedist ancestor spill, deficit closure, exact owners and older ledgers'
  end subroutine check_ancestor_register_close

  subroutine check_register_close()
    type(amr_patch_3d) :: patches(2)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: fg,mg,cg,padded,physical(3)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(reactive_eb_3d_config) :: cc
    type(eb_reactive_tree_3d) :: seed,expected
    type(amr_tree_state_3d) :: local,lt,kept,kt
    type(amr_tree_state_3d) :: scalar,serial_scalar,kept_scalar
    type(state_redist_register_3d) :: regs(3),initial_regs(3),baseline_regs(3),serial_regs(3),kept_regs(3)
    type(state_redist_register_3d) :: scalar_regs(3),scalar_serial_regs(3),scalar_kept_regs(3)
    type(eb_integrated_flux_3d) :: f(3),serial_flux(3),kept_flux(3)
    type(eb_integrated_flux_3d) :: scalar_flux(3),scalar_serial_flux(3),scalar_kept_flux(3)
    real(dp) :: cs(nv,4,4,4),ct(4,4,4),density(nv),before(nv),after(nv)
    integer :: p,mode,scenario,ownership(3)
    logical :: good
    cc%nx=4;cc%ny=4;cc%nz=4
    call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,fg,good)
    call require(good,'register close fine geometry')
    call coarsen_eb_geometry_3d(fg,2,mg,good)
    call require(good,'register close middle geometry')
    call coarsen_eb_geometry_3d(mg,2,cg,good)
    call require(good,'register close coarse geometry')
    call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patches(1),good)
    call require(good,'register close child')
    call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(2),good)
    call require(good,'register close descendant')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,mesh,good)
    call require(good,'register close mesh')
    call initialize_reactive_eb_density_sheet_3d(species,cc,cg,cs,ct,good)
    call require(good,'register close initial gas')
    call initialize_reactive_eb_tree_3d(species,mesh,fg,cs,ct,seed,good)
    call require(good,'register close initial tree')
    do p=1,3
      physical(p)=seed%nodes(p)%geometry
      call initialize_state_redist_tree_register_3d(mesh,p,nv,6,initial_regs(p),good,scale=.5_dp)
      call require(good,'register close initial masks')
      allocate(serial_flux(p)%x(nv,0:4,4,4),serial_flux(p)%y(nv,4,0:4,4),serial_flux(p)%z(nv,4,4,0:4))
      serial_flux(p)%x=0;serial_flux(p)%y=0;serial_flux(p)%z=0
    enddo
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'x',.37_dp,padded,good)
    call require(good,'register close padded geometry')
    call build_eb_tree_extended_geometry_3d(mesh,physical,padded,wide,good)
    call require(good,'register close extended geometry')
    density=1.e-7_dp*seed%nodes(1)%state(:,2,1,2)
    initial_regs(1)%coarse_delta(:,8,7,8)=density
    initial_regs(2)%fine_delta(:,8,6,7)=.2_dp*density
    initial_regs(2)%coarse_delta(:,9,7,7)=.3_dp*density
    initial_regs(3)%fine_delta(:,6,7,7)=.4_dp*density
    baseline_regs=initial_regs
    do scenario=2,1,-1
      initial_regs=baseline_regs
      do p=1,3
        serial_flux(p)%x=0;serial_flux(p)%y=0;serial_flux(p)%z=0
      end do
      if(scenario==2) then
        ! Equal and opposite whole-state packets at two exposed coarse cells.
        ! The first is overdrawn before reflux; the fine interface has ample
        ! inventory. No mass/species/energy clipping can satisfy this control.
        density=1.01_dp*seed%nodes(1)%state(:,3,1,2)
        initial_regs(1)%coarse_delta=0
        initial_regs(2)%fine_delta=0
        initial_regs(1)%coarse_delta(:,9,7,8)=-density
        initial_regs(1)%coarse_delta(:,9,10,8)=density
      end if
    expected=seed;serial_regs=initial_regs
    call synchronize_reactive_eb_tree_register_flux_3d(species,expected,serial_flux,wide,serial_regs,0,good)
    if(rank==0) print *, 'Register close scenario:',scenario,' serial accepted:',good
    if(.not.good) then
      block
        type(amr_tree_state_3d) :: diagnostic,spill,increments,source
        type(eb_integrated_flux_3d), allocatable :: replacements(:)
        real(dp) :: qp(reactive_nprim(size(species))),tg,sound,value(nv),budget(3,nv)
        real(dp), allocatable :: patch_increment(:,:,:,:)
        logical :: valid,has_spill,physical_ok
        integer :: node,ii,jj,kk
        call assemble_state_redist_interface_rereflux_3d(mesh,wide,initial_regs,nv,0,diagnostic,valid, &
          spill=spill,spill_active=has_spill)
        if(rank==0) print *, 'Register source assembly:',valid,' ancestor spill:',has_spill
        if(.not.valid) then
          call assemble_state_redist_interface_source_3d(mesh,wide,initial_regs,nv,0,source,valid)
          if(rank==0) print *, 'Register source restriction:',valid
          if(valid) then
            call assemble_state_redist_patch_rereflux_3d(wide(1),initial_regs(1),[4,4,4], &
              source%fields(1)%values,patch_increment,valid)
            if(rank==0) print *, 'Register cut redistribution:',valid
            if(valid) then
              call state_redist_patch_rereflux_budget_3d(wide(1),[4,4,4],source%fields(1)%values, &
                patch_increment,budget,valid)
              if(rank==0) print *, 'Register source/recipient density budget:',budget(:,1)
            endif
          endif
          valid=.false.
        endif
        if(valid) then
          call assemble_reactive_eb_tree_reflux_3d(species,seed,serial_flux,increments,replacements,valid, &
            0,1,rereflux_delta=diagnostic,rereflux_first_level=0)
          if(rank==0) print *, 'Register endpoint assembly:',valid
          if(valid) then
            do node=1,3
            do kk=1,4
            do jj=1,4
            do ii=1,4
              if(physical(node)%volume_fraction(ii,jj,kk)<=0) cycle
              value=seed%nodes(node)%state(:,ii,jj,kk)+increments%fields(node)%values(:,ii,jj,kk)
              call recover_reactive_eb_cell_3d(species,physical(node),[ii,jj,kk],value, &
                seed%nodes(node)%temperature(ii,jj,kk),qp,tg,sound,physical_ok)
              if(rank==0.and.(.not.physical_ok.or.any(value(6:)<0))) print *, &
                'Register endpoint refused:',node,ii,jj,kk,value(1),minval(value(6:)),physical_ok
            enddo
            enddo
            enddo
            enddo
          endif
        endif
      end block
    endif
    call require(good,'register close serial accepted control')
    if(scenario==2) then
      call require(any(serial_flux(2)%y/=0),'register deficit changes a conserved interface transfer')
      call integrate_reactive_eb_tree_3d(seed,before,good)
      call require(good,'register deficit initial integral')
      call integrate_reactive_eb_tree_3d(expected,after,good)
      call require(good,'register deficit accepted integral')
      call require(all(abs(after-before)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
        'register deficit conserves every composite component')
    end if
    do mode=1,3
      ownership=[0,1,mod(2,nranks)]
      if(mode==2) ownership=nranks-1-ownership
      if(mode==3) ownership=0
      allocate(local%fields(3),lt%fields(3))
      block
        type(state_redist_register_3d) :: empty
        do p=1,3
          regs(p)=empty
          if(ownership(p)/=rank) cycle
          regs(p)=initial_regs(p)
          local%fields(p)%values=seed%nodes(p)%state
          allocate(lt%fields(p)%values(1,4,4,4))
          lt%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          allocate(f(p)%x(nv,0:4,4,4),f(p)%y(nv,4,0:4,4),f(p)%z(nv,4,4,0:4))
          f(p)%x=0;f(p)%y=0;f(p)%z=0
        enddo
      end block
      call synchronize_mpi_eb_tree_register_flux_3d(species,mesh,ownership,physical,wide,local,lt,f,regs, &
        0,MPI_COMM_WORLD,good)
      call require(good,'register close MPI acceptance')
      do p=1,3
        if(ownership(p)/=rank) cycle
        call require(all(local%fields(p)%values==expected%nodes(p)%state).and. &
          all(lt%fields(p)%values(1,:,:,:)==expected%nodes(p)%temperature),'register close exact gas comparison')
        call require(all(f(p)%x==serial_flux(p)%x).and.all(f(p)%y==serial_flux(p)%y).and. &
          all(f(p)%z==serial_flux(p)%z),'register close exact flux comparison')
        if(allocated(regs(p)%coarse_delta)) then
          call require(all(regs(p)%coarse_delta==serial_regs(p)%coarse_delta),'consumed coarse and retained next ledger')
        endif
        if(allocated(regs(p)%fine_delta)) then
          call require(all(regs(p)%fine_delta==serial_regs(p)%fine_delta),'consumed fine and retained deeper ledger')
        endif
      enddo
      kept=local;kt=lt;kept_flux=f
      call integrate_reactive_eb_tree_3d(expected,before,good)
      call require(good,'register close inventory for impossible debt')
      ! After a successful repair this receiver can be nearly empty. A multiple
      ! of its local state is still fundable; remove twice the whole composite
      ! inventory instead to exercise a provably impossible conserved debt.
      if(rank==ownership(1)) regs(1)%coarse_delta(:,9,7,8)=-2._dp*before/ &
        (physical(1)%volume_fraction(3,1,2)*physical(1)%dx*physical(1)%dy*physical(1)%dz)
      kept_regs=regs
      call synchronize_mpi_eb_tree_register_flux_3d(species,mesh,ownership,physical,wide,local,lt,f,regs, &
        0,MPI_COMM_WORLD,good)
      call require(.not.good,'register close nonphysical gas collectively rejected')
      do p=1,3
        if(ownership(p)/=rank) cycle
        call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
          all(lt%fields(p)%values==kt%fields(p)%values),'register close gas rollback')
        call require(all(f(p)%x==kept_flux(p)%x).and.all(f(p)%y==kept_flux(p)%y).and. &
          all(f(p)%z==kept_flux(p)%z),'register close flux rollback')
        if(allocated(regs(p)%coarse_delta)) then
          call require(all(regs(p)%coarse_delta==kept_regs(p)%coarse_delta),'register close coarse rollback')
        endif
        if(allocated(regs(p)%fine_delta)) then
          call require(all(regs(p)%fine_delta==kept_regs(p)%fine_delta),'register close fine rollback')
        endif
        deallocate(f(p)%x,f(p)%y,f(p)%z)
      enddo
      deallocate(local%fields,lt%fields)
    enddo
    end do
    allocate(serial_scalar%fields(3))
    do p=1,3
      serial_scalar%fields(p)%values=seed%nodes(p)%state(5:5,:,:,:)
      call initialize_state_redist_tree_register_3d(mesh,p,1,6,scalar_serial_regs(p),good,scale=.5_dp)
      call require(good,'scalar register masks')
      allocate(scalar_serial_flux(p)%x(1,0:4,4,4),scalar_serial_flux(p)%y(1,4,0:4,4), &
        scalar_serial_flux(p)%z(1,4,4,0:4))
      scalar_serial_flux(p)%x=0;scalar_serial_flux(p)%y=0;scalar_serial_flux(p)%z=0
    enddo
    scalar_serial_regs(1)%coarse_delta(1,8,7,8)=1.e-7_dp*seed%nodes(1)%state(5,2,1,2)
    scalar_serial_regs(2)%fine_delta(1,8,6,7)=.2e-7_dp*seed%nodes(1)%state(5,2,1,2)
    scalar_serial_regs(2)%coarse_delta(1,9,7,7)=.3e-7_dp*seed%nodes(1)%state(5,2,1,2)
    scalar_serial_regs(3)%fine_delta(1,6,7,7)=.4e-7_dp*seed%nodes(1)%state(5,2,1,2)
    do mode=1,3
      ownership=[0,1,mod(2,nranks)]
      if(mode==2) ownership=nranks-1-ownership
      if(mode==3) ownership=0
      scalar_regs=scalar_serial_regs
      allocate(scalar%fields(3))
      do p=1,3
        if(ownership(p)/=rank) then
          block
            type(state_redist_register_3d) :: empty
            scalar_regs(p)=empty
          end block
          cycle
        endif
        scalar%fields(p)%values=seed%nodes(p)%state(5:5,:,:,:)
        allocate(scalar_flux(p)%x(1,0:4,4,4),scalar_flux(p)%y(1,4,0:4,4),scalar_flux(p)%z(1,4,4,0:4))
        scalar_flux(p)%x=0;scalar_flux(p)%y=0;scalar_flux(p)%z=0
      enddo
      ! The serial control is reset for each ownership arrangement.
      do p=1,3
        scalar_serial_flux(p)%x=0;scalar_serial_flux(p)%y=0;scalar_serial_flux(p)%z=0
      enddo
      block
        type(amr_tree_state_3d) :: expected_scalar
        type(state_redist_register_3d) :: expected_regs(3)
        expected_scalar=serial_scalar;expected_regs=scalar_serial_regs
        call synchronize_reactive_eb_tree_scalar_register_flux_3d(species,seed,expected_scalar, &
          scalar_serial_flux,wide,expected_regs,0,good)
        call require(good,'scalar register serial acceptance')
        call synchronize_mpi_eb_tree_scalar_register_flux_3d(species,mesh,ownership,physical,wide,scalar, &
          scalar_flux,scalar_regs,0,MPI_COMM_WORLD,good)
        call require(good,'scalar register MPI acceptance')
        do p=1,3
          if(ownership(p)/=rank) cycle
          call require(all(scalar%fields(p)%values==expected_scalar%fields(p)%values), &
            'scalar register serial/MPI field equality')
          call require(all(scalar_flux(p)%x==scalar_serial_flux(p)%x).and. &
            all(scalar_flux(p)%y==scalar_serial_flux(p)%y).and. &
            all(scalar_flux(p)%z==scalar_serial_flux(p)%z),'scalar register serial/MPI flux equality')
          if(allocated(scalar_regs(p)%coarse_delta)) call require( &
            all(scalar_regs(p)%coarse_delta==expected_regs(p)%coarse_delta),'scalar coarse ledger consumption')
          if(allocated(scalar_regs(p)%fine_delta)) call require( &
            all(scalar_regs(p)%fine_delta==expected_regs(p)%fine_delta),'scalar fine ledger consumption')
        enddo
      end block
      kept_scalar=scalar;scalar_kept_flux=scalar_flux;scalar_kept_regs=scalar_regs
      if(rank==ownership(1)) scalar_regs(1)%coarse_delta(1,8,7,8)=ieee_value(0._dp,ieee_quiet_nan)
      call synchronize_mpi_eb_tree_scalar_register_flux_3d(species,mesh,ownership,physical,wide,scalar, &
        scalar_flux,scalar_regs,0,MPI_COMM_WORLD,good)
      call require(.not.good,'scalar register invalid packet collectively rejected')
      do p=1,3
        if(ownership(p)/=rank) cycle
        call require(all(scalar%fields(p)%values==kept_scalar%fields(p)%values), &
          'scalar register field rollback')
        call require(all(scalar_flux(p)%x==scalar_kept_flux(p)%x).and. &
          all(scalar_flux(p)%y==scalar_kept_flux(p)%y).and. &
          all(scalar_flux(p)%z==scalar_kept_flux(p)%z),'scalar register flux rollback')
        deallocate(scalar_flux(p)%x,scalar_flux(p)%y,scalar_flux(p)%z)
      enddo
      deallocate(scalar%fields)
    enddo
  end subroutine check_register_close

  subroutine compare_interface_repair(mesh,ownership,metrics,state,temperature,faces,accepted, &
      boundary_faces,minimum_level,rereflux_delta,rereflux_first_level)
    type(amr_patch_tree_3d), intent(in) :: mesh
    integer, intent(in) :: ownership(:)
    type(eb_geometry_3d), intent(in) :: metrics(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_integrated_flux_3d), intent(inout) :: faces(:)
    logical, intent(out) :: accepted
    character(len=*), intent(in), optional :: boundary_faces(6)
    integer, intent(in), optional :: minimum_level,rereflux_first_level
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    type(amr_tree_state_3d), allocatable :: serial_delta
    type(eb_reactive_tree_3d) :: serial
    type(eb_integrated_flux_3d), allocatable :: serial_faces(:)
    integer :: p,d(3),status
    logical :: serial_ok
    integer :: level
    level=0
    if(present(minimum_level)) level=minimum_level
    if(present(rereflux_delta)) then
      allocate(serial_delta)
      allocate(serial_delta%fields(size(ownership)))
    endif
    serial%topology=mesh
    allocate(serial%nodes(size(ownership)),serial_faces(size(ownership)))
    do p=1,size(ownership)
      d=mesh%nodes(p)%shape;serial%nodes(p)%geometry=metrics(p)
      allocate(serial%nodes(p)%state(nv,d(1),d(2),d(3)),serial%nodes(p)%temperature(d(1),d(2),d(3)))
      allocate(serial_faces(p)%x(nv,0:d(1),d(2),d(3)),serial_faces(p)%y(nv,d(1),0:d(2),d(3)), &
        serial_faces(p)%z(nv,d(1),d(2),0:d(3)))
      if(rank==ownership(p)) then
        serial%nodes(p)%state=state%fields(p)%values
        serial%nodes(p)%temperature=temperature%fields(p)%values(1,:,:,:)
        serial_faces(p)=faces(p)
      endif
      call MPI_Bcast(serial%nodes(p)%state,size(serial%nodes(p)%state),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'serial control gas broadcast')
      call MPI_Bcast(serial%nodes(p)%temperature,size(serial%nodes(p)%temperature), &
        MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'serial control temperature broadcast')
      call MPI_Bcast(serial_faces(p)%x,size(serial_faces(p)%x),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'serial control x transfer broadcast')
      call MPI_Bcast(serial_faces(p)%y,size(serial_faces(p)%y),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'serial control y transfer broadcast')
      call MPI_Bcast(serial_faces(p)%z,size(serial_faces(p)%z),MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'serial control z transfer broadcast')
      if(present(rereflux_delta)) then
        allocate(serial_delta%fields(p)%values(nv,d(1),d(2),d(3)))
        if(rank==ownership(p)) serial_delta%fields(p)%values=rereflux_delta%fields(p)%values
        call MPI_Bcast(serial_delta%fields(p)%values,size(serial_delta%fields(p)%values), &
          MPI_DOUBLE_PRECISION,ownership(p),MPI_COMM_WORLD,status)
        call require(status==MPI_SUCCESS,'serial control register broadcast')
      endif
    enddo
    call repair_reactive_eb_tree_interface_budget_3d(species,serial,serial_faces,serial_ok,level,boundary_faces, &
      rereflux_delta=serial_delta,rereflux_first_level=rereflux_first_level)
    call repair_mpi_eb_tree_interface_budget_3d(species,mesh,ownership,metrics,state,temperature,faces, &
      MPI_COMM_WORLD,accepted,level,boundary_faces,rereflux_delta=rereflux_delta, &
      rereflux_first_level=rereflux_first_level)
    call require(accepted.eqv.serial_ok,'serial/MPI donor repair acceptance agrees')
    do p=1,size(ownership)
      if(rank/=ownership(p)) cycle
      call require(all(state%fields(p)%values==serial%nodes(p)%state), &
        'serial/MPI repaired and rolled-back conserved states agree exactly')
      call require(all(temperature%fields(p)%values(1,:,:,:)==serial%nodes(p)%temperature), &
        'serial/MPI repaired and rolled-back temperatures agree exactly')
      call require(all(faces(p)%x==serial_faces(p)%x).and.all(faces(p)%y==serial_faces(p)%y).and. &
        all(faces(p)%z==serial_faces(p)%z),'serial/MPI repaired and rolled-back path faces agree exactly')
    enddo
  end subroutine compare_interface_repair

  subroutine check_peer_coarse_donor()
    type(amr_patch_3d) :: patch(4)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: metrics(5)
    type(amr_tree_state_3d) :: state,temperature,pending
    type(eb_integrated_flux_3d), allocatable :: faces(:),canonical(:)
    real(dp) :: primitive(reactive_nprim(size(species))),vector(nv),tcell,sound,donor_before,target_before
    real(dp) :: volume,packet,target_start(nv),donor_start(nv),balance(nv),donor_face_start
    real(dp), allocatable :: covered_start(:,:,:,:)
    integer, allocatable :: ownership(:)
    integer :: p,d(3),j,k,seam,reach,cover,node_count,target_y,donor_y,fine_y,target_face,donor_face,mpi_status
    character(len=16) :: boundary_faces(6)
    logical :: good
    do seam=0,1
    do reach=1,2
    do cover=0,2
    if(cover>0.and.seam/=0) cycle
    if(cover==1.and.reach/=2) cycle
    if(cover==2.and.reach/=1) cycle
    node_count=4+merge(1,0,cover>0)
    boundary_faces='outflow'
    if(seam==1) boundary_faces(3:4)='periodic'
    call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,4,2,patch(1),good)
    call require(good,'peer coarse first patch')
    call initialize_amr_patch_3d(4,4,4,1,2,3,4,1,4,2,patch(2),good)
    call require(good,'peer coarse second patch')
    call initialize_amr_patch_3d(4,4,8,3,3,merge(1,3,seam==1),merge(2,4,seam==1),3,4,2,patch(3),good)
    call require(good,'peer coarse refined child')
    if(cover>0) then
      call initialize_amr_patch_3d(4,4,8,4,4,1,1,3,3,2,patch(4),good)
      call require(good,'covered coarse transit child')
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],patch,mesh,good)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patch(:3),mesh,good)
    endif
    call require(good,'peer coarse topology')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(1),good)
    call require(good,'peer coarse root geometry')
    call build_axis_plane_eb_geometry_3d(4,4,8,0._dp,.5_dp,0._dp,.5_dp,0._dp,1._dp,'x',-.1_dp,metrics(2),good)
    call require(good,'peer coarse first geometry')
    call build_axis_plane_eb_geometry_3d(4,4,8,0._dp,.5_dp,.5_dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(3),good)
    call require(good,'peer coarse second geometry')
    call build_axis_plane_eb_geometry_3d(2,4,4,.25_dp,.375_dp,merge(0._dp,.25_dp,seam==1), &
      merge(.25_dp,.5_dp,seam==1),.25_dp,.5_dp, &
      'x',-.1_dp,metrics(4),good)
    call require(good,'peer coarse refined geometry')
    if(cover>0) then
      call build_axis_plane_eb_geometry_3d(2,2,2,.375_dp,.5_dp,.5_dp,.625_dp,.25_dp,.375_dp, &
        'x',-.1_dp,metrics(5),good)
      call require(good,'covered coarse transit geometry')
    endif
    target_y=merge(1,4,seam==1);donor_y=merge(5-reach,reach,seam==1)
    fine_y=merge(1,3,seam==1);target_face=merge(0,4,seam==1);donor_face=merge(5-reach,reach-1,seam==1)
    allocate(ownership(node_count))
    ownership(:4)=[0,0,nranks-1,merge(1,0,nranks>2)]
    if(cover>0) ownership(5)=0
    allocate(state%fields(node_count),temperature%fields(node_count),faces(node_count))
    do p=1,node_count
      if(rank/=ownership(p)) cycle
      d=mesh%nodes(p)%shape
      allocate(state%fields(p)%values(nv,d(1),d(2),d(3)),temperature%fields(p)%values(1,d(1),d(2),d(3)))
      allocate(faces(p)%x(nv,0:d(1),d(2),d(3)),faces(p)%y(nv,d(1),0:d(2),d(3)), &
        faces(p)%z(nv,d(1),d(2),0:d(3)))
      faces(p)%x=0;faces(p)%y=0;faces(p)%z=0
      primitive=0;primitive(1)=merge(3._dp,1._dp,p==4);primitive(5)=1.e6_dp
      primitive(6)=merge(.1_dp,0._dp,p==1);primitive(size(primitive))=1-primitive(6)
      call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
      call require(good,'peer coarse initial gas')
      do k=1,d(3)
      do j=1,d(2)
        state%fields(p)%values(:,:,j,k)=spread(vector,2,d(1))
        temperature%fields(p)%values(1,:,j,k)=tcell
      enddo
      enddo
    enddo
    if(rank==ownership(2)) then
      primitive=0;primitive(1)=1;primitive(5)=1.e6_dp;primitive(6)=.1_dp
      primitive(size(primitive))=.9_dp
      call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
      call require(good,'peer coarse target gas')
      state%fields(2)%values(:,4,target_y,3)=vector
      temperature%fields(2)%values(1,4,target_y,3)=tcell
    endif
    if(rank==ownership(4)) then
      faces(4)%x(6,2,fine_y:fine_y+1,1:2)=-.05_dp
      faces(4)%x(nv,2,fine_y:fine_y+1,1:2)=.05_dp
    endif
    if(seam==0.and.reach==1.and.rank==ownership(3)) then
      ! Coarse peer canonicalization leaves a pending donor-cell divergence.
      faces(3)%y(6,4,0,3)=.001_dp
      faces(3)%y(nv,4,0,3)=-.001_dp
    endif
    call compare_interface_repair(mesh,ownership,metrics(:node_count),state,temperature,faces,good,boundary_faces,1)
    call require(.not.good,'empty peer coarse budget rejects without publishing changes')
    if(rank==ownership(3)) then
      primitive=0;primitive(1)=1;primitive(5)=1.2e6_dp;primitive(6)=.8_dp
      primitive(size(primitive))=.2_dp
      call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
      call require(good,'peer coarse donor gas')
      state%fields(3)%values(:,4,donor_y,3)=vector
      temperature%fields(3)%values(1,4,donor_y,3)=tcell
      donor_before=vector(6)
      donor_start=vector
      donor_face_start=faces(3)%y(6,4,donor_face,3)
    endif
    if(cover==2) then
    if(rank==ownership(5)) then
      primitive=0;primitive(1)=1;primitive(5)=1.2e6_dp;primitive(6)=.8_dp
      primitive(size(primitive))=.2_dp
      call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
      call require(good,'covered coarse physical donor gas')
      do k=1,2
      do j=1,2
        state%fields(5)%values(:,:,j,k)=spread(vector,2,2)
        temperature%fields(5)%values(1,:,j,k)=tcell
      enddo
      enddo
    endif
    endif
    if(rank==ownership(2)) then
      target_before=state%fields(2)%values(6,4,target_y,3)
      target_start=state%fields(2)%values(:,4,target_y,3)
    endif
    if(cover>0) then
      if(rank==ownership(5)) covered_start=state%fields(5)%values
    endif
    if(seam==0.and.reach==1) then
      call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics(:node_count),state,temperature,faces,pending, &
        canonical,MPI_COMM_WORLD,good,1,2,boundary_faces)
      call require(good,'pending coarse donor reflux assembly')
      if(cover==0.and.rank==ownership(3)) call require(pending%fields(3)%values(6,4,donor_y,3)/=0, &
        'peer donor has a pending canonical face correction')
    endif
    call compare_interface_repair(mesh,ownership,metrics(:node_count),state,temperature,faces,good,boundary_faces,1)
    call require(good,'peer coarse donor repairs empty refined interface')
    balance=0
    if(rank==ownership(2)) then
      volume=metrics(2)%dx*metrics(2)%dy*metrics(2)%dz
      balance=balance+(state%fields(2)%values(:,4,target_y,3)-target_start)*volume
      call require(state%fields(2)%values(6,4,target_y,3)>target_before,'peer coarse target receives H2')
      if(seam==0) call require(faces(2)%y(6,4,target_face,3)<0,'peer coarse target records packet')
      if(seam==1) call require(faces(2)%y(6,4,target_face,3)>0,'peer coarse periodic target records packet')
    endif
    if(rank==ownership(3)) then
      call require(state%fields(3)%values(6,4,donor_y,3)<donor_before,'peer coarse donor debited')
      if(seam==0) call require(faces(3)%y(6,4,donor_face,3)<0,'peer coarse mirrored face records packet')
      if(seam==1) call require(faces(3)%y(6,4,donor_face,3)>0,'peer coarse periodic mirror records packet')
      volume=metrics(3)%dx*metrics(3)%dy*metrics(3)%dz
      if(cover/=2) balance=balance+(state%fields(3)%values(:,4,donor_y,3)-donor_start)*volume
      packet=(donor_before-state%fields(3)%values(6,4,donor_y,3))*volume
      call require(abs(packet+real(1-2*seam,dp)*(faces(3)%y(6,4,donor_face,3)-donor_face_start)* &
        metrics(3)%dx*metrics(3)%dz)<64*epsilon(1._dp), &
        'peer coarse debit equals face packet')
    endif
    if(cover>0) then
    if(rank==ownership(5)) then
      volume=metrics(5)%dx*metrics(5)%dy*metrics(5)%dz
      balance=balance+sum(sum(sum(state%fields(5)%values-covered_start,dim=4),dim=3),dim=2)*volume
      if(cover==1) call require(any(abs(faces(5)%y(6,:,0,:))>0).and.any(abs(faces(5)%y(6,:,2,:))>0), &
        'covered coarse transit records descendant entry and exit faces')
      if(cover==2) call require(any(abs(faces(5)%y(6,:,0,:))>0).and. &
        all(state%fields(5)%values(6,:,:,:)<covered_start(6,:,:,:)), &
        'covered coarse donor debits physical descendants and records entry face')
    endif
    endif
    call MPI_Allreduce(MPI_IN_PLACE,balance,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_status)
    call require(mpi_status==MPI_SUCCESS.and.maxval(abs(balance))<64*epsilon(1._dp), &
      'peer coarse path conserves every endpoint component')
    call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics(:node_count),state,temperature,faces,MPI_COMM_WORLD, &
      good,1,2,boundary_faces=boundary_faces)
    call require(good,'peer coarse funded interface passes strict synchronization')
    deallocate(state%fields,temperature%fields,faces,ownership)
    if(allocated(covered_start)) deallocate(covered_start)
    enddo
    enddo
    enddo
  end subroutine check_peer_coarse_donor

  subroutine check_refined_interface_repair()
    type(amr_patch_3d) :: patch(2)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: metrics(3),finest,coarse_fine
    type(amr_tree_state_3d) :: state,temperature,kept,kept_t
    type(eb_integrated_flux_3d), allocatable :: transfers(:),kept_transfers(:)
    real(dp) :: primitive(reactive_nprim(size(species))),vector(nv),tcell,sound,initial(nv),final(nv),expected(nv)
    real(dp) :: volume,fuel,raw_vector(nv),background(nv),phi(0:16,0:16,0:16)
    integer :: ownership(3),p,d(3),mode,limited,i,j,k,status,geometry_case
    logical :: good
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),good)
    call require(good,'refined repair fine patch')
    call initialize_amr_patch_3d(4,8,8,4,4,5,6,3,4,2,patch(2),good)
    call require(good,'refined repair finest patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patch,mesh,good)
    call require(good,'refined repair three-level topology')
    do geometry_case=0,1
    if(geometry_case==0) then
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(1),good)
    call require(good,'refined repair root metrics')
    call build_axis_plane_eb_geometry_3d(4,8,8,0._dp,.5_dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(2),good)
    call require(good,'refined repair fine metrics')
    call build_axis_plane_eb_geometry_3d(2,4,4,.375_dp,.5_dp,.5_dp,.75_dp,.25_dp,.5_dp,'x',-.1_dp,metrics(3),good)
    call require(good,'refined repair finest metrics')
    else
      do k=0,16
      do j=0,16
      do i=0,16
        phi(i,j,k)=real(i,dp)/16+.2_dp*real(j,dp)/16-.53_dp
      enddo
      enddo
      enddo
      call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],finest,good)
      call require(good,'sloped finest geometry')
      call coarsen_eb_geometry_3d(finest,2,coarse_fine,good)
      call require(good,'sloped fine geometry')
      call coarsen_eb_geometry_3d(coarse_fine,2,metrics(1),good)
      call require(good,'sloped root geometry')
      call extract_eb_geometry_3d(coarse_fine,[1,1,1],[4,8,8],metrics(2),good)
      call require(good,'sloped fine patch metrics')
      call extract_eb_geometry_3d(finest,[7,9,5],[8,12,8],metrics(3),good)
      call require(good,'sloped finest patch metrics')
      call require(any(metrics(3)%volume_fraction>0.and.metrics(3)%volume_fraction<1), &
        'sloped propagation fixture contains actual cut cells')
    endif
    raw_vector=0;raw_vector(1)=.001_dp;raw_vector(2)=.002_dp;raw_vector(5)=.3_dp
    raw_vector(6)=-.05_dp;raw_vector(nv)=.051_dp
    do mode=1,2
    do limited=0,5
      if(limited==4.and.geometry_case/=0) cycle
      background=0
      if(limited==5) then
        background(1)=.03125_dp;background(2)=.00390625_dp;background(5)=1024._dp
        background(6)=.0078125_dp;background(nv)=.0234375_dp
      endif
      raw_vector=0;raw_vector(1)=.001_dp;raw_vector(2)=.002_dp;raw_vector(5)=.3_dp
      raw_vector(6)=-.05_dp;raw_vector(nv)=.051_dp
      if(limited==4) then
        raw_vector(6)=-.03125_dp;raw_vector(9)=-.03125_dp;raw_vector(nv)=.0635_dp
      endif
      ownership=[0,nranks-1,merge(1,0,nranks>2)]
      if(mode==2) ownership=0
      allocate(state%fields(3),temperature%fields(3),transfers(3))
      do p=1,3
        if(rank/=ownership(p)) cycle
        d=mesh%nodes(p)%shape
        allocate(state%fields(p)%values(nv,d(1),d(2),d(3)),temperature%fields(p)%values(1,d(1),d(2),d(3)))
        allocate(transfers(p)%x(nv,0:d(1),d(2),d(3)),transfers(p)%y(nv,d(1),0:d(2),d(3)), &
          transfers(p)%z(nv,d(1),d(2),0:d(3)))
        transfers(p)%x=0;transfers(p)%y=0;transfers(p)%z=0
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          fuel=merge(.1_dp,.4_dp,p==1)
          if(limited>=2.and.p==2) fuel=.8_dp
          if(p==3) then
            fuel=merge(.25_dp,.55_dp,i==1)
            if(limited>=1) fuel=merge(.05_dp,.75_dp,i==1)
            if(geometry_case==1) then
              fuel=merge(.5_dp,.7_dp,i==1)
              if(i==1.and.j==4) fuel=.23_dp
              if(limited>=1) fuel=merge(.05_dp,.9_dp,i==1)
            endif
          endif
          primitive=0;primitive(1)=1;primitive(5)=1.e6_dp
          if(geometry_case==1.and.p==3) primitive(5)=1.2e6_dp
          if(limited>=2.and.p==3) primitive(5)=1.2e6_dp
          if(limited>=2.and.p==2) then
            primitive(1)=2;primitive(5)=4.e6_dp
            if(limited==3.and.i==3) primitive(5)=2.4e6_dp
          endif
          primitive(6)=fuel;primitive(size(primitive))=1-fuel
          if(limited==4) then
            primitive(1)=1;primitive(5)=2.e6_dp
            if(p==1) primitive(5)=1.e6_dp
            primitive(6)=.0625_dp;primitive(9)=.0625_dp
            if(p==2) then
              primitive(5)=4.e6_dp
              primitive(6)=.5625_dp;primitive(9)=.1875_dp
              if(i<=2) then
                primitive(5)=2.e6_dp
                primitive(6)=.1875_dp;primitive(9)=.5625_dp
              endif
            endif
            if(p==3) then
              primitive(5)=1.2e6_dp
              primitive(6)=merge(.0625_dp,.3125_dp,i==1)
              primitive(9)=primitive(6)
            endif
            primitive(size(primitive))=1-(primitive(6)+primitive(9))
          endif
          call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
          call require(good,'refined repair actual gas fixture')
          state%fields(p)%values(:,i,j,k)=vector;temperature%fields(p)%values(1,i,j,k)=tcell
        enddo
        enddo
        enddo
        if(p==1.and.limited==5) transfers(p)%x(:,2,3,2)=background
        if(p==2) then
          do k=3,4
          do j=5,6
            transfers(p)%x(:,4,j,k)=raw_vector+background
          enddo
          enddo
        endif
        if(p==3) then
          do k=1,d(3)
          do j=1,d(2)
            transfers(p)%x(:,2,j,k)=raw_vector+background
          enddo
          enddo
        endif
      enddo
      call average_down_mpi_eb_tree_3d(species,mesh,ownership,metrics,state,temperature,MPI_COMM_WORLD,good)
      call require(good,'initial finest states average down before interface correction')
      if(geometry_case==1.and.limited==0) then
        if(rank==ownership(2)) call require(.23_dp<.2_dp/metrics(2)%volume_fraction(4,6,3), &
          'minimum concentration with a uniform parent increment would reject the sloped fixture')
        if(rank==ownership(3)) call require(.23_dp>.4_dp/sum(metrics(3)%volume_fraction(:,4,1)), &
          'actual volume-weighted column increment admits the sloped low-fuel cell')
      endif
      kept=state;kept_t=temperature;kept_transfers=transfers
      initial=0
      do p=1,3
        if(rank/=ownership(p)) cycle
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        d=mesh%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(p==1.and.i<=2) cycle
          if(p==2.and.i==4.and.j>=5.and.j<=6.and.k>=3.and.k<=4) cycle
          initial=initial+state%fields(p)%values(:,i,j,k)*metrics(p)%volume_fraction(i,j,k)*volume
        enddo
        enddo
        enddo
      enddo
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,transfers,good)
      if(limited==1) then
        call require(.not.good,'finest inventory bounds reject unsafe parent-average correction')
        do p=1,3
          if(rank/=ownership(p)) cycle
          call require(all(state%fields(p)%values==kept%fields(p)%values).and. &
            all(temperature%fields(p)%values==kept_t%fields(p)%values),'refined failed repair preserves every level gas')
          call require(all(transfers(p)%x==kept_transfers(p)%x).and.all(transfers(p)%y==kept_transfers(p)%y).and. &
            all(transfers(p)%z==kept_transfers(p)%z),'refined failed repair preserves every level face')
        enddo
      else
        call require(good,'refined repair propagates feasible correction to finest level')
        if(limited==2.and.geometry_case==0.and.rank==ownership(2)) call require( &
          all(state%fields(2)%values(1,3,5:6,3:4)>1.5_dp), &
          'covered recipient funding retains donor inventory beyond full equalization')
        if(limited==3.and.rank==ownership(2)) call require( &
          all(state%fields(2)%values(:,3,:,:)==kept%fields(2)%values(:,3,:,:)), &
          'future recipient EOS rejects cold donor and preserves its intermediate-cell inventory')
        if(limited==4.and.rank==ownership(2)) then
          call require(any(state%fields(2)%values(6,3,:,:)<kept%fields(2)%values(6,3,:,:)), &
            'H2 donor contributes while O2 remains the limiting recipient constraint')
          call require(any(state%fields(2)%values(9,2,:,:)<kept%fields(2)%values(9,2,:,:)), &
            'complementary O2 donor completes the covered recipient budget')
        endif
        if(rank==ownership(3)) then
          call require(any(state%fields(3)%values/=kept%fields(3)%values),'finest actual inventory changes')
          call require(all(transfers(3)%x(6,1,:,:)>0),'finest internal faces carry distributed correction')
          if(limited>=2) call require(any(transfers(3)%x(6,0,:,:)/=0).or. &
            any(transfers(3)%y(6,:,0,:)/=0).or.any(transfers(3)%y(6,:,4,:)/=0).or. &
            any(transfers(3)%z(6,:,:,0)/=0).or.any(transfers(3)%z(6,:,:,4)/=0), &
            'covered recipient receives donor supply through physical boundary faces')
          call require(all(transfers(3)%x(6,2,:,:)>kept_transfers(3)%x(6,2,:,:)), &
            'finest boundary faces carry distributed correction')
        endif
        call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,transfers,MPI_COMM_WORLD, &
          good,0,2)
        call require(good,'refined correction passes all-level strict synchronization')
        final=0
        do p=1,3
          if(rank/=ownership(p)) cycle
          call require(all(state%fields(p)%values(6:,:,:,:)>=0),'refined accepted species remain nonnegative')
          volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
          d=mesh%nodes(p)%shape
          do k=1,d(3)
          do j=1,d(2)
          do i=1,d(1)
            if(p==1.and.i<=2) cycle
            if(p==2.and.i==4.and.j>=5.and.j<=6.and.k>=3.and.k<=4) cycle
            final=final+state%fields(p)%values(:,i,j,k)*metrics(p)%volume_fraction(i,j,k)*volume
          enddo
          enddo
          enddo
        enddo
        final=final-initial
        call MPI_Allreduce(MPI_IN_PLACE,final,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        call require(status==MPI_SUCCESS,'refined inventory collective')
        call MPI_Allreduce(MPI_IN_PLACE,initial,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        call require(status==MPI_SUCCESS,'refined initial inventory collective')
        expected=raw_vector*.25_dp**2
        call require(maxval(abs(final-expected)/max(1._dp,abs(initial)))<64*epsilon(1._dp), &
          'refined correction conserves density, species, momentum and energy after average-down')
      endif
      deallocate(state%fields,temperature%fields,transfers)
    enddo
    enddo
    enddo
  end subroutine check_refined_interface_repair

  subroutine check_refined_corridor_repair()
    type(amr_patch_3d) :: patch(2)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: metrics(3),regular_metrics(3),finest,coarse_fine
    type(amr_tree_state_3d) :: state,temperature,kept,kept_t
    type(eb_integrated_flux_3d), allocatable :: transfers(:),kept_transfers(:)
    real(dp) :: primitive(reactive_nprim(size(species))),vector(nv),tcell,sound,initial(nv),final(nv),expected(nv)
    real(dp) :: volume,fuel,raw_vector(nv),phi(0:16,0:16,0:16)
    integer :: ownership(3),p,d(3),mode,limited,i,j,k,status,donor_lower,funded_parents,covered_y_high
    logical :: good
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),good)
    call require(good,'refined repair fine patch')
    call initialize_amr_patch_3d(4,8,8,3,3,5,6,3,4,2,patch(2),good)
    call require(good,'refined repair finest patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1,2],patch,mesh,good)
    call require(good,'refined repair three-level topology')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(1),good)
    call require(good,'refined repair root metrics')
    call build_axis_plane_eb_geometry_3d(4,8,8,0._dp,.5_dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(2),good)
    call require(good,'refined repair fine metrics')
    call build_axis_plane_eb_geometry_3d(2,4,4,.25_dp,.375_dp,.5_dp,.75_dp,.25_dp,.5_dp,'x',-.1_dp,metrics(3),good)
    call require(good,'refined repair finest metrics')
    regular_metrics=metrics
    raw_vector=0;raw_vector(1)=.001_dp;raw_vector(2)=.002_dp;raw_vector(5)=.3_dp
    raw_vector(6)=-.05_dp;raw_vector(nv)=.051_dp
    do mode=1,2
    do limited=0,12
      metrics=regular_metrics;donor_lower=3;covered_y_high=6
      if(limited==12) then
        donor_lower=2;covered_y_high=7
      endif
      if(limited==9) donor_lower=2
      if(limited==10) donor_lower=1
      call initialize_amr_patch_3d(4,8,8,donor_lower,3,5,covered_y_high,3,4,2,patch(2),good)
      call require(good,'covered donor depth')
      call initialize_amr_patch_tree_3d([4,4,4],[1,2],patch,mesh,good)
      call require(good,'covered donor topology')
      if(limited==9.or.limited==10) then
        call build_axis_plane_eb_geometry_3d(2*(4-donor_lower),4,4,real(donor_lower-1,dp)*.125_dp,.375_dp, &
          .5_dp,.75_dp,.25_dp,.5_dp,'x',-.1_dp,metrics(3),good)
        call require(good,'pooled covered donor metrics')
      endif
      if((limited>=2.and.limited<=4).or.limited==8.or.limited==11.or.limited==12) then
        do k=0,16
        do j=0,16
        do i=0,16
          phi(i,j,k)=real(i,dp)/16+.2_dp*real(j,dp)/16-.38_dp
          if(limited==4.or.limited==11.or.limited==12) phi(i,j,k)=real(j,dp)/16+.2_dp*real(i,dp)/16-.60_dp
          if(limited==12) phi(i,j,k)=real(j,dp)/16+.2_dp*real(i,dp)/16-.72_dp
        enddo
        enddo
        enddo
        call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],finest,good)
        call require(good,'turning corridor nodal geometry')
        call coarsen_eb_geometry_3d(finest,2,coarse_fine,good)
        call require(good,'turning corridor fine metrics')
        call coarsen_eb_geometry_3d(coarse_fine,2,metrics(1),good)
        call require(good,'turning corridor root metrics')
        call extract_eb_geometry_3d(coarse_fine,[1,1,1],[4,8,8],metrics(2),good)
        call require(good,'turning corridor fine patch metrics')
        call extract_eb_geometry_3d(finest,[2*donor_lower-1,9,5],[6,2*covered_y_high,8],metrics(3),good)
        call require(good,'turning corridor descendant metrics')
        call require(any(metrics(3)%volume_fraction>0.and.metrics(3)%volume_fraction<1), &
          'turning corridor includes physical cut descendants')
      endif
      ownership=[0,nranks-1,merge(1,0,nranks>2)]
      if(mode==2) ownership=0
      allocate(state%fields(3),temperature%fields(3),transfers(3))
      do p=1,3
        if(rank/=ownership(p)) cycle
        d=mesh%nodes(p)%shape
        allocate(state%fields(p)%values(nv,d(1),d(2),d(3)),temperature%fields(p)%values(1,d(1),d(2),d(3)))
        allocate(transfers(p)%x(nv,0:d(1),d(2),d(3)),transfers(p)%y(nv,d(1),0:d(2),d(3)), &
          transfers(p)%z(nv,d(1),d(2),0:d(3)))
        transfers(p)%x=0;transfers(p)%y=0;transfers(p)%z=0
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          fuel=.1_dp
          if(p==3.and.limited==2) fuel=.25_dp
          if(p==3.and.limited==3) fuel=.01_dp
          if(p==3.and.(limited==4.or.limited==11)) fuel=.001_dp
          if(p==3.and.(limited==5.or.limited==8)) fuel=.8_dp
          if(p==3.and.limited==6) fuel=merge(.05_dp,.75_dp,i==1)
          if(p==3.and.limited==7) fuel=merge(.25_dp,.55_dp,i==1)
          if(p==3.and.(limited==9.or.limited==10)) fuel=merge(.06_dp,.74_dp,mod(i,2)==1)
          if(p==2.and.i==2.and.j>=5.and.j<=6.and.k>=3.and.k<=4.and.(limited==4.or.limited==11)) fuel=.8_dp
          if(p==2.and.i==4.and.j==7.and.k>=3.and.k<=4.and.limited==4) fuel=.8_dp
          if(p==2.and.i==3.and.j==4.and.k>=3.and.k<=4.and.(limited==2.or.limited==3)) fuel=.8_dp
          if(p==2.and.i==2.and.j>=5.and.j<=6.and.k>=3.and.k<=4.and.limited==0) fuel=.8_dp
          if(p==3.and.limited==12) then
            fuel=.001_dp
            if(i<=2.and.j>=3.and.j<=4.and.k<=2) fuel=.8_dp
          endif
          primitive=0;primitive(1)=1;primitive(5)=1.e6_dp
          if(fuel>=.8_dp.or.(p==3.and.limited>=5.and.limited<=10)) primitive(5)=1.2e6_dp
          if(p==3.and.limited==9) primitive(5)=1.e6_dp
          primitive(6)=fuel;primitive(size(primitive))=1-fuel
          call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
          call require(good,'refined repair actual gas fixture')
          state%fields(p)%values(:,i,j,k)=vector;temperature%fields(p)%values(1,i,j,k)=tcell
        enddo
        enddo
        enddo
        if(p==2) then
          do k=3,4
          do j=5,6
            if(limited==12.and.j==5) cycle
            transfers(p)%x(:,4,j,k)=raw_vector
          enddo
          enddo
        endif
      enddo
      call average_down_mpi_eb_tree_3d(species,mesh,ownership,metrics,state,temperature,MPI_COMM_WORLD,good)
      call require(good,'initial finest states average down before interface correction')
      kept=state;kept_t=temperature;kept_transfers=transfers
      initial=0
      do p=1,3
        if(rank/=ownership(p)) cycle
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        d=mesh%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(p==1.and.i<=2) cycle
          if(p==2.and.i>=donor_lower.and.i<=3.and.j>=5.and.j<=covered_y_high.and.k>=3.and.k<=4) cycle
          initial=initial+state%fields(p)%values(:,i,j,k)*metrics(p)%volume_fraction(i,j,k)*volume
        enddo
        enddo
        enddo
      enddo
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,transfers,good)
      if(limited==1.or.limited==6) then
        call require(.not.good,'insufficient donor or descendant inventory rejects funding')
        do p=1,3
          if(rank/=ownership(p)) cycle
          call require(all(state%fields(p)%values==kept%fields(p)%values).and. &
            all(temperature%fields(p)%values==kept_t%fields(p)%values),'refined failed repair preserves every level gas')
          call require(all(transfers(p)%x==kept_transfers(p)%x).and.all(transfers(p)%y==kept_transfers(p)%y).and. &
            all(transfers(p)%z==kept_transfers(p)%z),'refined failed repair preserves every level face')
        enddo
      else
        call require(good,'donor funding crosses refined interior cells')
        if(rank==ownership(3).and.(limited==0.or.limited==2)) then
          call require(any(transfers(3)%x(6,1,:,:)/=0).and.any(transfers(3)%x(6,2,:,:)/=0), &
            'refined corridor carries packet on boundary and internal faces')
          if(limited==0) call require(any(transfers(3)%x(6,0,:,:)/=0), &
            'straight refined corridor carries packet on both boundary faces')
          if(limited==0) call require(all(state%fields(3)%values==kept%fields(3)%values), &
            'straight regular corridor leaves descendant inventories unchanged')
          if(limited==2) then
            call require(any(transfers(3)%y/=0),'turning refined corridor carries transverse packets')
            call require(any(state%fields(3)%values/=kept%fields(3)%values), &
              'turning refined corridor changes physical descendant inventories')
          endif
        endif
        if(limited>=5.and.limited<=10) then
          if(rank==ownership(3)) then
            call require(any(state%fields(3)%values(6,:,:,:)<kept%fields(3)%values(6,:,:,:)), &
              'covered donor debits physical descendant fuel inventory')
            call require(any(transfers(3)%x(6,1:mesh%nodes(3)%shape(1)-1,:,:)/=0).and. &
              any(transfers(3)%x(6,mesh%nodes(3)%shape(1),:,:)/=0), &
              'covered donor carries withdrawal on internal and incident faces')
            call require(all(transfers(3)%x(:,0,:,:)==0), &
              'covered donor leaves its untraversed outer boundary unchanged')
            if(limited==9.or.limited==10) then
              funded_parents=0
              do k=1,4,2
              do j=1,4,2
              do i=1,mesh%nodes(3)%shape(1),2
                if(sum(kept%fields(3)%values(6,i:i+1,j:j+1,k:k+1)- &
                  state%fields(3)%values(6,i:i+1,j:j+1,k:k+1))>64*epsilon(1._dp)) &
                  funded_parents=funded_parents+1
              enddo
              enddo
              enddo
              call require(funded_parents>=7,'pooling spends physical budgets from at least seven covered donors')
              if(limited==10) call require(minval(state%fields(3)%values(6,1::2,:,:))>1.e-10_dp, &
                'EOS-limited partial supply retains fuel before the positivity ceiling')
            endif
          endif
        endif
        if(limited==12.and.rank==ownership(3)) call require( &
          any(state%fields(3)%values(6,1:2,3:4,:)<kept%fields(3)%values(6,1:2,3:4,:)), &
          'covered donor supplies after independently inadmissible corridors are excluded')
        if(limited==3.and.rank==ownership(2)) call require( &
          any(state%fields(2)%values(6,3,4,3:4)<kept%fields(2)%values(6,3,4,3:4)), &
          'turning donor uses another path when its first covered turn is inadmissible')
        if(limited==3.or.limited==11) then
          if(limited==11.and.rank==ownership(2)) call require( &
            any(state%fields(2)%values(6,2,5:6,3:4)<kept%fields(2)%values(6,2,5:6,3:4)), &
            'same donor supplies through an alternative path around deficient descendants')
          if(rank==ownership(3)) then
            call require(all(state%fields(3)%values==kept%fields(3)%values), &
              'alternative route preserves physical descendant inventories')
            if(limited==11) call require(all(transfers(3)%x(:,:,1:2,:)==kept_transfers(3)%x(:,:,1:2,:)), &
              'alternative route preserves the inadmissible first-path descendant faces')
          endif
        endif
        if(limited==4) then
          if(rank==ownership(2)) call require( &
            any(state%fields(2)%values(6,2,5:6,3:4)<kept%fields(2)%values(6,2,5:6,3:4)).or. &
            any(state%fields(2)%values(6,4,7,3:4)<kept%fields(2)%values(6,4,7,3:4)), &
            'inadmissible refined route is bypassed before a feasible donor is spent')
          if(rank==ownership(3)) then
            call require(all(state%fields(3)%values(:,:,1:2,:)==kept%fields(3)%values(:,:,1:2,:)), &
              'rejected route trial preserves deficient descendant inventory')
            call require(all(transfers(3)%x(:,:,1:2,:)==kept_transfers(3)%x(:,:,1:2,:)), &
              'rejected route trial preserves deficient descendant faces')
          endif
        endif
        call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,transfers,MPI_COMM_WORLD, &
          good,0,2)
        call require(good,'refined corridor passes all-level strict synchronization')
        final=0
        do p=1,3
          if(rank/=ownership(p)) cycle
          call require(all(state%fields(p)%values(6:,:,:,:)>=0),'refined accepted species remain nonnegative')
          call require(all(state%fields(p)%values(1,:,:,:)>0),'refined donor accepted density stays strictly positive')
          volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
          d=mesh%nodes(p)%shape
          do k=1,d(3)
          do j=1,d(2)
          do i=1,d(1)
            if(p==1.and.i<=2) cycle
            if(p==2.and.i>=donor_lower.and.i<=3.and.j>=5.and.j<=covered_y_high.and.k>=3.and.k<=4) cycle
            final=final+state%fields(p)%values(:,i,j,k)*metrics(p)%volume_fraction(i,j,k)*volume
          enddo
          enddo
          enddo
        enddo
        final=final-initial
        call MPI_Allreduce(MPI_IN_PLACE,final,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        call require(status==MPI_SUCCESS,'refined inventory collective')
        call MPI_Allreduce(MPI_IN_PLACE,initial,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
        call require(status==MPI_SUCCESS,'refined initial inventory collective')
        expected=raw_vector*.25_dp**2
        if(limited==12) expected=raw_vector*metrics(2)%dy*metrics(2)%dz*sum(metrics(2)%x_face_fraction(4,6,3:4))
        call require(maxval(abs(final-expected)/max(1._dp,abs(initial)))<64*epsilon(1._dp), &
          'refined corridor conserves density, species, momentum and energy after average-down')
      endif
      deallocate(state%fields,temperature%fields,transfers)
    enddo
    enddo
  end subroutine check_refined_corridor_repair

  subroutine check_cross_patch_repair()
    type(amr_patch_3d) :: patch(2)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: metrics(3)
    type(amr_tree_state_3d) :: state,temperature,kept,kept_t,pending
    type(eb_integrated_flux_3d), allocatable :: transfers(:),kept_transfers(:),canonical(:)
    real(dp) :: primitive(reactive_nprim(size(species))),vector(nv),tcell,sound,initial(nv),final(nv),expected(nv)
    real(dp) :: volume,lower,upper
    integer :: ownership(3),p,d(3),mode,periodic_case,target_y,donor_y,donor_x,j,k,status
    character(len=16) :: boundary_faces(6)
    logical :: good
    call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,4,2,patch(1),good)
    call require(good,'cross-patch target geometry')
    call initialize_amr_patch_3d(4,4,4,1,2,3,4,1,4,2,patch(2),good)
    call require(good,'cross-patch donor geometry')
    call initialize_amr_patch_tree_3d([4,4,4],[1,1],patch,mesh,good)
    call require(good,'cross-patch topology')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(1),good)
    call require(good,'cross-patch root metrics')
    do p=2,3
      lower=real(p-2,dp)*.5_dp;upper=lower+.5_dp
      call build_axis_plane_eb_geometry_3d(4,4,8,0._dp,.5_dp,lower,upper,0._dp,1._dp,'x',-.1_dp,metrics(p),good)
      call require(good,'cross-patch fine metrics')
    enddo
    do mode=1,4
    do periodic_case=0,1
      ownership=[0,nranks-1,merge(1,0,nranks>2)]
      if(mode==2) ownership=0
      boundary_faces='outflow'
      target_y=3;donor_y=1;donor_x=4
      if(mode==4) donor_x=2
      if(periodic_case==1) then
        boundary_faces(3:4)='periodic';target_y=1;donor_y=3
      endif
      allocate(state%fields(3),temperature%fields(3),transfers(3))
      do p=1,3
        if(rank/=ownership(p)) cycle
        d=mesh%nodes(p)%shape
        allocate(state%fields(p)%values(nv,d(1),d(2),d(3)),temperature%fields(p)%values(1,d(1),d(2),d(3)))
        allocate(transfers(p)%x(nv,0:d(1),d(2),d(3)),transfers(p)%y(nv,d(1),0:d(2),d(3)), &
          transfers(p)%z(nv,d(1),d(2),0:d(3)))
        transfers(p)%x=0;transfers(p)%y=0;transfers(p)%z=0
        primitive=0;primitive(1)=1;primitive(5)=1.e6_dp
        primitive(6)=merge(.1_dp,0._dp,p==1);primitive(size(primitive))=1-primitive(6)
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'cross-patch initial gas')
        do k=1,d(3)
        do j=1,d(2)
          state%fields(p)%values(:,:,j,k)=spread(vector,2,d(1))
          temperature%fields(p)%values(1,:,j,k)=tcell
        enddo
        enddo
      enddo
      if(rank==ownership(2)) then
        transfers(2)%x(6,4,target_y:target_y+1,3:4)=-.05_dp
        transfers(2)%x(nv,4,target_y:target_y+1,3:4)=.05_dp
      endif
      if(mode==3) then
        ! Raw patch faces can disagree before synchronization. Their divergence
        ! must be included in the available donor/recipient inventory once.
        p=merge(2,3,periodic_case==1)
        if(rank==ownership(p)) then
          transfers(p)%y(6,4,0,3:4)=.001_dp
          transfers(p)%y(nv,4,0,3:4)=-.001_dp
          if(periodic_case==1) then
            primitive=0;primitive(1)=1;primitive(5)=1.e6_dp;primitive(6)=.008_dp
            primitive(size(primitive))=1-primitive(6)
            call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
            call require(good,'periodic raw recipient gas')
            do k=3,4
              state%fields(p)%values(:,4,1,k)=vector
              temperature%fields(p)%values(1,4,1,k)=tcell
            enddo
          endif
        endif
      endif
      kept=state;kept_t=temperature;kept_transfers=transfers
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,transfers,good,boundary_faces)
      call require(.not.good,'cross-patch exhausted inventory rejects transaction')
      do p=1,3
        if(rank/=ownership(p)) cycle
        call require(all(state%fields(p)%values==kept%fields(p)%values).and. &
          all(temperature%fields(p)%values==kept_t%fields(p)%values),'cross-patch failure preserves gas')
        call require(all(transfers(p)%x==kept_transfers(p)%x).and.all(transfers(p)%y==kept_transfers(p)%y).and. &
          all(transfers(p)%z==kept_transfers(p)%z),'cross-patch failure preserves all face transfers')
      enddo
      if(rank==ownership(3)) then
        primitive=0;primitive(1)=1;primitive(5)=1.2e6_dp;primitive(6)=.8_dp
        primitive(size(primitive))=.2_dp
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'cross-patch donor gas')
        do k=3,4
        do j=donor_y,donor_y+1
          state%fields(3)%values(:,donor_x,j,k)=vector
          temperature%fields(3)%values(1,donor_x,j,k)=tcell
        enddo
        enddo
      endif
      call assemble_mpi_eb_tree_reflux_3d(species,mesh,ownership,metrics,state,temperature,transfers,pending, &
        canonical,MPI_COMM_WORLD,good,0,1,boundary_faces)
      call require(good,'cross-patch canonical inventory baseline')
      initial=0
      do p=1,3
        if(rank/=ownership(p)) cycle
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        if(p==1) then
          initial=initial+sum(sum(sum(state%fields(p)%values(:,3:4,:,:),dim=4),dim=3),dim=2)*volume
        else
          initial=initial+sum(sum(sum(state%fields(p)%values+pending%fields(p)%values,dim=4),dim=3),dim=2)*volume
        endif
      enddo
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,transfers,good,boundary_faces)
      call require(good,'cross-patch donor repairs exhausted fine interface')
      if(rank==ownership(2)) then
        j=merge(0,4,periodic_case==1)
        if(mode==4) then
          call require(any(abs(transfers(2)%y(6,:,j,3:4))>0),'target patch records turning boundary route')
        else
          call require(all(abs(transfers(2)%y(6,4,j,3:4))>0),'target patch records boundary exchange')
        endif
      endif
      if(rank==ownership(3)) then
        j=merge(4,0,periodic_case==1)
        if(mode==4) then
          call require(any(abs(transfers(3)%y(6,:,j,3:4))>0),'remote donor patch records turning boundary route')
        else
          call require(all(abs(transfers(3)%y(6,4,j,3:4))>0),'remote donor patch records boundary exchange')
        endif
      endif
      call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,transfers,MPI_COMM_WORLD, &
        good,0,1,boundary_faces=boundary_faces)
      call require(good,'cross-patch funded gas passes strict synchronization')
      final=0
      do p=1,3
        if(rank/=ownership(p)) cycle
        call require(all(state%fields(p)%values(6:,:,:,:)>=0),'cross-patch accepted species nonnegative')
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        if(p==1) then
          final=final+sum(sum(sum(state%fields(p)%values(:,3:4,:,:),dim=4),dim=3),dim=2)*volume
        else
          final=final+sum(sum(sum(state%fields(p)%values,dim=4),dim=3),dim=2)*volume
        endif
      enddo
      final=final-initial
      call MPI_Allreduce(MPI_IN_PLACE,final,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'cross-patch inventory collective')
      call MPI_Allreduce(MPI_IN_PLACE,initial,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'cross-patch initial inventory collective')
      expected=0;expected(6)=-.2_dp*.25_dp**3;expected(nv)=-expected(6)
      call require(maxval(abs(final-expected)/max(1._dp,abs(initial)))<64*epsilon(1._dp), &
        'cross-patch ordinary and periodic paths conserve every component')
      deallocate(state%fields,temperature%fields,transfers)
    enddo
    enddo
  end subroutine check_cross_patch_repair

  subroutine check_interface_repair()
    type(amr_patch_3d) :: patch(1)
    type(amr_patch_tree_3d) :: mesh
    type(eb_geometry_3d) :: metrics(2),regular_metrics(2),global_fine
    type(amr_tree_state_3d) :: state,temperature,kept,kept_t,trial,trial_t
    type(eb_integrated_flux_3d), allocatable :: faces(:),kept_faces(:),trial_faces(:)
    real(dp) :: primitive(reactive_nprim(size(species))),vector(nv),tcell,sound,initial(nv),final(nv)
    real(dp) :: total(nv),expected(nv),volume,correction,barrier_phi(0:8,0:8,0:8)
    integer :: ownership(2),p,mode,d(3),i,j,k,status,donor_case,donor_x,donor_y,search_a,search_s,search_d,coordinate(3)
    logical :: good,common_support
    character(len=16) :: repair_faces(6)
    call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),good)
    call require(good,'joint repair child patch')
    call initialize_amr_patch_tree_3d([4,4,4],[1],patch,mesh,good)
    call require(good,'joint repair topology')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(1),good)
    call require(good,'joint repair coarse geometry')
    call build_axis_plane_eb_geometry_3d(4,8,8,0._dp,.5_dp,0._dp,1._dp,0._dp,1._dp,'x',-.1_dp,metrics(2),good)
    call require(good,'joint repair fine geometry')
    regular_metrics=metrics
    do mode=1,2
      ownership=[0,nranks-1]
      if(mode==2) ownership=0
      allocate(state%fields(2),temperature%fields(2),faces(2))
      initial=0
      do p=1,2
        if(ownership(p)/=rank) cycle
        d=mesh%nodes(p)%shape
        allocate(state%fields(p)%values(nv,d(1),d(2),d(3)),temperature%fields(p)%values(1,d(1),d(2),d(3)))
        primitive=0;primitive(1)=1;primitive(5)=1e6_dp
        primitive(6)=merge(.1_dp,.4_dp,p==1);primitive(size(primitive))=1-primitive(6)
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'joint repair physical gas state')
        do k=1,d(3)
        do j=1,d(2)
          state%fields(p)%values(:,:,j,k)=spread(vector,2,d(1))
        enddo
        enddo
        temperature%fields(p)%values=tcell
        allocate(faces(p)%x(nv,0:d(1),d(2),d(3)),faces(p)%y(nv,d(1),0:d(2),d(3)), &
          faces(p)%z(nv,d(1),d(2),0:d(3)))
        faces(p)%x=0;faces(p)%y=0;faces(p)%z=0
        if(p==2) then
          faces(p)%x(6,4,5:6,3:4)=-.05_dp
          faces(p)%x(nv,4,5:6,3:4)=.05_dp
        endif
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        if(p==1) then
          initial=initial+sum(sum(sum(state%fields(p)%values(:,3:4,:,:),dim=4),dim=3),dim=2)*volume
        else
          initial=initial+sum(sum(sum(state%fields(p)%values,dim=4),dim=3),dim=2)*volume
        endif
      enddo
      kept=state;kept_t=temperature;kept_faces=faces
      call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,faces,MPI_COMM_WORLD, &
        good,0,1)
      call require(.not.good,'unreconciled interface rejects negative coarse fuel')
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,faces,good)
      call require(good,'owner-local coarse/fine interface repair')
      do p=1,2
        call require(allocated(state%fields(p)%values).eqv.(ownership(p)==rank),'repair preserves idle-rank ownership')
        if(ownership(p)/=rank) cycle
        call require(all(state%fields(p)%values(6:,:,:,:)>=0),'repair publishes nonnegative fine species')
        if(p==1) call require(all(state%fields(p)%values==kept%fields(p)%values), &
          'repair defers coarse update to strict reflux')
        if(p==2) call require(all(faces(p)%x(6,4,5:6,3:4)>kept_faces(p)%x(6,4,5:6,3:4)), &
          'repair updates the matching fine boundary transfer')
      enddo
      call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,faces,MPI_COMM_WORLD, &
        good,0,1)
      call require(good,'repaired interface passes strict reflux acceptance')
      final=0
      do p=1,2
        if(ownership(p)/=rank) cycle
        call require(all(state%fields(p)%values(6:,:,:,:)>=0),'synchronized interface has nonnegative species')
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        if(p==1) then
          final=final+sum(sum(sum(state%fields(p)%values(:,3:4,:,:),dim=4),dim=3),dim=2)*volume
        else
          final=final+sum(sum(sum(state%fields(p)%values,dim=4),dim=3),dim=2)*volume
        endif
      enddo
      total=final-initial
      call MPI_Allreduce(MPI_IN_PLACE,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      expected=0;expected(6)=-.2_dp*.25_dp**3;expected(nv)=-expected(6)
      call require(status==MPI_SUCCESS.and.maxval(abs(total-expected))<64*epsilon(1._dp), &
        'joint face/state repair conserves the full composite reflux inventory')
      state=kept;temperature=kept_t;faces=kept_faces
      if(rank==ownership(2)) then
        primitive=0;primitive(1)=1;primitive(5)=1e6_dp;primitive(size(primitive))=1
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'exhausted fine fuel fixture')
        d=mesh%nodes(2)%shape
        do k=1,d(3)
        do j=1,d(2)
          state%fields(2)%values(:,:,j,k)=spread(vector,2,d(1))
        enddo
        enddo
        temperature%fields(2)%values=tcell
      endif
      kept=state;kept_t=temperature
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,faces,good)
      call require(.not.good,'incompatible coarse/fine interface requires wider donors')
      do p=1,2
        if(ownership(p)/=rank) cycle
        call require(all(state%fields(p)%values==kept%fields(p)%values).and. &
          all(temperature%fields(p)%values==kept_t%fields(p)%values),'failed joint repair preserves gas and temperature')
        call require(all(faces(p)%x==kept_faces(p)%x).and.all(faces(p)%y==kept_faces(p)%y).and. &
          all(faces(p)%z==kept_faces(p)%z),'failed joint repair preserves every face')
      enddo
      ! An exposed coarse neighbor can fund the same empty fine interface.
      state=kept;temperature=kept_t;faces=kept_faces
      if(rank==ownership(2)) then
        primitive=0;primitive(1)=3;primitive(5)=1.e6_dp
        primitive(size(primitive))=1
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'coarse donor fixture has a solvent budget on the fine side')
        d=mesh%nodes(2)%shape
        do k=1,d(3)
        do j=1,d(2)
          state%fields(2)%values(:,:,j,k)=spread(vector,2,d(1))
        enddo
        enddo
        temperature%fields(2)%values=tcell
      endif
      if(rank==ownership(1)) then
        primitive=0;primitive(1)=1;primitive(5)=1.2e6_dp
        primitive(6)=.8_dp;primitive(size(primitive))=.2_dp
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'exposed coarse donor gas fixture')
        state%fields(1)%values(:,3,2,2)=vector
        temperature%fields(1)%values(1,3,2,2)=tcell
      endif
      initial=0
      do p=1,2
        if(ownership(p)/=rank) cycle
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        d=mesh%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(p==1.and.i<=2) cycle
          initial=initial+state%fields(p)%values(:,i,j,k)*volume
        enddo
        enddo
        enddo
      enddo
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,faces,good)
      call require(good,'exposed coarse donor repairs empty fine interface')
      if(rank==ownership(1)) call require(state%fields(1)%values(6,3,2,2)<vector(6), &
        'coarse donor supplies H2 before strict reflux')
      if(rank==ownership(1)) call require(faces(1)%y(6,3,2,2)>0, &
        'coarse donor packet is recorded on the shared interior face')
      call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,faces,MPI_COMM_WORLD, &
        good,0,1)
      call require(good,'coarse-funded interface passes strict reflux')
      final=0
      do p=1,2
        if(ownership(p)/=rank) cycle
        call require(all(state%fields(p)%values(6:,:,:,:)>=0),'coarse-funded species remain nonnegative')
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        d=mesh%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(p==1.and.i<=2) cycle
          final=final+state%fields(p)%values(:,i,j,k)*volume
        enddo
        enddo
        enddo
      enddo
      total=final-initial
      call MPI_Allreduce(MPI_IN_PLACE,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS.and.maxval(abs(total-expected))<64*epsilon(1._dp), &
        'coarse-funded interface conserves every reflux inventory')
      ! Exercise both an inward donor and a two-face tangential path with
      ! shared intermediate faces. Each path must have zero internal inventory.
      do donor_case=1,7
      state=kept;temperature=kept_t;faces=kept_faces
      repair_faces="outflow"
      metrics=regular_metrics
      if(donor_case==6) then
        do k=0,8
        do j=0,8
        do i=0,8
          barrier_phi(i,j,k)=max(abs(real(i,dp)/8-.3125_dp),abs(real(j,dp)/8-.5625_dp), &
            abs(real(k,dp)/8-.3125_dp))-.066_dp
        enddo
        enddo
        enddo
        call build_level_set_eb_geometry_3d(barrier_phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],global_fine,good)
        call require(good,'closed obstacle fine geometry')
        call coarsen_eb_geometry_3d(global_fine,2,metrics(1),good)
        call require(good,'closed obstacle coarse geometry')
        call extract_eb_geometry_3d(global_fine,[1,1,1],[4,8,8],metrics(2),good)
        call require(good,'closed obstacle fine patch geometry')
        call require(metrics(2)%volume_fraction(3,5,3)==0.and.metrics(2)%x_face_fraction(3,5,3)==0, &
          'obstacle blocks the direct donor route with a covered cell and closed face')
      endif
      donor_y=merge(5,3,donor_case==1)
      if(donor_case==3) then
        donor_y=7;repair_faces(3:4)="periodic"
        if(rank==ownership(2)) then
          faces(2)%x(6,4,:,:)=0;faces(2)%x(nv,4,:,:)=0
          faces(2)%x(6,4,1:2,3:4)=-.05_dp
          faces(2)%x(nv,4,1:2,3:4)=.05_dp
        endif
      endif
      donor_x=merge(3,4,donor_case==1)
      if(rank==ownership(2)) then
        primitive=0;primitive(1)=1;primitive(5)=1.2e6_dp;primitive(6)=.8_dp
        primitive(size(primitive))=.2_dp
        call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
        call require(good,'interior fuel donor gas fixture')
        if(donor_case<=3) then
        do k=3,4
        do j=donor_y,donor_y+1
          state%fields(2)%values(:,donor_x,j,k)=vector
          temperature%fields(2)%values(1,donor_x,j,k)=tcell
        enddo
        enddo
        else if(donor_case==7) then
          primitive(6)=.3_dp;primitive(size(primitive))=.7_dp
          call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
          call require(good,'pooled weak donor gas fixture')
          do k=1,3
          do j=1,3
            state%fields(2)%values(:,2,j,k)=vector;temperature%fields(2)%values(1,2,j,k)=tcell
          enddo
          enddo
          call require(maxval(state%fields(2)%values(6,:,:,:))<.4_dp, &
            'every individual donor lacks the convex pair inventory for one exhausted target')
        else if(donor_case==6) then
          state%fields(2)%values(:,2,5,3)=vector;temperature%fields(2)%values(1,2,5,3)=tcell
          state%fields(2)%values(:,2,6,4)=vector;temperature%fields(2)%values(1,2,6,4)=tcell
        else
          if(donor_case==5) then
            primitive(6)=.41_dp;primitive(size(primitive))=.59_dp
            call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
            call require(good,'limited early donor gas fixture')
          endif
          state%fields(2)%values(:,3,5,3)=vector;temperature%fields(2)%values(1,3,5,3)=tcell
          state%fields(2)%values(:,2,6,3)=vector;temperature%fields(2)%values(1,2,6,3)=tcell
          if(donor_case==4) then
            state%fields(2)%values(:,4,6,6)=vector;temperature%fields(2)%values(1,4,6,6)=tcell
          else
            primitive(5)=1.e6_dp;primitive(6)=.01_dp;primitive(size(primitive))=.99_dp
            call reactive_primitive_to_conserved(species,primitive,vector,tcell,sound,good)
            call require(good,'insufficient last donor gas fixture')
          endif
          state%fields(2)%values(:,4,3,4)=vector;temperature%fields(2)%values(1,4,3,4)=tcell
        endif
      endif
      if(donor_case==4.and.rank==ownership(2)) then
        do search_a=1,3
        do search_s=-1,1,2
        do search_d=1,7
          common_support=.true.
          do k=3,4
          do j=5,6
            coordinate=[4,j,k];coordinate(search_a)=coordinate(search_a)+search_s*search_d
            if(any(coordinate<1).or.any(coordinate>mesh%nodes(2)%shape)) then
              common_support=.false.;cycle
            endif
            if(state%fields(2)%values(6,coordinate(1),coordinate(2),coordinate(3))<=0) common_support=.false.
          enddo
          enddo
          call require(.not.common_support,'mixed fixture has no common translated donor support for all subfaces')
        enddo
        enddo
        enddo
      endif
      initial=0
      do p=1,2
        if(ownership(p)/=rank) cycle
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        d=mesh%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(p==1.and.i<=2) cycle
          initial=initial+state%fields(p)%values(:,i,j,k)*metrics(p)%volume_fraction(i,j,k)*volume
        enddo
        enddo
        enddo
      enddo
      trial=state;trial_t=temperature;trial_faces=faces
      call compare_interface_repair(mesh,ownership,metrics,state,temperature,faces,good,repair_faces)
      if(donor_case==5) then
        call require(.not.good,'insufficient later donor rejects earlier private exchanges')
        do p=1,2
          if(rank/=ownership(p)) cycle
          call require(all(state%fields(p)%values==trial%fields(p)%values).and. &
            all(temperature%fields(p)%values==trial_t%fields(p)%values),'partial funding failure rolls back all gas')
          call require(all(faces(p)%x==trial_faces(p)%x).and.all(faces(p)%y==trial_faces(p)%y).and. &
            all(faces(p)%z==trial_faces(p)%z),'partial funding failure rolls back all path faces')
        enddo
        cycle
      endif
      call require(good,'connected interior donors repair exhausted fine boundary cells')
      if(rank==ownership(2)) then
        if(donor_case==1) call require(all(faces(2)%x(6,3,5:6,3:4)>0), &
          'inward donor exchange publishes its internal face flux')
        if(donor_case==2) call require(all(faces(2)%y(6,4,4:5,3:4)>0), &
          'tangential donor exchange publishes every shared internal path face')
        if(donor_case==4) then
          call require(all(faces(2)%x(6,2:3,6,3)>0),'independent inward route spans two faces')
          call require(faces(2)%y(6,4,4,4)>0,'independent tangential route funds another subface')
        endif
        if(donor_case==7) then
          call require(sum(merge(1,0,state%fields(2)%values(6,2,1:3,1:3)<.3_dp))>=2, &
            'pooled repair debits multiple individually insufficient donor cells')
        endif
        if(donor_case==6) then
          call require(state%fields(2)%values(6,2,5,3)<.8_dp,'connected route reaches donor behind the obstacle')
          call require(all(faces(2)%x(:,3,5,3)==0),'connected route leaves the closed obstacle face untouched')
          call require(any(faces(2)%y/=0).or.any(faces(2)%z/=0),'connected route records its turns on transverse faces')
        endif
        if(donor_case==3) then
          call require(all(faces(2)%y(6,4,0,3:4)>0),'periodic donor exchange crosses low seam')
          call require(all(faces(2)%y(:,4,0,3:4)==faces(2)%y(:,4,8,3:4)), &
            'periodic donor exchange publishes identical packets on both seam faces')
        endif
      endif
      call synchronize_mpi_eb_tree_flux_3d(species,mesh,ownership,metrics,state,temperature,faces,MPI_COMM_WORLD, &
        good,0,1,boundary_faces=repair_faces)
      call require(good,'interior-funded interface passes strict reflux acceptance')
      final=0
      do p=1,2
        if(ownership(p)/=rank) cycle
        call require(all(state%fields(p)%values(6:,:,:,:)>=0),'interior donor and recipient species stay nonnegative')
        volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz
        d=mesh%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(p==1.and.i<=2) cycle
          final=final+state%fields(p)%values(:,i,j,k)*metrics(p)%volume_fraction(i,j,k)*volume
        enddo
        enddo
        enddo
      enddo
      total=final-initial
      call MPI_Allreduce(MPI_IN_PLACE,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS,'interior exchange inventory collective')
      call MPI_Allreduce(MPI_IN_PLACE,initial,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(status==MPI_SUCCESS.and.maxval(abs(total-expected)/max(1._dp,abs(initial)))<64*epsilon(1._dp), &
        'interior-funded reflux conserves species, density, momentum and energy')
      enddo
      metrics=regular_metrics
      deallocate(state%fields,temperature%fields,faces)
    enddo
  end subroutine check_interface_repair

  subroutine check_mpi_rereflux()
    type(amr_tree_state_3d) :: local,lt,initial,it,correction,serial_correction,increment,serial_increment,kept
    type(eb_geometry_3d) :: geometry(count)
    type(eb_reactive_tree_3d) :: expected
    type(rotating_frame_3d) :: test_frame
    type(eb_integrated_flux_3d), allocatable :: lf(:),expected_flux(:),replacement(:),serial_replacement(:)
    integer :: owners(count),p,mode,first
    logical :: good
    test_frame=frame
    if(.not.rotating) test_frame%angular_speed=0
    allocate(serial_correction%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      serial_correction%fields(p)%values=1.e-7_dp*saved%nodes(p)%state
    enddo
    do mode=1,4
      first=0
      owners=[(mod(p-1,nranks),p=1,count)]
      if(mode==2) owners=0
      if(mode>=3) owners=nranks-1-owners
      if(mode==4) then
        first=1;serial_correction%fields(1)%values=0
      endif
      allocate(local%fields(count),lt%fields(count),correction%fields(count),lf(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=saved%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
        correction%fields(p)%values=serial_correction%fields(p)%values
        lf(p)=raw(p)
      enddo
      initial=local;it=lt;expected=saved;expected_flux=raw
      call assemble_reactive_eb_tree_reflux_3d(species,saved,raw,serial_increment,serial_replacement,good, &
        first,2,boundary_faces=faces,frame=test_frame,rereflux_delta=serial_correction)
      call require(good,'serial nested rereflux assembly')
      call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,local,lt,lf,increment,replacement, &
        MPI_COMM_WORLD,good,first,2,boundary_faces=faces,frame=test_frame,rereflux_delta=correction)
      call require(good,'MPI nested rereflux assembly')
      do p=1,count
        call require(allocated(increment%fields(p)%values).eqv.(owners(p)==rank),'rereflux retains owner-local allocation')
        if(owners(p)/=rank) cycle
        call require(all(increment%fields(p)%values==serial_increment%fields(p)%values), &
          'serial/MPI nested rereflux increments agree exactly')
      enddo
      call synchronize_reactive_eb_tree_flux_3d(species,expected,expected_flux,good,first,2, &
        boundary_faces=faces,frame=test_frame,rereflux_delta=serial_correction)
      call require(good,'serial nested rereflux EOS control')
      call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
        first,2,boundary_faces=faces,frame=test_frame,rereflux_delta=correction)
      call require(good,'MPI nested rereflux EOS acceptance')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==expected%nodes(p)%state).and. &
          all(lt%fields(p)%values(1,:,:,:)==expected%nodes(p)%temperature), &
          'serial/MPI nested rereflux states and temperatures agree exactly')
        call require(all(lf(p)%x==expected_flux(p)%x).and.all(lf(p)%y==expected_flux(p)%y).and. &
          all(lf(p)%z==expected_flux(p)%z),'serial/MPI nested rereflux effective fluxes agree exactly')
      enddo
      local=initial;lt=it
      do p=1,count
        if(owners(p)==rank) lf(p)=raw(p)
      enddo
      if(rank==owners(count)) correction%fields(count)%values=-2*saved%nodes(count)%state
      call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
        first,2,boundary_faces=faces,frame=test_frame,rereflux_delta=correction)
      call require(.not.good,'MPI descendant rereflux rejection is collective')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values).and. &
          all(lt%fields(p)%values==it%fields(p)%values),'MPI nested rereflux all gas rollback')
        call require(all(lf(p)%x==raw(p)%x).and.all(lf(p)%y==raw(p)%y).and. &
          all(lf(p)%z==raw(p)%z),'MPI nested rereflux all flux rollback')
      enddo
      if(rank==owners(count)) correction%fields(count)%values=ieee_value(0._dp,ieee_quiet_nan)
      kept=increment
      call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,local,lt,lf,increment,replacement, &
        MPI_COMM_WORLD,good,first,2,boundary_faces=faces,frame=test_frame,rereflux_delta=correction)
      call require(.not.good,'MPI invalid rereflux register is collective')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(increment%fields(p)%values==kept%fields(p)%values),'MPI invalid register assembly rollback')
      enddo
      if(rank==0) then
        call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
          first,2,boundary_faces=faces,frame=test_frame)
      else
        call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
          first,2,boundary_faces=faces,frame=test_frame,rereflux_delta=correction)
      endif
      call require(.not.good,'MPI rank-disagreeing rereflux presence rejected before communication')
      deallocate(local%fields,lt%fields,correction%fields,lf)
    enddo
  end subroutine check_mpi_rereflux

  subroutine check_mpi_sync()
    type(amr_tree_state_3d) :: local,lt,initial,it,increment,serial_increment,kept_increment
    type(amr_tree_state_3d) :: scalar,expected_scalar
    type(eb_geometry_3d) :: geometry(count)
    type(eb_reactive_tree_3d) :: expected
    type(eb_integrated_flux_3d), allocatable :: lf(:),expected_flux(:),saved_flux(:)
    type(eb_integrated_flux_3d), allocatable :: replacement(:),serial_replacement(:),kept_replacement(:)
    type(eb_integrated_flux_3d), allocatable :: scalar_flux(:),expected_scalar_flux(:)
    integer :: owners(count),p,mode,first,last,failure
    logical :: good
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
    end do
    do mode=1,3
      owners=[(mod(p-1,nranks),p=1,count)]
      if (mode==2) owners=0
      first=0;last=2
      if (mode==3) first=1
      expected=tree;expected_flux=flux
      if (mode==3) then
        expected=saved;expected_flux=raw
        if(rotating) then
          call synchronize_reactive_eb_tree_flux_3d(species,expected,expected_flux,good,first,last,boundary_faces=faces,frame=frame)
        else
          call synchronize_reactive_eb_tree_flux_3d(species,expected,expected_flux,good,first,last,boundary_faces=faces)
        end if
        call require(good,'serial level-scoped sync')
      end if
      allocate(local%fields(count),lt%fields(count),lf(count))
      do p=1,count
        if (owners(p)/=rank) cycle
        local%fields(p)%values=saved%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
        lf(p)=raw(p)
      end do
      initial=local;it=lt
      if(rotating) then
        call assemble_reactive_eb_tree_reflux_3d(species,saved,raw,serial_increment,serial_replacement,good, &
          first,last,boundary_faces=faces,frame=frame)
      else
        call assemble_reactive_eb_tree_reflux_3d(species,saved,raw,serial_increment,serial_replacement,good, &
          first,last,boundary_faces=faces)
      endif
      call require(good,'serial signed composite reflux assembly')
      if(rotating) then
        call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,local,lt,lf,increment,replacement, &
          MPI_COMM_WORLD,good,first,last,boundary_faces=faces,frame=frame)
      else
        call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,local,lt,lf,increment,replacement, &
          MPI_COMM_WORLD,good,first,last,boundary_faces=faces)
      endif
      call require(good,'MPI signed composite reflux assembly')
      do p=1,count
        call require(allocated(increment%fields(p)%values).eqv.(owners(p)==rank),'owner-local signed increments')
        call require(allocated(replacement(p)%x).eqv.(owners(p)==rank),'owner-local assembled faces')
        if(owners(p)/=rank) cycle
        call require(all(increment%fields(p)%values==serial_increment%fields(p)%values), &
          'exact serial/MPI signed reflux increments')
        call require(all(replacement(p)%x==serial_replacement(p)%x).and. &
          all(replacement(p)%y==serial_replacement(p)%y).and.all(replacement(p)%z==serial_replacement(p)%z), &
          'exact serial/MPI assembled replacements')
        call require(all(replacement(p)%x==expected_flux(p)%x).and.all(replacement(p)%y==expected_flux(p)%y).and. &
          all(replacement(p)%z==expected_flux(p)%z),'assembled faces equal accepted synchronization faces')
        call require(all(local%fields(p)%values==initial%fields(p)%values).and. &
          all(lt%fields(p)%values==it%fields(p)%values),'assembly leaves gas and temperature unchanged')
        call require(all(lf(p)%x==raw(p)%x).and.all(lf(p)%y==raw(p)%y).and.all(lf(p)%z==raw(p)%z), &
          'assembly leaves input integrated flux unchanged')
        if(topology%nodes(p)%level<first) call require(all(increment%fields(p)%values==0), &
          'unselected level has zero assembled reflux increment')
      enddo

      if(rotating) then
        call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
          first,last,boundary_faces=faces,frame=frame)
      else
        call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
          first,last,boundary_faces=faces)
      end if
      call require(good,'MPI flux sync')
      do p=1,count
        call require(allocated(lf(p)%x).eqv.(owners(p)==rank),'owner-local flux')
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==expected%nodes(p)%state),'exact serial reflux state')
        call require(all(lt%fields(p)%values(1,:,:,:)==expected%nodes(p)%temperature),'exact serial reflux T')
        call require(all(lf(p)%x==expected_flux(p)%x).and.all(lf(p)%y==expected_flux(p)%y).and. &
          all(lf(p)%z==expected_flux(p)%z),'exact serial effective flux')
      end do
      allocate(expected_scalar%fields(count),scalar%fields(count),scalar_flux(count),expected_scalar_flux(count))
      do p=1,count
        expected_scalar%fields(p)%values=saved%nodes(p)%state(5:5,:,:,:)
        allocate(expected_scalar_flux(p)%x(1,0:4,4,4),expected_scalar_flux(p)%y(1,4,0:4,4), &
          expected_scalar_flux(p)%z(1,4,4,0:4))
        expected_scalar_flux(p)%x=raw(p)%x(5:5,:,:,:)
        expected_scalar_flux(p)%y=raw(p)%y(5:5,:,:,:)
        expected_scalar_flux(p)%z=raw(p)%z(5:5,:,:,:)
        if(owners(p)/=rank) cycle
        scalar%fields(p)=expected_scalar%fields(p)
        scalar_flux(p)=expected_scalar_flux(p)
      enddo
      call synchronize_reactive_eb_tree_scalar_flux_3d(species,saved,expected_scalar,expected_scalar_flux, &
        good,first,last,boundary_faces=faces)
      call require(good,'serial scalar flux sync')
      call synchronize_mpi_eb_tree_scalar_flux_3d(species,topology,owners,geometry,scalar,scalar_flux, &
        MPI_COMM_WORLD,good,first,last,boundary_faces=faces)
      call require(good,'MPI scalar flux sync')
      do p=1,count
        call require(allocated(scalar%fields(p)%values).eqv.(owners(p)==rank),'owner-local scalar field')
        if(owners(p)/=rank) cycle
        call require(all(scalar%fields(p)%values==expected_scalar%fields(p)%values), &
          'exact serial/MPI scalar reflux field')
        call require(all(scalar_flux(p)%x==expected_scalar_flux(p)%x).and. &
          all(scalar_flux(p)%y==expected_scalar_flux(p)%y).and. &
          all(scalar_flux(p)%z==expected_scalar_flux(p)%z),'exact serial/MPI scalar reflux flux')
      enddo
      deallocate(expected_scalar%fields,scalar%fields,scalar_flux,expected_scalar_flux)
      do failure=1,merge(merge(5,4,allocated(faces)),merge(4,3,allocated(faces)),rotating)
        local=initial;lt=it
        do p=1,count
          if (owners(p)==rank) lf(p)=raw(p)
        end do
        if (owners(count)==rank) then
          if (failure==1) lf(count)%x(1,0,1,1)=ieee_value(1._dp,ieee_quiet_nan)
          if (failure==2) lf(count)%x(1,4,:,:)=1.e20_dp
        end if
        saved_flux=lf;first=0
        if (failure==3.and.rank==0) first=1
        if(allocated(faces)) then
          if(failure==4.and.rank==0) faces(2*mod(axis,3)+1:2*mod(axis,3)+2)='outflow'
        end if
        if(rotating) then
          if(((failure==4.and..not.allocated(faces)).or.failure==5).and.rank==0) frame%angular_speed=101._dp
          call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
            first,2,boundary_faces=faces,frame=frame)
          frame%angular_speed=100._dp
        else
          call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,geometry,local,lt,lf,MPI_COMM_WORLD,good, &
            first,2,boundary_faces=faces)
        end if
        call require(.not.good,'collective flux, EOS, level or periodic label rejection')
        kept_increment=increment;kept_replacement=replacement
        if(rotating) then
          ! Restore the rank-disagreeing frame for the assembly audit too.
          if(((failure==4.and..not.allocated(faces)).or.failure==5).and.rank==0) frame%angular_speed=101._dp
          call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,local,lt,lf,increment,replacement, &
            MPI_COMM_WORLD,good,first,2,boundary_faces=faces,frame=frame)
          frame%angular_speed=100._dp
        else
          call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,geometry,local,lt,lf,increment,replacement, &
            MPI_COMM_WORLD,good,first,2,boundary_faces=faces)
        endif
        if(failure==2) then
          call require(good,'assembly retains finite nonphysical reflux for composite reconciliation')
        else
          call require(.not.good,'assembly collectively rejects invalid flux or metadata')
          do p=1,count
            if(owners(p)/=rank) cycle
            call require(all(increment%fields(p)%values==kept_increment%fields(p)%values), &
              'failed assembly retains prior signed increments')
            call require(all(replacement(p)%x==kept_replacement(p)%x).and. &
              all(replacement(p)%y==kept_replacement(p)%y).and.all(replacement(p)%z==kept_replacement(p)%z), &
              'failed assembly retains prior replacement faces')
          enddo
        endif

        if(allocated(faces)) then
          faces='periodic';faces(2*axis-1:2*axis)='outflow'
        end if
        do p=1,count
          if (owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values),'U rollback')
          call require(all(lt%fields(p)%values==it%fields(p)%values),'T rollback')
          call require(all(transfer(lf(p)%x,[0_int64])==transfer(saved_flux(p)%x,[0_int64])).and. &
            all(transfer(lf(p)%y,[0_int64])==transfer(saved_flux(p)%y,[0_int64])).and. &
            all(transfer(lf(p)%z,[0_int64])==transfer(saved_flux(p)%z,[0_int64])),'bitwise flux rollback')
        end do
      end do
      deallocate(local%fields,lt%fields,lf)
    end do
  end subroutine check_mpi_sync
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_flux_3d
