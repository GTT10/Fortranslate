program test_mpi_eb_tree_hydro_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_is_finite
  use precision_mod, only: dp
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d
  use mpi_eb_tree_redistribution_3d_mod, only: compute_state_redist_mpi_level_stage_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  use eb_geometry_3d_mod, only: extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: fill_reactive_eb_tree_halo_3d
  use mpi_eb_tree_full_3d_mod, only: mpi_eb_mol_interval_3d,reconstruct_mpi_eb_mol_interval_3d, &
    advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d,replay_mpi_eb_mol_held_interval_3d
  use dual_energy_mod, only: dual_energy_policy
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use mpi_eb_tree_hydro_3d_mod, only: compute_mol_mpi_eb_tree_hydro_3d,fill_mpi_eb_flux_allowance_halo_3d, &
    apply_mpi_eb_flux_correction_3d
  use eb_reactive_hydro_3d_mod, only: limit_reactive_eb_flux_correction_3d, &
    correct_reactive_eb_spatial_hydro_flux_3d
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d,convert_rotating_energy_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_hydro_3d, &
    build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_hydro_3d_mod, only: advance_mpi_eb_tree_hydro_3d,compute_sdc_mpi_eb_tree_hydro_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d,repair_mpi_eb_tree_interface_budget_3d, &
    assemble_mpi_eb_tree_reflux_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use amr_patch_tree_remap_3d_mod, only: covered_amr=>covered
  use eb_reactive_patch_tree_3d_mod, only: compute_sdc_eb_tree_hydro_3d,compute_mol_eb_tree_hydro_3d
  use eb_reactive_patch_tree_3d_mod, only: fill_sdc_eb_tree_source_halo_3d
  use eb_reactive_patch_tree_3d_mod, only: apply_eb_tree_flux_correction_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_tree_source_halo_3d, &
    fill_sdc_mpi_eb_particle_predictor_halo_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,extended_finest
  type(eb_geometry_3d), allocatable :: extended(:),physical(:)
  type(eb_reactive_tree_3d) :: tree,saved
  type(rotating_frame_3d), allocatable :: frame
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),inlet(:,:),inlet_t(:)
  character(len=16), allocatable :: labels(:)
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),primitive(:)
  real(dp) :: residual,coefficient,et,sound,position(3),potential
  character(len=32) :: test_mode,solver
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,rank,nranks,status,i,j,k
  logical :: ok
  integer :: held_clock_reactor_calls=0,held_clock_refusals=0
  call get_command_argument(1,test_mode)
  solver='hllc'
  if(index(test_mode,'ctu')==1) solver='pelec'
  coefficient=0
  if(test_mode=='viscosity'.or.index(test_mode,'ctu')==1) coefficient=.037_dp
  if(test_mode=='default') coefficient=.1_dp
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),wall(nv),primitive(nv))
  if(test_mode=='flux_budget') then
    call check_spatial_flux_correction()
    call check_flux_allowance_exchange()
    if(rank==0) print *, 'PASS: MPI flux donor allowances, shared correction, periodic seams and rollback'
    call MPI_Finalize(status)
    stop
  endif
  do axis=1,3
  do count=4,5
  if((index(test_mode,'ctu')==1.or.test_mode=='mol'.or.test_mode=='mol_reflux').and.(axis/=3.or.count/=5)) cycle
  if((test_mode=='rotating_mol'.or.test_mode=='mol_halo').and.(axis/=3.or.count/=5)) cycle
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
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  if(test_mode=='rotating_mol') then
    allocate(frame)
    frame%axis=1;frame%angular_speed=100._dp;frame%axis_point=[.5_dp,.5_dp,.5_dp]
    do k=1,4
    do j=1,4
    do i=1,4
      if(coarse%volume_fraction(i,j,k)<=0) cycle
      call reactive_conserved_to_primitive(species,u(:,i,j,k),t(i,j,k),primitive,et,sound,ok)
      call require(ok,'rotating root primitive')
      primitive(5)=primitive(5)*1000._dp/et
      call reactive_primitive_to_conserved(species,primitive,u(:,i,j,k),t(i,j,k),sound,ok)
      call require(ok,'rotating root EOS')
      position=(real([i,j,k],dp)-.5_dp)/4
      call rotating_potential_3d(frame,position,potential,ok);call require(ok,'root potential')
      u(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
    end do
    end do
    end do
  end if
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok,frame)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  if(index(test_mode,'ctu')==1.or.test_mode=='rotating_mol'.or.test_mode=='mol'.or.test_mode=='mol_reflux') then
    allocate(physical(count))
    do node=1,count
      physical(node)=tree%nodes(node)%geometry
      if(test_mode=='rotating_mol') cycle
      do k=1,4
      do j=1,4
      do i=1,4
        if(physical(node)%volume_fraction(i,j,k)<=0) cycle
        call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k), &
          tree%nodes(node)%temperature(i,j,k),primitive,et,sound,ok)
        call require(ok,'CTU primitive fixture')
        primitive(5)=primitive(5)*1000._dp/et
        call reactive_primitive_to_conserved(species,primitive,tree%nodes(node)%state(:,i,j,k), &
          tree%nodes(node)%temperature(i,j,k),sound,ok)
        call require(ok,'CTU fixture EOS')
      end do
      end do
      end do
    end do
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      axes(axis),.31_dp,extended_finest,ok)
    call require(ok,'extended finest')
    call build_eb_tree_extended_geometry_3d(topology,physical,extended_finest,extended,ok)
    call require(ok,'extended cache')
  end if
  call average_down_reactive_eb_tree_3d(species,tree,ok,frame)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  if(test_mode=='mol_halo') then
    call check_mol_spatial_halos()
    if(rank==0) print *, 'PASS: MPI prepared hydro spatial halos and StateRedist stage'
    call MPI_Finalize(status)
    stop
  endif
  if(test_mode=='mol'.or.test_mode=='mol_reflux') then
    call check_mol_only()
    deallocate(u,t,physical,extended)
    cycle
  end if
  if(test_mode=='rotating_mol') then
    call check_rotating_mol()
    deallocate(u,t,frame,physical,extended)
    cycle
  end if
  if(index(test_mode,'inflow')>0) then
    allocate(labels(6),inlet(nv,6),inlet_t(6))
    labels='outflow';labels(1)='inflow';labels(4)='inflow'
    call reactive_conserved_to_primitive(species,tree%nodes(1)%state(:,1,1,3), &
      tree%nodes(1)%temperature(1,1,3),primitive,et,sound,ok)
    call require(ok,'inlet primitive')
    primitive(2:4)=[100._dp,-50._dp,0._dp];primitive(5)=.8_dp*primitive(5)
    do i=1,6
      call reactive_primitive_to_conserved(species,primitive,inlet(:,i),inlet_t(i),sound,ok)
      call require(ok,'inlet EOS')
    end do
  end if
  if(index(test_mode,'periodic')>0) then
    allocate(labels(6));labels='periodic';labels(5:6)='outflow'
  end if
  if(index(test_mode,'wall')>0) then
    allocate(labels(6));labels='outflow';labels(1:2)='slip_wall';labels(3:4)='no_slip_wall'
  end if
  if(allocated(extended)) call check_sdc_hydro()
  boundary=17._dp;wall=19._dp
  call advance_reactive_eb_tree_hydro_3d(species,solver,1.e-8_dp,tree,boundary,wall,ok, &
    difmag=coefficient,extended_geometry=extended, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
  call require(ok,'tree hydro step')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary-wall)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'hydro conservation including boundary and exposed wall')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep hydro update')
  call check_mpi_hydro()
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_hydro_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_mol_spatial_halos()
    type(eb_reactive_tree_3d) :: fixture
    type(eb_geometry_3d), allocatable :: geo(:),wide(:)
    type(eb_geometry_3d) :: padded
    type(amr_tree_state_3d) :: local,lt,rates,halos,old,serial_rates,serial_halos,serial_old, &
      kept_rates,kept_halos,kept_old,redist,serial_aux_rate,mpi_aux_rate,serial_aux_halos,mpi_aux_halos, &
      serial_held,local_held,serial_aux_old,mpi_aux_old,serial_held_end,local_held_end
    type(eb_integrated_flux_3d), allocatable :: fluxes(:),serial_fluxes(:),kept_fluxes(:), &
      serial_aux_flux(:),mpi_aux_flux(:)
    type(state_redist_register_3d), allocatable :: regs(:)
    type(rotating_frame_3d), allocatable :: spin
    real(dp), allocatable :: phi
    real(dp) :: q(reactive_nprim(size(species))),pos(3),temp_cell,speed
    integer :: policy_id, owners(count),p,x,y,z,c,ownership,rotation_id,depth,solver_id
    character(len=32) :: faces(6)
    logical :: good
    fixture=saved
    allocate(geo(count))
    do p=1,count
      geo(p)=fixture%nodes(p)%geometry
    enddo
    call build_axis_plane_eb_geometry_3d(64,64,64,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'z',.31_dp,padded,good)
    call require(good,'MPI hydro halo padded geometry')
    call build_eb_tree_extended_geometry_3d(topology,geo,padded,wide,good)
    call require(good,'MPI hydro halo prepared geometry')
    do rotation_id=0,1
      if(rotation_id==1) then
        allocate(spin,phi)
        spin%axis=3;spin%angular_speed=100;spin%axis_point=[.5_dp,.5_dp,.5_dp]
      endif
      do p=1,count
        do z=1,4
        do y=1,4
        do x=1,4
          if(geo(p)%volume_fraction(x,y,z)<=0) cycle
          pos=[geo(p)%x_lower+(x-.5_dp)*geo(p)%dx,geo(p)%y_lower+(y-.5_dp)*geo(p)%dy, &
            geo(p)%z_lower+(z-.5_dp)*geo(p)%dz]
          q=0;q(1)=1+.01_dp*pos(1);q(2:4)=.1_dp*pos;q(5)=135000
          do c=1,size(species)
            if(species(c)%name=='N2') q(c+5)=1
          enddo
          if(allocated(spin)) then
            call rotating_potential_3d(spin,pos,phi,good)
            call require(good,'MPI hydro halo initial potential')
          endif
          call reactive_primitive_to_conserved(species,q,fixture%nodes(p)%state(:,x,y,z),temp_cell,speed,good,phi)
          call require(good,'MPI hydro halo smooth state')
          fixture%nodes(p)%temperature(x,y,z)=temp_cell
        enddo
        enddo
        enddo
      enddo
      call average_down_reactive_eb_tree_3d(species,fixture,good,spin)
      call require(good,'MPI hydro halo restriction')
      allocate(serial_held%fields(count))
      do p=1,count
        allocate(serial_held%fields(p)%values(1,4,4,4))
        serial_held%fields(p)%values=0
        where(geo(p)%volume_fraction>0) &
          serial_held%fields(p)%values(1,:,:,:)=fixture%nodes(p)%state(5,:,:,:)+real(p,dp)
      enddo
      faces='outflow'
      if(rotation_id==1) faces(1:4)='periodic'
      do ownership=0,2
        owners=0
        if(ownership==1) owners=[(modulo(p-1,nranks),p=1,count)]
        if(ownership==2) owners=[(modulo(count-p,nranks),p=1,count)]
        allocate(local%fields(count),lt%fields(count),local_held%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          local%fields(p)%values=fixture%nodes(p)%state
          allocate(lt%fields(p)%values(1,4,4,4))
          lt%fields(p)%values(1,:,:,:)=fixture%nodes(p)%temperature
          local_held%fields(p)%values=serial_held%fields(p)%values
        enddo
        do solver_id=0,1
        do depth=-1,2
          if(depth==-1) then
            call compute_mol_eb_tree_hydro_3d(species,fixture,wide,1.e-9_dp,2,serial_rates,serial_fluxes,good, &
              boundary_faces=faces,lax_friedrichs=solver_id==1,frame=spin, &
              spatial_rate_halos=serial_halos,old_state_halos=serial_old)
          else
            call compute_mol_eb_tree_hydro_3d(species,fixture,wide,1.e-9_dp,2,serial_rates,serial_fluxes,good, &
              boundary_faces=faces,target_depth=depth,lax_friedrichs=solver_id==1,frame=spin, &
              spatial_rate_halos=serial_halos,old_state_halos=serial_old)
          endif
          call require(good,'MPI hydro halo serial control')
          if(depth==-1) then
            call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
              rates,fluxes,MPI_COMM_WORLD,good,boundary_faces=faces,lax_friedrichs=solver_id==1,frame=spin, &
              spatial_rate_halos=halos,old_state_halos=old)
          else
            call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
              rates,fluxes,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=depth, &
              lax_friedrichs=solver_id==1,frame=spin,spatial_rate_halos=halos,old_state_halos=old)
          endif
          call require(good,'MPI hydro halo collective computation')
          if(rotation_id==0.and.solver_id==1.and.depth==-1) then
            call compute_mol_eb_tree_hydro_3d(species,fixture,wide,1.e-9_dp,2,serial_rates,serial_fluxes,good, &
              boundary_faces=faces,lax_friedrichs=.true.,spatial_rate_halos=serial_halos, &
              old_state_halos=serial_old,internal_energy_rate=serial_aux_rate,internal_energy_flux=serial_aux_flux, &
              internal_energy_spatial_rate_halos=serial_aux_halos,internal_energy_state=serial_held, &
              internal_energy_old_state_halos=serial_aux_old)
            call require(good,'serial retained energy LF spatial control')
            call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
              rates,fluxes,MPI_COMM_WORLD,good,boundary_faces=faces,lax_friedrichs=.true., &
              spatial_rate_halos=halos,old_state_halos=old, &
              internal_energy_rate=mpi_aux_rate,internal_energy_flux=mpi_aux_flux, &
              internal_energy_spatial_rate_halos=mpi_aux_halos,internal_energy_state=local_held, &
              internal_energy_old_state_halos=mpi_aux_old)
            call require(good,'MPI retained energy LF spatial evaluation')
            do p=1,count
              if(owners(p)/=rank) then
                call require(.not.allocated(mpi_aux_rate%fields(p)%values), &
                  'MPI retained energy absent on remote owner')
                call require(.not.allocated(mpi_aux_halos%fields(p)%values), &
                  'MPI retained energy halo absent on remote owner')
                call require(.not.allocated(mpi_aux_old%fields(p)%values), &
                  'MPI retained old energy halo absent on remote owner')
                cycle
              endif
              call require(maxval(abs(mpi_aux_halos%fields(p)%values-serial_aux_halos%fields(p)%values)/ &
                max(1._dp,abs(serial_aux_halos%fields(p)%values)))<1.e-12_dp, &
                'MPI retained energy LF halo matches serial')
              call require(maxval(abs(mpi_aux_old%fields(p)%values-serial_aux_old%fields(p)%values)/ &
                max(1._dp,abs(serial_aux_old%fields(p)%values)))<1.e-12_dp, &
                'MPI retained energy old halo matches serial')
              call require(maxval(abs(mpi_aux_rate%fields(p)%values-serial_aux_rate%fields(p)%values)/ &
                max(1._dp,abs(serial_aux_rate%fields(p)%values)))<1.e-12_dp, &
                'MPI retained energy LF rate matches serial')
              call require(maxval(abs(mpi_aux_flux(p)%x-serial_aux_flux(p)%x)/ &
                max(1._dp,abs(serial_aux_flux(p)%x)))<1.e-12_dp.and. &
                maxval(abs(mpi_aux_flux(p)%y-serial_aux_flux(p)%y)/ &
                max(1._dp,abs(serial_aux_flux(p)%y)))<1.e-12_dp.and. &
                maxval(abs(mpi_aux_flux(p)%z-serial_aux_flux(p)%z)/ &
                max(1._dp,abs(serial_aux_flux(p)%z)))<1.e-12_dp, &
                'MPI retained energy LF flux matches serial')
            enddo
          endif
          do p=1,count
            if(owners(p)/=rank) then
              call require(.not.allocated(halos%fields(p)%values).and..not.allocated(old%fields(p)%values), &
                'MPI hydro halo no remote state/source allocation')
              cycle
            endif
            call require(maxval(abs(rates%fields(p)%values-serial_rates%fields(p)%values)/ &
              max(1._dp,abs(serial_rates%fields(p)%values)))<1.e-12_dp,'MPI hydro halo valid serial equality')
            if(depth>=0.and.topology%nodes(p)%level/=depth) then
              call require(.not.allocated(halos%fields(p)%values),'MPI hydro halo inactive level absent')
              cycle
            endif
            call require(maxval(abs(fluxes(p)%x-serial_fluxes(p)%x)/max(1._dp,abs(serial_fluxes(p)%x)))<1.e-12_dp.and. &
              maxval(abs(fluxes(p)%y-serial_fluxes(p)%y)/max(1._dp,abs(serial_fluxes(p)%y)))<1.e-12_dp.and. &
              maxval(abs(fluxes(p)%z-serial_fluxes(p)%z)/max(1._dp,abs(serial_fluxes(p)%z)))<1.e-12_dp, &
              'MPI hydro halo integrated flux serial equality')
            call require(all(lbound(old%fields(p)%values)==1),'MPI hydro halo state indices')
            call require(maxval(abs(old%fields(p)%values-serial_old%fields(p)%values)/ &
              max(1._dp,abs(serial_old%fields(p)%values)))<1.e-12_dp,'MPI hydro halo old serial equality')
            call require(maxval(abs(halos%fields(p)%values-serial_halos%fields(p)%values)/ &
              max(1._dp,abs(serial_halos%fields(p)%values)))<1.e-12_dp,'MPI hydro halo entire source serial equality')
            call require(all(halos%fields(p)%values(:,7:10,7:10,7:10)==rates%fields(p)%values), &
              'MPI hydro halo valid source after peer limiting')
          enddo
          if(depth>=0) then
            redist=rates
            call compute_state_redist_mpi_level_stage_3d(topology,owners,wide,old,halos,nv,depth,1.e-9_dp, &
              redist,regs,MPI_COMM_WORLD,good,boundary_faces=faces,max_order=2)
            call require(good,'MPI hydro halo feeds StateRedist level')
          endif
        enddo
        enddo
        if(rotation_id==0.and.ownership==1) then
          serial_held_end=serial_held;local_held_end=local_held
          do p=1,count
            if(topology%nodes(p)%level==2) cycle
            serial_held_end%fields(p)%values=serial_held_end%fields(p)%values+1._dp
            if(owners(p)==rank) local_held_end%fields(p)%values=local_held_end%fields(p)%values+1._dp
          enddo
          do policy_id=0,1
          call compute_mol_eb_tree_hydro_3d(species,fixture,wide,1.e-9_dp,2,serial_rates,serial_fluxes,good, &
            boundary_faces=faces,target_depth=2,ancestor_start=fixture,ancestor_end=fixture, &
            ancestor_fraction=[.5_dp,.5_dp],internal_energy_rate=serial_aux_rate, &
            internal_energy_flux=serial_aux_flux,internal_energy_state=serial_held, &
            internal_energy_old_state_halos=serial_aux_old, &
            ancestor_internal_energy_start=serial_held,ancestor_internal_energy_end=serial_held_end, &
            separate_momentum=policy_id==1)
          call require(good,'serial retained energy coarse-time halo control')
          call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
            rates,fluxes,MPI_COMM_WORLD,good,boundary_faces=faces,target_depth=2, &
            ancestor_start=local,ancestor_end=local,temperature_start=lt,temperature_end=lt, &
            ancestor_fraction=[.5_dp,.5_dp],internal_energy_rate=mpi_aux_rate, &
            internal_energy_flux=mpi_aux_flux,internal_energy_state=local_held, &
            internal_energy_old_state_halos=mpi_aux_old, &
            ancestor_internal_energy_start=local_held,ancestor_internal_energy_end=local_held_end, &
            separate_momentum=policy_id==1)
          call require(good,'MPI retained energy coarse-time halo')
          do p=1,count
            if(owners(p)/=rank.or.topology%nodes(p)%level/=2) then
              call require(.not.allocated(mpi_aux_old%fields(p)%values), &
                'MPI retained energy inactive coarse-time halo absent')
              cycle
            endif
            call require(maxval(abs(mpi_aux_old%fields(p)%values-serial_aux_old%fields(p)%values)/ &
              max(1._dp,abs(serial_aux_old%fields(p)%values)))<1.e-12_dp, &
              'MPI retained energy coarse-time halo matches serial')
          enddo
          enddo
        endif
        kept_rates=rates;kept_halos=halos;kept_old=old;kept_fluxes=fluxes
        if(nranks>1) then
        call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
          rates,fluxes,MPI_COMM_WORLD,good,internal_energy_rate=mpi_aux_rate, &
          internal_energy_flux=mpi_aux_flux,internal_energy_state=local_held, &
          internal_energy_old_state_halos=mpi_aux_old,separate_momentum=rank==0)
        call require(.not.good,'MPI hydro halo rejects mismatched momentum policy')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(rates%fields(p)%values==kept_rates%fields(p)%values), &
            'MPI hydro momentum policy refusal preserves rates')
          if(allocated(fluxes(p)%x)) then
            call require(all(fluxes(p)%x==kept_fluxes(p)%x).and.all(fluxes(p)%y==kept_fluxes(p)%y).and. &
              all(fluxes(p)%z==kept_fluxes(p)%z),'MPI hydro momentum policy refusal preserves fluxes')
          endif
        enddo
        endif
        if(owners(count)==rank) local%fields(count)%values(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
        call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
          rates,fluxes,MPI_COMM_WORLD,good,boundary_faces=faces,frame=spin,spatial_rate_halos=halos,old_state_halos=old)
        call require(.not.good,'MPI hydro halo remote EOS failure')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(rates%fields(p)%values==kept_rates%fields(p)%values),'MPI hydro halo rates rollback')
          call require(allocated(fluxes(p)%x).eqv.allocated(kept_fluxes(p)%x),'MPI hydro halo flux allocation rollback')
          if(allocated(fluxes(p)%x)) then
            call require(all(fluxes(p)%x==kept_fluxes(p)%x).and.all(fluxes(p)%y==kept_fluxes(p)%y).and. &
              all(fluxes(p)%z==kept_fluxes(p)%z),'MPI hydro halo flux rollback')
          endif
          if(allocated(old%fields(p)%values)) then
            call require(all(old%fields(p)%values==kept_old%fields(p)%values),'MPI hydro halo snapshot rollback')
            call require(all(halos%fields(p)%values==kept_halos%fields(p)%values),'MPI hydro halo RHS rollback')
          endif
        enddo
        if(owners(count)==rank) local%fields(count)%values=fixture%nodes(count)%state
        if(nranks>1) then
          if(rank==0) then
            call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
              rates,fluxes,MPI_COMM_WORLD,good,spatial_rate_halos=halos)
          else
            call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,geo,wide,local,lt,1.e-9_dp,2, &
              rates,fluxes,MPI_COMM_WORLD,good)
          endif
          call require(.not.good,'MPI hydro halo mismatched output presence')
        endif
        deallocate(local%fields,lt%fields,local_held%fields)
      enddo
      deallocate(serial_held%fields)
    enddo
  end subroutine check_mol_spatial_halos

  subroutine check_spatial_flux_correction()
    type(eb_geometry_3d) :: g
    type(eb_integrated_flux_3d) :: baseline,flux,saved_flux
    real(dp) :: cells(nv,4,4,4),source(nv,4,4,4),rate(nv,4,4,4),saved_rate(nv,4,4,4), &
      f(nv),residual(nv),dt,factor,total(nv),endpoint(nv)
    integer :: c
    logical :: good
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,4._dp,0._dp,4._dp,0._dp,4._dp,'x',-.1_dp,g,good)
    call require(good,'spatial correction regular geometry')
    allocate(baseline%x(nv,0:4,4,4),baseline%y(nv,4,0:4,4),baseline%z(nv,4,4,0:4))
    baseline%x=0;baseline%y=0;baseline%z=0;flux=baseline
    cells=0;cells(1,:,:,:)=1;cells(6,:,:,:)=.1_dp;cells(7,:,:,:)=.9_dp
    source=0;source(6,:,:,:)=.1_dp;source(7,:,:,:)=-.1_dp
    dt=.1_dp;f=0;f(2:5)=[2._dp,-3._dp,4._dp,7._dp];f(6)=2;f(7)=-2
    residual=0;residual(2:5)=[3._dp,-2._dp,1._dp,5._dp]
    do c=1,nv
      rate(c,:,:,:)=residual(c)
    enddo
    rate(:,2,2,2)=rate(:,2,2,2)-f;rate(:,3,2,2)=rate(:,3,2,2)+f
    flux%x(:,2,2,2)=dt*f
    call correct_reactive_eb_spatial_hydro_flux_3d(g,cells,source,dt,rate,flux,baseline,good)
    call require(good,'whole spatial CTU correction accepted')
    factor=((.1_dp+dt*.1_dp)/(dt*2))*(1-16*epsilon(1._dp))
    call require(maxval(abs(flux%x(:,2,2,2)-dt*factor*f))<1.e-15_dp,'source-aware common face factor')
    call require(maxval(abs(rate(:,2,2,2)-(residual-factor*f)))<1.e-14_dp.and. &
      maxval(abs(rate(:,3,2,2)-(residual+factor*f)))<1.e-14_dp,'spatial rate retains wall residual')
    endpoint=cells(:,2,2,2)+dt*(rate(:,2,2,2)+source(:,2,2,2))
    call require(endpoint(1)>0.and.all(endpoint(6:)>=0),'corrected spatial donor remains admissible')
    do c=1,nv
      total(c)=sum(rate(c,:,:,:))
    enddo
    call require(maxval(abs(total-64*residual))<1.e-13_dp,'internal spatial face preserves every inventory')
    saved_rate=rate;saved_flux=flux;source(6,:,:,:)=-2
    call correct_reactive_eb_spatial_hydro_flux_3d(g,cells,source,dt,rate,flux,baseline,good)
    call require(.not.good,'inadmissible low-order source rejected')
    call require(all(rate==saved_rate).and.all(flux%x==saved_flux%x).and.all(flux%y==saved_flux%y).and. &
      all(flux%z==saved_flux%z),'spatial flux/rate joint rollback')
    ! Extracted recursive-spray failure: removal of the whole baseline inflow
    ! hides a much smaller orthogonal outgoing face in the summed allowance.
    ! Both flows remain represented when the limited endpoint is reconstructed.
    baseline%x=0;baseline%y=0;baseline%z=0
    cells=0;cells(1,:,:,:)=1;cells(nv,:,:,:)=1
    dt=9.3847656249999995e-9_dp
    baseline%y(6,2,0,1)=6.4911469960607036e-67_dp
    baseline%y(nv,2,0,1)=-baseline%y(6,2,0,1)
    flux=baseline;flux%y(:,2,0,1)=0
    flux%z(6,2,1,0)=-4.3133348252984133e-124_dp
    flux%z(nv,2,1,0)=-flux%z(6,2,1,0)
    call require(dt*baseline%y(6,2,0,1)+dt*abs(flux%z(6,2,1,0))== &
      dt*baseline%y(6,2,0,1),'small outgoing face is below summed allowance resolution')
    call limit_reactive_eb_flux_correction_3d(g,cells,baseline,flux,dt,good)
    call require(good,'saturated allowance retains orthogonal trace outflow')
    endpoint=cells(:,2,1,1)+dt*(flux%y(:,2,0,1)+flux%z(:,2,1,0))
    call require(all(endpoint(6:)>=0).and.endpoint(6)>0,'saturated trace endpoint is physical without clipping')
    call require(flux%z(6,2,1,0)<0,'represented trace outflow is retained')
    call require(flux%y(6,2,0,1)+flux%y(nv,2,0,1)==0.and. &
      flux%z(6,2,1,0)+flux%z(nv,2,1,0)==0,'common allowance preserves paired species fluxes')
  end subroutine check_spatial_flux_correction
  subroutine check_flux_allowance_exchange()
    type(amr_patch_3d) :: siblings(2)
    type(amr_patch_tree_3d) :: layout
    type(amr_tree_state_3d) :: allowance,halo,saved_halo
    type(amr_tree_state_3d) :: patch_state,patch_source,patch_rate,saved_rate
    type(eb_reactive_tree_3d) :: serial_tree
    type(amr_tree_state_3d) :: serial_source,serial_rate
    type(eb_integrated_flux_3d), allocatable :: serial_low(:),serial_faces(:)
    type(eb_geometry_3d) :: g
    type(eb_geometry_3d) :: metrics(3)
    type(eb_integrated_flux_3d) :: baseline,raw,limited
    type(eb_integrated_flux_3d), allocatable :: low_faces(:),high_faces(:),saved_faces(:)
    type(amr_tree_state_3d) :: held_rate,serial_held_rate,saved_held_rate
    type(eb_integrated_flux_3d), allocatable :: held_low(:),held_faces(:),serial_held_low(:),serial_held_faces(:)
    type(eb_integrated_flux_3d), allocatable :: saved_held_faces(:)
    type(mpi_eb_mol_interval_3d) :: held_history
    type(amr_tree_state_3d) :: reconstructed_gas,reconstructed_held,kept_reconstructed_gas,kept_reconstructed_held
    type(eb_integrated_flux_3d), allocatable :: changed_gas_faces(:),changed_held_faces(:)
    real(dp) :: gas_face_change,held_face_change,expected_held
    real(dp) :: held_change,held_total,raw_held_flux
    real(dp) :: cells(nv,4,4,4),f(nv),local_min(nv),local_max(nv),minimum(nv),maximum(nv)
    real(dp) :: change(nv),total(nv),endpoint(nv),dt,theta
    integer :: ownership(3),p,n(3),face,cell_index
    logical :: good,wrap(3)
    call require(nranks>=2,'flux allowance test needs two or more ranks')
    call initialize_amr_patch_3d(4,2,2,1,2,1,2,1,2,2,siblings(1),good)
    call require(good,'allowance left patch')
    call initialize_amr_patch_3d(4,2,2,3,4,1,2,1,2,2,siblings(2),good)
    call require(good,'allowance right patch')
    call initialize_amr_patch_tree_3d([4,2,2],[1,1],siblings,layout,good)
    call require(good,'allowance tree')
    ownership=[0,0,1]
    allocate(allowance%fields(3))
    if(rank==0) allocate(allowance%fields(1)%values(nv-4,4,2,2),source=1._dp)
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,4._dp,0._dp,4._dp,0._dp,4._dp,'x',-.1_dp,g,good)
    call require(good,'allowance donor geometry')
    allocate(baseline%x(nv,0:4,4,4),baseline%y(nv,4,0:4,4),baseline%z(nv,4,4,0:4))
    baseline%x=0;baseline%y=0;baseline%z=0
    dt=.1_dp;f=0;f(2:5)=[2._dp,-3._dp,4._dp,7._dp];f(6)=1;f(7)=-1
    cells=0;cells(1,:,:,:)=1
    if(rank==0.or.rank==1) then
      p=rank+2;face=merge(4,0,rank==0);cell_index=merge(4,1,rank==0)
      if(rank==0) then
        cells(6,:,:,:)=.25_dp;cells(7,:,:,:)=.75_dp
      else
        cells(6,:,:,:)=.96_dp;cells(7,:,:,:)=.04_dp
      endif
      raw=baseline;raw%x(:,face,2,2)=f;limited=raw
      call limit_reactive_eb_flux_correction_3d(g,cells,baseline,limited,dt,good, &
        alpha_output=allowance%fields(p)%values)
      call require(good,'owner-local donor allowance')
    endif
    wrap=.false.
    call fill_mpi_eb_flux_allowance_halo_3d(layout,ownership,allowance,nv-4,halo,MPI_COMM_WORLD,good,wrap)
    call require(good,'distributed donor exchange including idle ranks')
    local_min=huge(1._dp);local_max=-huge(1._dp);change=0
    if(rank==0.or.rank==1) then
      call require(all(lbound(halo%fields(p)%values)==[1,0,0,0]),'donor halo face indices')
      limited=raw
      call limit_reactive_eb_flux_correction_3d(g,cells,baseline,limited,dt,good, &
        peer_alpha=halo%fields(p)%values)
      call require(good,'shared correction consumes MPI donor allowance')
      local_min=limited%x(:,face,2,2);local_max=local_min
      theta=(.04_dp/dt)*(1-16*epsilon(1._dp))
      call require(maxval(abs(local_min-theta*f))<1e-15_dp,'remote loss controls the common face factor')
      change=merge(-dt,dt,rank==0)*local_min
      endpoint=cells(:,cell_index,2,2)+change
      call require(all(endpoint(6:)>=0),'both distributed donors retain nonnegative species')
    endif
    call MPI_Allreduce(local_min,minimum,nv,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS,'shared minimum reduction')
    call MPI_Allreduce(local_max,maximum,nv,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.all(minimum==maximum),'opposite MPI faces match exactly')
    call MPI_Allreduce(change,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.all(total==0),'distributed correction preserves every inventory')
    if(rank==0) then
      allowance%fields(2)%values(1,1,1,1)=.321_dp
      allowance%fields(2)%values(2,2,4,2)=.271_dp
    endif
    if(rank==1) allowance%fields(3)%values(1,4,1,1)=.123_dp
    wrap=.true.
    call fill_mpi_eb_flux_allowance_halo_3d(layout,ownership,allowance,nv-4,halo,MPI_COMM_WORLD,good,wrap)
    call require(good,'periodic donor exchange')
    if(rank==0) then
      call require(halo%fields(2)%values(1,0,1,1)==.123_dp,'periodic lower peer donor')
      call require(halo%fields(2)%values(2,2,0,2)==.271_dp,'periodic same-owner donor')
    endif
    if(rank==1) call require(halo%fields(3)%values(1,5,1,1)==.321_dp,'periodic upper peer donor')
    saved_halo=halo
    if(rank==1) allowance%fields(3)%values(3,1,2,2)=-1
    call fill_mpi_eb_flux_allowance_halo_3d(layout,ownership,allowance,nv-4,halo,MPI_COMM_WORLD,good,wrap)
    call require(.not.good,'invalid remote allowance is collectively rejected')
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(halo%fields(p)%values==saved_halo%fields(p)%values),'failed exchange preserves existing halo')
    enddo
    if(rank==1) allowance%fields(3)%values(3,1,2,2)=1
    wrap=rank==0
    call fill_mpi_eb_flux_allowance_halo_3d(layout,ownership,allowance,nv-4,halo,MPI_COMM_WORLD,good,wrap)
    call require(.not.good,'rank-disagreeing periodic allowance policy rejected')
    do p=1,3
      if(ownership(p)/=rank) cycle
      n=layout%nodes(p)%shape
      call require(all(shape(halo%fields(p)%values)==[nv-4,n+2]).and. &
        all(halo%fields(p)%values==saved_halo%fields(p)%values),'policy rejection preserves halo payload and bounds')
    enddo
    call build_axis_plane_eb_geometry_3d(4,2,2,0._dp,8._dp,0._dp,4._dp,0._dp,4._dp,'x',-.1_dp,metrics(1),good)
    call require(good,'rate correction root metrics')
    metrics(2)=g
    call build_axis_plane_eb_geometry_3d(4,4,4,4._dp,8._dp,0._dp,4._dp,0._dp,4._dp,'x',-.1_dp,metrics(3),good)
    call require(good,'rate correction right metrics')
    allocate(patch_state%fields(3),patch_source%fields(3),patch_rate%fields(3),low_faces(3),high_faces(3))
    do p=1,3
      if(ownership(p)/=rank) cycle
      n=layout%nodes(p)%shape
      allocate(patch_state%fields(p)%values(nv,n(1),n(2),n(3)), &
        patch_source%fields(p)%values(nv,n(1),n(2),n(3)),patch_rate%fields(p)%values(nv,n(1),n(2),n(3)))
      patch_state%fields(p)%values=0;patch_state%fields(p)%values(1,:,:,:)=1
      patch_state%fields(p)%values(7,:,:,:)=1
      patch_source%fields(p)%values=0;patch_rate%fields(p)%values=0
      allocate(low_faces(p)%x(nv,0:n(1),n(2),n(3)),low_faces(p)%y(nv,n(1),0:n(2),n(3)), &
        low_faces(p)%z(nv,n(1),n(2),0:n(3)))
      low_faces(p)%x=0;low_faces(p)%y=0;low_faces(p)%z=0;high_faces(p)=low_faces(p)
      if(p==1) cycle
      patch_state%fields(p)%values=cells
      face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
      high_faces(p)%x(:,face,2,2)=dt*f
      patch_rate%fields(p)%values(:,cell_index,2,2)=merge(-1._dp,1._dp,p==2)*f
    enddo
    wrap=.false.
    serial_tree%topology=layout
    allocate(serial_tree%nodes(3),serial_source%fields(3),serial_rate%fields(3),serial_low(3),serial_faces(3))
    do p=1,3
      n=layout%nodes(p)%shape;serial_tree%nodes(p)%geometry=metrics(p)
      allocate(serial_tree%nodes(p)%state(nv,n(1),n(2),n(3)),serial_tree%nodes(p)%temperature(n(1),n(2),n(3)), &
        serial_source%fields(p)%values(nv,n(1),n(2),n(3)),serial_rate%fields(p)%values(nv,n(1),n(2),n(3)))
      serial_tree%nodes(p)%state=0;serial_tree%nodes(p)%state(1,:,:,:)=1
      serial_tree%nodes(p)%state(7,:,:,:)=1;serial_tree%nodes(p)%temperature=300
      serial_source%fields(p)%values=0;serial_rate%fields(p)%values=0
      allocate(serial_low(p)%x(nv,0:n(1),n(2),n(3)),serial_low(p)%y(nv,n(1),0:n(2),n(3)), &
        serial_low(p)%z(nv,n(1),n(2),0:n(3)))
      serial_low(p)%x=0;serial_low(p)%y=0;serial_low(p)%z=0;serial_faces(p)=serial_low(p)
      if(p==1) cycle
      serial_tree%nodes(p)%state(6,:,:,:)=merge(.25_dp,.96_dp,p==2)
      serial_tree%nodes(p)%state(7,:,:,:)=merge(.75_dp,.04_dp,p==2)
      face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
      serial_faces(p)%x(:,face,2,2)=dt*f
      serial_rate%fields(p)%values(:,cell_index,2,2)=merge(-1._dp,1._dp,p==2)*f
    enddo
    call apply_eb_tree_flux_correction_3d(serial_tree,serial_source,serial_low,dt,serial_rate,serial_faces,good, &
      wrap,target_depth=1)
    call require(good,'serial complete forward correction')
    call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
      patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1)
    call require(good,'MPI correction publishes matching rates and integrated flux')
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(patch_rate%fields(p)%values==serial_rate%fields(p)%values).and. &
        all(high_faces(p)%x==serial_faces(p)%x).and.all(high_faces(p)%y==serial_faces(p)%y).and. &
        all(high_faces(p)%z==serial_faces(p)%z),'serial/MPI rate and flux correction match exactly')
    enddo
    change=0
    do p=2,3
      if(ownership(p)/=rank) cycle
      face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
      endpoint=patch_state%fields(p)%values(:,cell_index,2,2)+dt*patch_rate%fields(p)%values(:,cell_index,2,2)
      call require(all(endpoint(6:)>=0),'corrected forward rate retains nonnegative species')
      call require(maxval(abs(dt*patch_rate%fields(p)%values(:,cell_index,2,2)- &
        merge(-1._dp,1._dp,p==2)*high_faces(p)%x(:,face,2,2)))<1e-16_dp, &
        'corrected rate equals the integrated face divergence')
      change=change+dt*patch_rate%fields(p)%values(:,cell_index,2,2)
    enddo
    call MPI_Allreduce(change,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.maxval(abs(total))<1e-16_dp,'corrected MPI rates conserve every inventory')
    saved_rate=patch_rate;saved_faces=high_faces
    call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
      patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1)
    call require(good,'already admissible correction can be reapplied')
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(patch_rate%fields(p)%values==saved_rate%fields(p)%values).and. &
        all(high_faces(p)%x==saved_faces(p)%x).and.all(high_faces(p)%y==saved_faces(p)%y).and. &
        all(high_faces(p)%z==saved_faces(p)%z),'unchanged correction is bitwise identity')
    enddo
    if(rank==1) patch_source%fields(3)%values(6,1,2,2)=-100
    call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
      patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1)
    call require(.not.good,'invalid complete source budget rejects correction collectively')
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(patch_rate%fields(p)%values==saved_rate%fields(p)%values).and. &
        all(high_faces(p)%x==saved_faces(p)%x).and.all(high_faces(p)%y==saved_faces(p)%y).and. &
        all(high_faces(p)%z==saved_faces(p)%z),'failed correction preserves rates and integrated flux together')
    enddo
    ! The held energy is a distinct signed face family, not a total-energy
    ! projection. Its spatial rate includes an independent volume term that
    ! must survive the correction of the shared species-limited face.
    allocate(held_rate%fields(3),serial_held_rate%fields(3),held_low(3),held_faces(3), &
      serial_held_low(3),serial_held_faces(3))
    raw_held_flux=-17.125_dp
    do p=1,3
      n=layout%nodes(p)%shape
      allocate(serial_held_rate%fields(p)%values(1,n(1),n(2),n(3)),source=3._dp)
      allocate(serial_held_low(p)%x(1,0:n(1),n(2),n(3)), &
        serial_held_low(p)%y(1,n(1),0:n(2),n(3)),serial_held_low(p)%z(1,n(1),n(2),0:n(3)))
      serial_held_low(p)%x=0;serial_held_low(p)%y=0;serial_held_low(p)%z=0
      serial_held_faces(p)=serial_held_low(p)
      if(p>1) then
        face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
        serial_faces(p)=serial_low(p);serial_faces(p)%x(:,face,2,2)=dt*f
        serial_rate%fields(p)%values=0
        serial_rate%fields(p)%values(:,cell_index,2,2)=merge(-1._dp,1._dp,p==2)*f
        serial_held_faces(p)%x(1,face,2,2)=dt*raw_held_flux
        serial_held_rate%fields(p)%values(1,cell_index,2,2)=3+merge(-1._dp,1._dp,p==2)*raw_held_flux
      endif
      if(ownership(p)/=rank) cycle
      held_rate%fields(p)=serial_held_rate%fields(p)
      held_low(p)=serial_held_low(p);held_faces(p)=serial_held_faces(p)
      patch_source%fields(p)%values=0
      high_faces(p)=serial_faces(p);patch_rate%fields(p)=serial_rate%fields(p)
    enddo
    call apply_eb_tree_flux_correction_3d(serial_tree,serial_source,serial_low,dt,serial_rate,serial_faces,good, &
      wrap,target_depth=1,internal_energy_baseline=serial_held_low,internal_energy_rate=serial_held_rate, &
      internal_energy_flux=serial_held_faces)
    call require(good,'serial held face follows complete gas correction')
    call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
      patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1,internal_energy_baseline=held_low, &
      internal_energy_rate=held_rate,internal_energy_flux=held_faces)
    call require(good,'MPI held face follows collective complete gas correction')
    theta=(.04_dp/dt)*(1-16*epsilon(1._dp));held_change=0
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(patch_rate%fields(p)%values==serial_rate%fields(p)%values).and. &
        all(high_faces(p)%x==serial_faces(p)%x).and.all(high_faces(p)%y==serial_faces(p)%y).and. &
        all(high_faces(p)%z==serial_faces(p)%z),'paired limiter preserves exact gas serial/MPI result')
      call require(all(held_rate%fields(p)%values==serial_held_rate%fields(p)%values).and. &
        all(held_faces(p)%x==serial_held_faces(p)%x).and.all(held_faces(p)%y==serial_held_faces(p)%y).and. &
        all(held_faces(p)%z==serial_held_faces(p)%z),'paired limiter exact owner/idle held rates and faces')
      if(p==1) then
        call require(all(held_rate%fields(p)%values==3),'paired limiter preserves inactive scalar level')
        cycle
      endif
      face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
      call require(abs(held_faces(p)%x(1,face,2,2)-theta*dt*raw_held_flux)<1.e-15_dp, &
        'signed held face uses independently calculated gas donor coefficient')
      call require(abs(dt*(held_rate%fields(p)%values(1,cell_index,2,2)-3)- &
        merge(-1._dp,1._dp,p==2)*held_faces(p)%x(1,face,2,2))<1.e-15_dp, &
        'held rate keeps volume term and matches corrected face divergence')
      held_change=held_change+dt*(held_rate%fields(p)%values(1,cell_index,2,2)-3)
    enddo
    call MPI_Allreduce(held_change,held_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.abs(held_total)<1.e-15_dp,'paired held correction conserves face inventory')
    saved_rate=patch_rate;saved_faces=high_faces;saved_held_rate=held_rate;saved_held_faces=held_faces
    if(rank==1) held_low(3)%x(1,0,2,2)=ieee_value(0._dp,ieee_quiet_nan)
    call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
      patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1,internal_energy_baseline=held_low, &
      internal_energy_rate=held_rate,internal_energy_flux=held_faces)
    call require(.not.good,'remote invalid held baseline rejects entire gas/scalar transaction')
    if(rank==1) held_low(3)%x(1,0,2,2)=0
    if(rank==nranks-1) then
      call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
        patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1)
    else
      call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
        patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=1,internal_energy_baseline=held_low, &
        internal_energy_rate=held_rate,internal_energy_flux=held_faces)
    endif
    call require(.not.good,'owner/idle held-option mismatch rejects collectively')
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(patch_rate%fields(p)%values==saved_rate%fields(p)%values).and. &
        all(high_faces(p)%x==saved_faces(p)%x).and.all(high_faces(p)%y==saved_faces(p)%y).and. &
        all(high_faces(p)%z==saved_faces(p)%z),'paired rejection keeps gas rates and faces')
      call require(all(held_rate%fields(p)%values==saved_held_rate%fields(p)%values).and. &
        all(held_faces(p)%x==saved_held_faces(p)%x).and.all(held_faces(p)%y==saved_held_faces(p)%y).and. &
        all(held_faces(p)%z==saved_held_faces(p)%z),'paired rejection keeps held rates and faces')
    enddo
    ! A manufactured accepted interval retains independent scalar volume
    ! terms. Replacing a shared gas/held face has an analytic signed change
    ! on its two cells and must not alter the archived reset ledger.
    held_history%depth=1;held_history%interval=dt;held_history%has_internal_energy=.true.
    held_history%initial_energy_correction=1.25_dp;held_history%numerical_energy_correction=2.5_dp
    allocate(held_history%forcing%fields(3),held_history%spatial_flux(3), &
      held_history%initial_internal_energy%fields(3),held_history%local_internal_energy%fields(3), &
      held_history%old_internal_energy_spatial_rate%fields(3),held_history%internal_energy_spatial_rate%fields(3), &
      held_history%old_internal_energy_flux(3),held_history%internal_energy_flux(3))
    do p=2,3
      if(ownership(p)/=rank) cycle
      held_history%forcing%fields(p)=patch_rate%fields(p);held_history%spatial_flux(p)=high_faces(p)
      held_history%initial_internal_energy%fields(p)=held_rate%fields(p)
      held_history%initial_internal_energy%fields(p)%values=7
      held_history%local_internal_energy%fields(p)=held_rate%fields(p)
      held_history%local_internal_energy%fields(p)%values=8
      held_history%old_internal_energy_spatial_rate%fields(p)=held_rate%fields(p)
      held_history%internal_energy_spatial_rate%fields(p)=held_rate%fields(p)
      held_history%old_internal_energy_flux(p)=held_faces(p);held_history%internal_energy_flux(p)=held_faces(p)
    enddo
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,held_history,high_faces,nv, &
      reconstructed_gas,MPI_COMM_WORLD,good,internal_energy_flux=held_faces,internal_energy_rate=reconstructed_held)
    call require(good,'paired interval identity reconstruction')
    do p=1,3
      if(ownership(p)/=rank) cycle
      if(p==1) then
        call require(all(reconstructed_gas%fields(p)%values==0).and.all(reconstructed_held%fields(p)%values==0), &
          'paired reconstruction keeps inactive source levels zero')
      else
        call require(all(reconstructed_gas%fields(p)%values==patch_rate%fields(p)%values).and. &
          all(reconstructed_held%fields(p)%values==held_rate%fields(p)%values), &
          'paired archive retains exact admitted spatial rate and volume term')
      endif
    enddo
    changed_gas_faces=high_faces;changed_held_faces=held_faces
    do p=2,3
      if(ownership(p)/=rank) cycle
      face=merge(4,0,p==2)
      changed_gas_faces(p)%x(5,face,2,2)=changed_gas_faces(p)%x(5,face,2,2)+.25_dp*dt
      changed_held_faces(p)%x(1,face,2,2)=changed_held_faces(p)%x(1,face,2,2)-.5_dp*dt
    enddo
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,held_history,changed_gas_faces,nv, &
      reconstructed_gas,MPI_COMM_WORLD,good,internal_energy_flux=changed_held_faces, &
      internal_energy_rate=reconstructed_held)
    call require(good,'paired interval changed-face reconstruction')
    held_change=0;change=0
    do p=2,3
      if(ownership(p)/=rank) cycle
      face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
      gas_face_change=(changed_gas_faces(p)%x(5,face,2,2)-high_faces(p)%x(5,face,2,2))/dt
      held_face_change=(changed_held_faces(p)%x(1,face,2,2)-held_faces(p)%x(1,face,2,2))/dt
      endpoint=patch_rate%fields(p)%values(:,cell_index,2,2)
      endpoint(5)=endpoint(5)+merge(-1._dp,1._dp,p==2)*gas_face_change
      expected_held=held_rate%fields(p)%values(1,cell_index,2,2)+merge(-1._dp,1._dp,p==2)*held_face_change
      call require(all(reconstructed_gas%fields(p)%values(:,cell_index,2,2)==endpoint).and. &
        reconstructed_held%fields(p)%values(1,cell_index,2,2)==expected_held, &
        'paired changed face has exact independently signed gas/scalar cell increment')
      change=change+dt*(reconstructed_gas%fields(p)%values(:,cell_index,2,2)- &
        patch_rate%fields(p)%values(:,cell_index,2,2))
      held_change=held_change+dt*(reconstructed_held%fields(p)%values(1,cell_index,2,2)- &
        held_rate%fields(p)%values(1,cell_index,2,2))
    enddo
    call MPI_Allreduce(change,total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.maxval(abs(total))<1.e-15_dp,'paired changed gas face conserves inventory')
    call MPI_Allreduce(held_change,held_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.abs(held_total)<1.e-15_dp,'paired changed held face conserves inventory')
    call require(held_history%initial_energy_correction==1.25_dp.and. &
      held_history%numerical_energy_correction==2.5_dp,'face reconstruction preserves independent reset ledger')
    kept_reconstructed_gas=reconstructed_gas;kept_reconstructed_held=reconstructed_held
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,held_history,changed_gas_faces,nv, &
      reconstructed_gas,MPI_COMM_WORLD,good)
    call require(.not.good,'held interval rejects gas-only reconstruction')
    if(rank==nranks-1) held_history%numerical_energy_correction=2.75_dp
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,held_history,changed_gas_faces,nv, &
      reconstructed_gas,MPI_COMM_WORLD,good,internal_energy_flux=changed_held_faces, &
      internal_energy_rate=reconstructed_held)
    call require(.not.good,'owner/idle interval correction mismatch rejects collectively')
    held_history%numerical_energy_correction=2.5_dp
    if(rank==1) held_history%initial_internal_energy%fields(3)%values(1,1,2,2)= &
      ieee_value(0._dp,ieee_quiet_nan)
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,held_history,changed_gas_faces,nv, &
      reconstructed_gas,MPI_COMM_WORLD,good,internal_energy_flux=changed_held_faces, &
      internal_energy_rate=reconstructed_held)
    call require(.not.good,'remote invalid held interval rejects both reconstructed families')
    do p=1,3
      if(ownership(p)/=rank) cycle
      call require(all(reconstructed_gas%fields(p)%values==kept_reconstructed_gas%fields(p)%values).and. &
        all(reconstructed_held%fields(p)%values==kept_reconstructed_held%fields(p)%values), &
        'paired reconstruction failure retains all owner outputs')
    enddo
    ! An absent species freezes this interior face. The original intensive
    ! rate and the integrated correction have different rounding histories;
    ! their subtraction used to leave a negative species at zero inventory.
    f=0;f(6)=23._dp/997._dp;f(7)=-f(6)
    serial_faces(1)=serial_low(1);serial_rate%fields(1)%values=0
    serial_faces(1)%x(:,2,1,1)=dt*f
    serial_rate%fields(1)%values(:,2,1,1)=-f/2
    serial_rate%fields(1)%values(:,3,1,1)=f/2
    serial_rate%fields(1)%values(2,:,:,:)=3
    if(rank==0) then
      high_faces(1)=low_faces(1);patch_rate%fields(1)%values=serial_rate%fields(1)%values
      high_faces(1)%x(:,2,1,1)=dt*f
    endif
    call apply_eb_tree_flux_correction_3d(serial_tree,serial_source,serial_low,dt,serial_rate,serial_faces,good, &
      wrap,target_depth=0)
    call require(good,'zero-inventory serial correction avoids cancellation residual')
    call apply_mpi_eb_flux_correction_3d(layout,ownership,metrics,patch_state,patch_source,low_faces,dt, &
      patch_rate,high_faces,MPI_COMM_WORLD,good,wrap,target_depth=0)
    call require(good,'zero-inventory MPI correction avoids cancellation residual')
    if(rank==0) then
      call require(all(high_faces(1)%x==0).and.all(patch_rate%fields(1)%values(6:,:,:,:)==0), &
        'frozen species face has exactly zero flux and divergence')
      call require(all(patch_rate%fields(1)%values(2,:,:,:)==3),'correction preserves independent wall momentum')
      call require(all(patch_rate%fields(1)%values==serial_rate%fields(1)%values).and. &
        all(high_faces(1)%x==serial_faces(1)%x),'zero-inventory serial/MPI corrections match exactly')
    endif
    call check_held_interval_capture(layout,ownership,metrics)
  end subroutine check_flux_allowance_exchange

  ! Pure N2 through the valid fixed mechanism has a linear inert trajectory.
  ! This checks real interval capture separately from the detailed
  ! chemistry qualification and the manufactured face-change oracle above.
  subroutine check_held_interval_capture(layout,ownership,metrics)
    type(amr_patch_tree_3d), intent(in) :: layout
    integer, intent(in) :: ownership(:)
    type(eb_geometry_3d), intent(in) :: metrics(:)
    type(eb_geometry_3d), allocatable :: padded(:)
    type(elementary_reaction), allocatable :: reaction_set(:)
    type(gas_transport_species), allocatable :: coefficients(:)
    type(amr_tree_state_3d) :: gas,temp,held,restored,restored_held
    type(mpi_eb_mol_interval_3d) :: archive,kept_archive
    type(mpi_eb_mol_interval_3d), allocatable :: clocks(:),kept_clocks(:)
    type(amr_tree_state_3d) :: seed_gas,seed_temp,seed_held,kept_gas,kept_temp,kept_held
    real(dp) :: before(nv),after(nv),accepted_correction
    integer(int64) :: level_steps(2)
    integer :: iterations,r,total_refusals
    type(eb_integrated_flux_3d), allocatable :: kept_faces(:)
    type(eb_integrated_flux_3d), allocatable :: gas_faces(:),held_faces(:)
    type(dual_energy_policy) :: policy,rank_policy
    real(dp) :: pv(nv),uniform(nv),uniform_t,sound,stored,b(nv),w(nv),ch(nv),ext(nv),theta,correction
    integer :: p,n(3),sp
    logical :: good
    character(len=16) :: faces(6)
    faces='outflow';pv=0;pv(1)=.31_dp
    do sp=1,size(species)
      select case(trim(species(sp)%name))
      case('N2');pv(5+sp)=1
      end select
    enddo
    pv(2:4)=[2._dp,-1._dp,.5_dp]
    pv(5)=pv(1)*1000*sum(pv(6:)*8.31446261815324e3_dp/species%molecular_weight)
    call reactive_primitive_to_conserved(species,pv,uniform,uniform_t,sound,good)
    call require(good,'held interval native linear fixture EOS')
    stored=uniform(5)-.5_dp*sum(uniform(2:4)**2)/uniform(1)
    call load_h2o2_elementary_transport(coefficients,good)
    call require(good,'held interval native linear fixture transport')
    call load_h2o2_elementary_mechanism(reaction_set,good)
    call require(good,'held interval valid fixed mechanism')
    allocate(gas%fields(size(ownership)),temp%fields(size(ownership)), &
      held%fields(size(ownership)),padded(size(ownership)))
    do p=1,size(ownership)
      n=layout%nodes(p)%shape
      associate(g=>metrics(p))
      call build_axis_plane_eb_geometry_3d(n(1)+12,n(2)+12,n(3)+12, &
        g%x_lower-6*g%dx,g%x_lower+(n(1)+6)*g%dx, &
        g%y_lower-6*g%dy,g%y_lower+(n(2)+6)*g%dy, &
        g%z_lower-6*g%dz,g%z_lower+(n(3)+6)*g%dz,'x',-1000._dp,padded(p),good)
      end associate
      call require(good,'held interval regular grow6 geometry')
      if(ownership(p)/=rank) cycle
      allocate(gas%fields(p)%values(nv,n(1),n(2),n(3)),temp%fields(p)%values(1,n(1),n(2),n(3)), &
        held%fields(p)%values(1,n(1),n(2),n(3)))
      do sp=1,nv
        gas%fields(p)%values(sp,:,:,:)=uniform(sp)
      enddo
      gas%fields(p)%values(5,:,:,:)=uniform(5)-11
      temp%fields(p)%values=uniform_t;held%fields(p)%values=stored
    enddo
    policy%allow_negative_energy=.false.;policy%eta2=2;correction=7
    seed_gas=gas;seed_temp=temp;seed_held=held
    call advance_mpi_eb_tree_sdc_3d(species,reaction_set,coefficients,layout,ownership,metrics,gas,temp,padded, &
      1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp, &
      faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,target_depth=1,integrated_flux=gas_faces, &
      mol_mode=.true.,internal_energy=held,internal_energy_flux=held_faces,energy_policy=policy, &
      energy_correction=correction,interval_history=archive)
    call require(good,'native linear held common interval capture')
    call require(archive%has_internal_energy.and.archive%depth==1.and.archive%internal_energy_policy%eta2==2, &
      'accepted held archive retains source level and reset policy')
    call require(archive%initial_energy_correction>0.and.correction==7+archive%numerical_energy_correction, &
      'accepted held archive retains actual initial repair and numerical ledger')
    do p=1,size(ownership)
      call require(allocated(archive%initial_internal_energy%fields(p)%values).eqv. &
        (ownership(p)==rank.and.layout%nodes(p)%level==1),'held interval owner/level sparse archive')
      if(ownership(p)/=rank.or.layout%nodes(p)%level/=1) cycle
      call require(all(archive%initial_internal_energy%fields(p)%values==stored).and. &
        all(archive%local_internal_energy%fields(p)%values==held%fields(p)%values), &
        'accepted held archive retains reactor initial energy and published endpoint')
      call require(all(archive%internal_energy_flux(p)%x==held_faces(p)%x).and. &
        all(archive%internal_energy_flux(p)%y==held_faces(p)%y).and. &
        all(archive%internal_energy_flux(p)%z==held_faces(p)%z),'accepted held archive retains actual spatial faces')
    enddo
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,archive,gas_faces,nv,restored, &
      MPI_COMM_WORLD,good,internal_energy_flux=held_faces,internal_energy_rate=restored_held)
    call require(good,'accepted native held interval can reconstruct both actual source families')
    do p=1,size(ownership)
      if(ownership(p)/=rank.or.layout%nodes(p)%level/=1) cycle
      call require(all(restored%fields(p)%values==archive%forcing%fields(p)%values).and. &
        all(restored_held%fields(p)%values==archive%internal_energy_spatial_rate%fields(p)%values), &
        'native interval identity reconstruction preserves admitted gas and held rates exactly')
    enddo
    call check_native_held_interval_replay(layout,ownership,metrics,reaction_set,archive,seed_gas,seed_temp, &
      seed_held,gas_faces,held_faces)
    kept_archive=archive;rank_policy=policy
    if(rank==nranks-1) rank_policy%eta2=3
    call advance_mpi_eb_tree_sdc_3d(species,reaction_set,coefficients,layout,ownership,metrics,gas,temp,padded, &
      1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp, &
      faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,target_depth=1,integrated_flux=gas_faces, &
      mol_mode=.true.,internal_energy=held,internal_energy_flux=held_faces,energy_policy=rank_policy, &
      energy_correction=correction,interval_history=archive)
    call require(.not.good.and.correction==7+kept_archive%numerical_energy_correction, &
      'held policy rejection retains prior interval correction')
    do p=1,size(ownership)
      if(ownership(p)/=rank.or.layout%nodes(p)%level/=1) cycle
      call require(all(archive%initial_internal_energy%fields(p)%values== &
        kept_archive%initial_internal_energy%fields(p)%values).and. &
        all(archive%local_internal_energy%fields(p)%values==kept_archive%local_internal_energy%fields(p)%values), &
        'held policy rejection preserves accepted archive payload')
    enddo
    ! Reusing the destination for a gas-only interval must clear every old
    ! held payload rather than retaining a stale policy or correction record.
    call advance_mpi_eb_tree_sdc_3d(species,reaction_set,coefficients,layout,ownership,metrics,gas,temp,padded, &
      1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],2.e-7_dp,1.e-12_dp, &
      faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good,target_depth=1,integrated_flux=gas_faces, &
      mol_mode=.true.,interval_history=archive)
    call require(good.and..not.archive%has_internal_energy,'gas-only interval replaces previous held archive')
    call require(.not.allocated(archive%initial_internal_energy%fields).and. &
      .not.allocated(archive%local_internal_energy%fields).and. &
      .not.allocated(archive%old_internal_energy_spatial_rate%fields).and. &
      .not.allocated(archive%internal_energy_spatial_rate%fields).and. &
      .not.allocated(archive%old_internal_energy_flux).and..not.allocated(archive%internal_energy_flux).and. &
      .not.allocated(archive%clock_internal_energy%fields).and.archive%clock_energy_correction==0.and. &
      archive%initial_energy_correction==0.and.archive%numerical_energy_correction==0.and. &
      archive%forcing_energy_correction==0, &
      'gas-only replacement clears all held state/face/rate and correction payload')
    call reconstruct_mpi_eb_mol_interval_3d(layout,ownership,metrics,archive,gas_faces,nv,restored,MPI_COMM_WORLD,good)
    call require(good,'gas-only interval reconstruction remains compatible after held archive replacement')
    ! A real recursive native chemistry run must archive every owner-local
    ! level at clock start, while the common interval payload stays sparse.
    ! Pure N2 makes this a short analytical chemistry fixture; no mock solver.
    do iterations=1,2
      gas=seed_gas;temp=seed_temp;held=seed_held;correction=7
      call integrate_mpi_eb_tree_3d(layout,ownership,metrics,gas,before,MPI_COMM_WORLD,good)
      call require(good,'held recursive inert initial inventory')
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reaction_set,coefficients,layout,ownership,metrics,gas,temp, &
        1.e-10_dp,iterations,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        2.e-7_dp,1.e-12_dp,padded,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good, &
        mol_mode=.true.,internal_energy=held,internal_energy_flux=held_faces,energy_policy=policy, &
        energy_correction=correction,interval_history=clocks,level_steps=level_steps)
      call require(good,'native inert held recursive interval capture')
      call require(all(level_steps==[1_int64,2_int64]).and.size(clocks)==3, &
        'held recursive inert native clock counts')
      call require(.not.any(clocks%replayed),'inert held recursive clock does not manufacture a replay')
      call integrate_mpi_eb_tree_3d(layout,ownership,metrics,gas,after,MPI_COMM_WORLD,good)
      call require(good,'held recursive inert final inventory')
      after=after-before-b-w-ch-ext;after(5)=after(5)-(correction-7)
      call require(maxval(abs(after)/max(1._dp,abs(before)))<3.e-14_dp, &
        'held recursive inert inventory separates numeric repairs from chemistry')
      call require(clocks(1)%clock_energy_correction>7.and. &
        ieee_is_finite(clocks(1)%clock_energy_correction), &
        'full clock correction prefix includes initial composite repair')
      do r=1,size(clocks)
        call require(clocks(r)%has_internal_energy.and. &
          allocated(clocks(r)%clock_initial%fields).and. &
          allocated(clocks(r)%clock_temperature%fields).and. &
          allocated(clocks(r)%clock_reaction_source%fields).and. &
          allocated(clocks(r)%clock_internal_energy%fields), &
          'held recursive archive includes complete gas and scalar clock')
        do p=1,size(ownership)
          call require(allocated(clocks(r)%clock_internal_energy%fields(p)%values).eqv. &
            (ownership(p)==rank),'held recursive full clock includes inactive levels and excludes remote owners')
          if(ownership(p)/=rank) cycle
          call require(all(ieee_is_finite(clocks(r)%clock_internal_energy%fields(p)%values)), &
            'held recursive full clock scalar is finite')
          if(layout%nodes(p)%level/=clocks(r)%depth) cycle
          call require(all(clocks(r)%clock_internal_energy%fields(p)%values== &
            clocks(r)%initial_internal_energy%fields(p)%values), &
            'held recursive full clock anchors actual reactor initial internal energy')
        enddo
      enddo
      ! All-rank policy refusal must retain the accepted outer clock family,
      ! correction prefix, gas, temperature and scalar as one transaction.
      kept_clocks=clocks;kept_gas=gas;kept_temp=temp;kept_held=held;accepted_correction=correction
      rank_policy=policy
      if(rank==nranks-1) rank_policy%eta2=3
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reaction_set,coefficients,layout,ownership,metrics,gas,temp, &
        1.e-10_dp,iterations,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        2.e-7_dp,1.e-12_dp,padded,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good, &
        mol_mode=.true.,internal_energy=held,internal_energy_flux=held_faces,energy_policy=rank_policy, &
        energy_correction=correction,interval_history=clocks,level_steps=level_steps)
      call require(.not.good.and.correction==accepted_correction, &
        'held recursive idle/owner policy refusal keeps accepted correction')
      do p=1,size(ownership)
        if(ownership(p)/=rank) cycle
        call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
          all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
          all(held%fields(p)%values==kept_held%fields(p)%values), &
          'held recursive policy refusal keeps all three accepted state fields')
        do r=1,size(clocks)
          call require(clocks(r)%clock_energy_correction==kept_clocks(r)%clock_energy_correction.and. &
            all(clocks(r)%clock_internal_energy%fields(p)%values== &
              kept_clocks(r)%clock_internal_energy%fields(p)%values), &
            'held recursive policy refusal retains accepted full clock prefix and scalar')
        enddo
      enddo
    enddo
    ! Refuse a native-admitted child reaction after the parent has completed.
    ! Initial and parent repairs must stay private, and the provider refusal
    ! must not be retried as a failed physical closure.
    gas=seed_gas;temp=seed_temp;held=seed_held;correction=7
    kept_gas=gas;kept_temp=temp;kept_held=held;kept_faces=held_faces
    held_clock_reactor_calls=0;held_clock_refusals=0
    call advance_mpi_eb_tree_sdc_subcycle_3d(species,reaction_set,coefficients,layout,ownership,metrics,gas,temp, &
      1.e-10_dp,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
      2.e-7_dp,1.e-12_dp,padded,faces,b,w,ch,ext,theta,MPI_COMM_WORLD,good, &
      mol_mode=.true.,internal_energy=held,internal_energy_flux=held_faces,energy_policy=policy, &
      energy_correction=correction,interval_history=clocks,level_steps=level_steps, &
      evaluate_reaction=held_clock_refusing_reactor)
    call require(.not.good.and.correction==7,'late held recursive provider refusal discards all private numeric repairs')
    call MPI_Allreduce(held_clock_refusals,total_refusals,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(total_refusals==1,'held recursive native-admitted provider refusal occurs exactly once')
    if(rank==1) call require(held_clock_reactor_calls==2,'held recursive child provider refusal is not retried')
    do p=1,size(ownership)
      if(ownership(p)/=rank) cycle
      call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
        all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
        all(held%fields(p)%values==kept_held%fields(p)%values), &
        'late held recursive provider refusal keeps the original gas/T/scalar')
      call require(all(held_faces(p)%x==kept_faces(p)%x).and.all(held_faces(p)%y==kept_faces(p)%y).and. &
        all(held_faces(p)%z==kept_faces(p)%z),'late held recursive provider refusal keeps previous accepted scalar faces')
      do r=1,size(clocks)
        call require(clocks(r)%clock_energy_correction==kept_clocks(r)%clock_energy_correction.and. &
          all(clocks(r)%clock_internal_energy%fields(p)%values== &
            kept_clocks(r)%clock_internal_energy%fields(p)%values), &
          'late held recursive provider refusal keeps previous accepted clock history')
      enddo
    enddo
    if(rank==0) print *, 'PASS: native inert held common/recursive archives, paired reconstruction and rollback'
  end subroutine check_held_interval_capture

  ! The same atomic native replay kernel is called by failed-clock recursion.
  ! Change a matched interior face pair with distinct signed gas/held impulses:
  ! scalar spatial reconstruction must not become reactor energy forcing.
  subroutine check_native_held_interval_replay(layout,ownership,metrics,reactions,accepted,seed_gas,seed_temp, &
      seed_held,original_gas_faces,original_held_faces)
    type(amr_patch_tree_3d), intent(in) :: layout
    integer, intent(in) :: ownership(:)
    type(eb_geometry_3d), intent(in) :: metrics(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(mpi_eb_mol_interval_3d), intent(in) :: accepted
    type(amr_tree_state_3d), intent(in) :: seed_gas,seed_temp,seed_held
    type(eb_integrated_flux_3d), intent(in) :: original_gas_faces(:),original_held_faces(:)
    type(mpi_eb_mol_interval_3d) :: record,kept_record
    type(amr_tree_state_3d) :: gas,temp,held,ir,old,old_held,kept_gas,kept_temp,kept_held,kept_ir
    type(amr_tree_state_3d) :: heat,kept_heat
    type(eb_integrated_flux_3d), allocatable :: gas_faces(:),held_faces(:)
    real(dp) :: chemical(nv),correction,kept_correction,v0(nv),v1(nv),expected(nv),expected_ie,scalar_prediction
    real(dp) :: gas_impulse,held_impulse,volume,area,projection_change,expected_projection
    real(dp) :: baseline_ie,baseline_scalar,rounded_scalar
    integer :: p,n(3),face,cell_index,total_refusals
    logical :: good
    call require(allocated(accepted%initial_reaction_source%fields), &
      'ordinary accepted interval archives its actual initial reaction source')
    gas=seed_gas;temp=seed_temp;held=seed_held
    allocate(ir%fields(size(ownership)),heat%fields(size(ownership)))
    do p=1,size(ownership)
      if(ownership(p)/=rank) cycle
      n=layout%nodes(p)%shape
      allocate(ir%fields(p)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
      allocate(heat%fields(p)%values(1,n(1),n(2),n(3)),source=27._dp)
      if(layout%nodes(p)%level/=accepted%depth) cycle
      gas%fields(p)=accepted%initial%fields(p);temp%fields(p)=accepted%initial_temperature%fields(p)
      held%fields(p)=accepted%initial_internal_energy%fields(p);ir%fields(p)=accepted%initial_reaction_source%fields(p)
    enddo
    old=gas;old_held=held;record=accepted
    gas_faces=original_gas_faces;held_faces=original_held_faces
    do p=2,3
      if(ownership(p)/=rank) cycle
      face=merge(4,0,p==2)
      gas_faces(p)%x(5,face,2,2)=gas_faces(p)%x(5,face,2,2)+2.5e6_dp*accepted%interval
      ! An independently changed scalar face on one side leaves a nonzero
      ! integrated gas-to-scalar projection, exposing a stale archive value.
      if(p==2) held_faces(p)%x(1,face,2,2)=held_faces(p)%x(1,face,2,2)-5.e6_dp*accepted%interval
    enddo
    call integrate_mpi_eb_tree_3d(layout,ownership,metrics,gas,v0,MPI_COMM_WORLD,good)
    call require(good,'native changed-pair replay initial independent inventory')
    correction=17;chemical=123
    call replay_mpi_eb_mol_held_interval_3d(species,reactions,layout,ownership,metrics,record,gas_faces,held_faces, &
      gas,temp,held,ir,correction,chemical,2.e-7_dp,1.e-12_dp,MPI_COMM_WORLD,good,heat_release=heat)
    call require(good.and.record%replayed,'native changed-pair interval replay accepts its actual replacement solve')
    do p=1,size(ownership)
      call require(allocated(record%internal_energy_flux(p)%x).eqv. &
        (ownership(p)==rank.and.layout%nodes(p)%level==accepted%depth), &
        'native replay publishes only owner/level sparse accepted scalar faces')
      if(ownership(p)/=rank) cycle
      call require(all(record%clock_internal_energy%fields(p)%values==old_held%fields(p)%values), &
        'native replay full clock retains old scalar on every owner level')
      call require(all(record%clock_heat%fields(p)%values==27), &
        'native replay full clock retains the actual complete heat history before replacement')
      if(layout%nodes(p)%level==accepted%depth) then
        call require(all(heat%fields(p)%values==0),'native inert replacement publishes actual zero chemical heat')
      else
        call require(all(heat%fields(p)%values==27),'native replay keeps heat history on inactive owner levels')
      endif
    enddo
    projection_change=0
    do p=2,3
      if(ownership(p)/=rank) cycle
      face=merge(4,0,p==2);cell_index=merge(4,1,p==2)
      volume=metrics(p)%dx*metrics(p)%dy*metrics(p)%dz*metrics(p)%volume_fraction(cell_index,2,2)
      area=metrics(p)%dy*metrics(p)%dz*metrics(p)%x_face_fraction(face,2,2)
      gas_impulse=merge(-1._dp,1._dp,p==2)*area/volume* &
        (gas_faces(p)%x(5,face,2,2)-original_gas_faces(p)%x(5,face,2,2))
      held_impulse=merge(-1._dp,1._dp,p==2)*area/volume* &
        (held_faces(p)%x(1,face,2,2)-original_held_faces(p)%x(1,face,2,2))
      expected=old%fields(p)%values(:,cell_index,2,2)+ &
        accepted%interval*accepted%forcing%fields(p)%values(:,cell_index,2,2)
      baseline_ie=expected(5)-.5_dp*sum(expected(2:4)**2)/expected(1)
      baseline_scalar=old_held%fields(p)%values(1,cell_index,2,2)+ &
        accepted%interval*accepted%internal_energy_spatial_rate%fields(p)%values(1,cell_index,2,2)
      expected(5)=expected(5)+gas_impulse
      expected_ie=expected(5)-.5_dp*sum(expected(2:4)**2)/expected(1)
      ! Compare the rounded cell updates, rather than their infinitesimal
      ! face difference: the finite impulse is small beside the held energy.
      rounded_scalar=baseline_scalar+held_impulse
      projection_change=projection_change+volume*((expected_ie-baseline_ie)-(rounded_scalar-baseline_scalar))
      call require(maxval(abs(gas%fields(p)%values(:,cell_index,2,2)-expected)/max(1._dp,abs(expected)))<3.e-14_dp, &
        'native paired replay gas endpoint matches independent signed gas-face impulse')
      call require(abs(held%fields(p)%values(1,cell_index,2,2)-expected_ie)<3.e-14_dp*max(1._dp,abs(expected_ie)), &
        'native paired replay UEINT endpoint uses projected gas forcing and held old initial energy')
      scalar_prediction=old_held%fields(p)%values(1,cell_index,2,2)+held_impulse
      call require(abs(held%fields(p)%values(1,cell_index,2,2)-scalar_prediction)>.5_dp*abs(gas_impulse), &
        'native replay does not substitute independently signed held spatial RHS for reactor energy forcing')
      call require(abs(accepted%interval*(record%internal_energy_spatial_rate%fields(p)%values(1,cell_index,2,2)- &
        accepted%internal_energy_spatial_rate%fields(p)%values(1,cell_index,2,2))-held_impulse)< &
        3.e-14_dp*max(1._dp,abs(held_impulse)), &
        'native replay archive retains the distinct accepted scalar spatial impulse')
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,projection_change,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS,'independent scalar projection integral')
    expected_projection=accepted%internal_energy_projection_change+projection_change
    call require(abs(record%internal_energy_projection_change-expected_projection)< &
      3.e-14_dp*max(1._dp,abs(expected_projection)), &
      'MOL replay recomputes projection from independently changed gas and scalar face impulses')
    call require(abs(record%internal_energy_projection_change-accepted%internal_energy_projection_change)> &
      .5_dp*abs(projection_change),'accepted MOL replay replaces its previous scalar projection')
    call integrate_mpi_eb_tree_3d(layout,ownership,metrics,gas,v1,MPI_COMM_WORLD,good)
    call require(good,'native changed-pair replay final independent inventory')
    v1=v1-v0-chemical;v1(5)=v1(5)-(correction-17)
    call require(maxval(abs(v1)/max(1._dp,abs(v0)))<3.e-14_dp, &
      'native paired replay separates chemical/numeric ledgers and conserves internal interface inventory')
    kept_gas=gas;kept_temp=temp;kept_held=held;kept_ir=ir;kept_record=record;kept_correction=correction
    kept_heat=heat
    ! Projection metadata is owner-neutral, including the idle rank. An
    ! inconsistent archive must refuse before invoking any reactor provider.
    if(rank==nranks-1) record%internal_energy_projection_change=record%internal_energy_projection_change+1
    expected_projection=record%internal_energy_projection_change
    chemical=321;held_clock_reactor_calls=0
    call replay_mpi_eb_mol_held_interval_3d(species,reactions,layout,ownership,metrics,record,gas_faces,held_faces, &
      gas,temp,held,ir,correction,chemical,2.e-7_dp,1.e-12_dp,MPI_COMM_WORLD,good, &
      evaluate_reaction=held_clock_refusing_reactor,heat_release=heat)
    call require(.not.good.and.held_clock_reactor_calls==0.and.all(chemical==321).and. &
      record%internal_energy_projection_change==expected_projection.and.correction==kept_correction, &
      'owner/idle scalar projection mismatch refuses before chemistry and preserves the supplied archive')
    record=kept_record

    ! A mismatched prefix on an idle rank must refuse the entire replay before
    ! any native provider call, keeping that rank's own input prefix unchanged.
    if(rank==nranks-1) correction=correction+1
    kept_correction=correction;chemical=321;held_clock_reactor_calls=0
    call replay_mpi_eb_mol_held_interval_3d(species,reactions,layout,ownership,metrics,record,gas_faces,held_faces, &
      gas,temp,held,ir,correction,chemical,2.e-7_dp,1.e-12_dp,MPI_COMM_WORLD,good, &
      evaluate_reaction=held_clock_refusing_reactor,heat_release=heat)
    call require(.not.good.and.correction==kept_correction.and.all(chemical==321).and.held_clock_reactor_calls==0, &
      'native replay owner/idle numeric-prefix mismatch refuses before chemistry and keeps all caller ledgers')
    ! Refuse after native admission, after private selected-clock repairs have
    ! already happened. The new kernel must publish none of that private work.
    correction=17
    do p=1,size(ownership)
      if(ownership(p)/=rank.or.layout%nodes(p)%level/=accepted%depth) cycle
      gas%fields(p)%values(5,:,:,:)=gas%fields(p)%values(5,:,:,:)-33
    enddo
    kept_gas=gas;held_clock_reactor_calls=0;held_clock_refusals=0
    call replay_mpi_eb_mol_held_interval_3d(species,reactions,layout,ownership,metrics,record,gas_faces,held_faces, &
      gas,temp,held,ir,correction,chemical,2.e-7_dp,1.e-12_dp,MPI_COMM_WORLD,good, &
      evaluate_reaction=held_clock_refusing_reactor,heat_release=heat)
    call require(.not.good.and.correction==17.and.all(chemical==321), &
      'native replay provider refusal discards private clock resets and chemical/numeric changes')
    call MPI_Allreduce(held_clock_refusals,total_refusals,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(total_refusals==1,'native replay provider refused one input proved admissible by native chemistry')
    if(rank==1) call require(held_clock_reactor_calls==2,'native interval provider refusal is not retried')
    call require(record%internal_energy_projection_change==kept_record%internal_energy_projection_change, &
      'refused MOL replay retains its accepted scalar projection')
    do p=1,size(ownership)
      if(ownership(p)/=rank) cycle
      call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
        all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
        all(held%fields(p)%values==kept_held%fields(p)%values).and.all(ir%fields(p)%values==kept_ir%fields(p)%values), &
        'native replay provider refusal preserves gas/T/held/reaction history atomically')
      call require(all(heat%fields(p)%values==kept_heat%fields(p)%values), &
        'native replay provider refusal retains complete chemical heat history')
      call require(record%clock_energy_correction==kept_record%clock_energy_correction.and. &
        all(record%clock_internal_energy%fields(p)%values==kept_record%clock_internal_energy%fields(p)%values), &
        'native replay provider refusal preserves accepted full clock and numeric prefix')
    enddo
    ! Optional heat output can be omitted on a later accepted replacement.
    ! Its clock must describe that actual invocation, clearing the previous
    ! full heat payload rather than silently retaining an unrelated old clock.
    call replay_mpi_eb_mol_held_interval_3d(species,reactions,layout,ownership,metrics,record,gas_faces,held_faces, &
      gas,temp,held,ir,correction,chemical,2.e-7_dp,1.e-12_dp,MPI_COMM_WORLD,good)
    call require(good.and..not.allocated(record%clock_heat%fields), &
      'accepted native replay without heat output replaces its previous heat clock payload')
    if(rank==0) print *, 'PASS: atomic native held interval replay, independent paired impulses and provider rollback'
  end subroutine check_native_held_interval_replay

  subroutine held_clock_refusing_reactor(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,good,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    held_clock_reactor_calls=held_clock_reactor_calls+1
    call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,good,maximum_steps=maximum_steps)
    if(good.and.rank==1.and.held_clock_reactor_calls==2) then
      held_clock_refusals=held_clock_refusals+1
      species_density=-1;temperature=-1;completed_steps=0;good=.false.
    endif
  end subroutine held_clock_refusing_reactor


  subroutine check_mol_only()
    type(amr_tree_state_3d) :: local,lt
    integer :: owners(count),p
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(local%fields(count),lt%fields(count))
    do p=1,count
      if(owners(p)/=rank) cycle
      local%fields(p)%values=saved%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
    end do
    call check_mol_sources(saved,owners,local,lt)
  end subroutine check_mol_only
  subroutine check_rotating_mol()
    type(eb_geometry_3d) :: three(count)
    type(amr_tree_state_3d) :: local,lt,rates
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    type(eb_integrated_flux_3d) :: reference_flux
    type(rotating_frame_3d) :: mismatched
    real(dp) :: ghost(nv,10,10,10),gt(10,10,10),vf(10,10,10),expected(nv,4,4,4),b(nv),w(nv)
    real(dp) :: prescribed(nv,6),tin(6),mach(6),ordinary_energy
    character(len=16) :: faces(6)
    integer :: owners(count),p,mode
    logical :: good
    owners=[(mod(p-1,nranks),p=1,count)]
    allocate(local%fields(count),lt%fields(count))
    do p=1,count
      call extract_eb_geometry_3d(extended(p),[4,4,4],[13,13,13],three(p),good)
      call require(good,'rotating three-layer geometry')
      if(owners(p)/=rank) cycle
      local%fields(p)%values=saved%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
      call fill_reactive_eb_tree_halo_3d(species,saved,p,ghost,gt,good, &
        halo_width=3,target_volume_fraction=three(p)%volume_fraction,ghost_volume_fraction=vf, &
        restrict_snapshot=.false.,frame=frame)
      call require(good,'rotating serial MOL donors before MPI')
      call compute_mol_eb_hydro_3d(species,three(p),ghost,gt,1e-8_dp,2,expected,reference_flux,b,w,good,frame=frame)
      call require(good,'rotating serial MOL rate before MPI')
    end do
    if(nranks>1) then
      mismatched=frame
      if(rank==1) mismatched%angular_speed=mismatched%angular_speed+1._dp
      call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,2, &
        rates,flux,MPI_COMM_WORLD,good,frame=mismatched)
      call require(.not.good,'rank-mismatched MOL frame rejected')
      call require(.not.allocated(rates%fields).and..not.allocated(flux),'MOL mismatch rollback')
    end if
    prescribed(:,1)=saved%nodes(1)%state(:,1,1,3)
    position=[physical(1)%x_lower,physical(1)%y_lower,physical(1)%z_lower]+ &
      (real([1,1,3],dp)-.5_dp)*[physical(1)%dx,physical(1)%dy,physical(1)%dz]
    ordinary_energy=prescribed(5,1)
    call convert_rotating_energy_3d(frame,position,prescribed(1,1),prescribed(5,1),.false.,ordinary_energy,good)
    call require(good,'ordinary inlet energy')
    prescribed(5,1)=ordinary_energy
    do p=2,6
      prescribed(:,p)=prescribed(:,1)
    end do
    tin=saved%nodes(1)%temperature(1,1,3)
    faces='outflow';faces(1)='inflow';mach=-1._dp;mach(1)=.1_dp
    do mode=1,3
      select case(mode)
      case(1)
        call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,2, &
          rates,flux,MPI_COMM_WORLD,good,frame=frame)
      case(2)
        call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,2, &
          rates,flux,MPI_COMM_WORLD,good,boundary_faces=faces,inflow_state=prescribed,inflow_temperature=tin,frame=frame)
      case(3)
        call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,2, &
          rates,flux,MPI_COMM_WORLD,good,inflow_mach=mach,inflow_temperature=tin,frame=frame)
      end select
      call require(good,'rotating distributed MOL rate')
    do p=1,count
      if(owners(p)/=rank) cycle
      select case(mode)
      case(1)
        call fill_reactive_eb_tree_halo_3d(species,saved,p,ghost,gt,good, &
          halo_width=3,target_volume_fraction=three(p)%volume_fraction,ghost_volume_fraction=vf, &
          restrict_snapshot=.false.,frame=frame)
      case(2)
        call fill_reactive_eb_tree_halo_3d(species,saved,p,ghost,gt,good, &
          halo_width=3,target_volume_fraction=three(p)%volume_fraction,ghost_volume_fraction=vf, &
          restrict_snapshot=.false.,boundary_faces=faces,inflow_state=prescribed,inflow_temperature=tin,frame=frame)
      case(3)
        call fill_reactive_eb_tree_halo_3d(species,saved,p,ghost,gt,good, &
          halo_width=3,target_volume_fraction=three(p)%volume_fraction,ghost_volume_fraction=vf, &
          restrict_snapshot=.false.,inflow_mach=mach,inflow_temperature=tin,frame=frame)
      end select
      call require(good,'rotating serial MOL donors')
      call compute_mol_eb_hydro_3d(species,three(p),ghost,gt,1e-8_dp,2,expected,reference_flux,b,w,good,frame=frame)
      call require(good,'rotating serial MOL rate')
      call require(maxval(abs(rates%fields(p)%values-expected)/max(1._dp,abs(expected)))<1e-12_dp, &
        'rotating serial/MPI MOL rate')
      call require(maxval(abs(flux(p)%x-reference_flux%x)/max(1._dp,abs(reference_flux%x)))<1e-12_dp.and. &
        maxval(abs(flux(p)%y-reference_flux%y)/max(1._dp,abs(reference_flux%y)))<1e-12_dp.and. &
        maxval(abs(flux(p)%z-reference_flux%z)/max(1._dp,abs(reference_flux%z)))<1e-12_dp, &
        'rotating serial/MPI MOL flux')
    end do
    end do
  end subroutine check_rotating_mol
  subroutine check_mol_sources(stage_tree,owners,local,lt)
    type(eb_reactive_tree_3d), intent(in) :: stage_tree
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: local,lt
    type(eb_reactive_tree_3d) :: lf_tree
    type(eb_geometry_3d) :: three(count)
    type(amr_tree_state_3d) :: rates,prior,lf_local,lf_temperature,serial_rates,candidate,candidate_temperature,saved_temperature
    type(amr_tree_state_3d) :: signed_endpoint,signed_increment
    type(eb_integrated_flux_3d), allocatable :: flux(:),previous(:),serial_flux(:)
    type(eb_integrated_flux_3d), allocatable :: audited_flux(:)
    type(eb_integrated_flux_3d) :: reference_flux
    real(dp) :: ghost(nv,10,10,10),gt(10,10,10),vf(10,10,10),expected(nv,4,4,4),b(nv),w(nv)
    real(dp) :: qp(reactive_nprim(size(species))),temp_cell,speed,cell_budget,face_budget
    real(dp) :: left_face(nv,4,4),right_face(nv,4,4)
    character(len=16) :: seam_faces(6)
    integer :: p,wrong_order,h2,n2,x,y,z,changed,global_changed,depth,repairs,negative_endpoints,global_negative
    real(dp) :: signed_inventory(nv),accepted_inventory(nv),face_change(nv),observed_change(nv)
    real(dp) :: face_residual,global_face_residual,inventory_residual
    logical :: good
    do p=1,count
      call extract_eb_geometry_3d(extended(p),[4,4,4],[13,13,13],three(p),good)
      call require(good,'MOL three-layer geometry')
    end do
    call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,2, &
      rates,flux,MPI_COMM_WORLD,good,boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(good,'distributed MOL rate')
    do p=1,count
      call require(allocated(rates%fields(p)%values).eqv.(owners(p)==rank),'owner-local MOL rate')
      call require(allocated(flux(p)%x).eqv.(owners(p)==rank),'owner-local MOL flux')
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,stage_tree,p,ghost,gt,good, &
        halo_width=3,target_volume_fraction=three(p)%volume_fraction,ghost_volume_fraction=vf, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t,restrict_snapshot=.false.)
      call require(good,'serial MOL donors')
      call compute_mol_eb_hydro_3d(species,three(p),ghost,gt,1e-8_dp,2,expected,reference_flux,b,w,good)
      call require(good,'serial MOL rate')
      call require(maxval(abs(rates%fields(p)%values-expected)/max(1._dp,abs(expected)))<1e-12_dp, &
        'serial distributed MOL rate agreement')
      call require(maxval(abs(flux(p)%x-reference_flux%x)/max(1._dp,abs(reference_flux%x)))<1e-12_dp.and. &
        maxval(abs(flux(p)%y-reference_flux%y)/max(1._dp,abs(reference_flux%y)))<1e-12_dp.and. &
        maxval(abs(flux(p)%z-reference_flux%z)/max(1._dp,abs(reference_flux%z)))<1e-12_dp, &
        'serial distributed MOL flux agreement')
    end do
    prior=rates;previous=flux;wrong_order=2
    if(rank==nranks-1) wrong_order=0
    call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,wrong_order, &
      rates,flux,MPI_COMM_WORLD,good,boundary_faces=labels)
    call require(.not.good,'collective MOL invalid control')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(rates%fields(p)%values==prior%fields(p)%values),'MOL rate rollback')
      call require(all(flux(p)%x==previous(p)%x).and.all(flux(p)%y==previous(p)%y).and. &
        all(flux(p)%z==previous(p)%z),'MOL flux rollback')
    end do
    call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,local,lt,1e-8_dp,2, &
      rates,flux,MPI_COMM_WORLD,good,boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t, &
      lax_friedrichs=.true.)
    call require(good,'distributed LF MOL rate and peer coefficient exchange')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(ieee_is_finite(rates%fields(p)%values)),'distributed LF MOL finite rates')
    end do
    h2=0;n2=0
    do p=1,size(species)
      if(species(p)%name=='H2') h2=p
      if(species(p)%name=='N2') n2=p
    end do
    call require(h2>0.and.n2>0,'LF peer species fixture')
    lf_tree=stage_tree
    do p=1,count
      do z=1,4
      do y=1,4
      do x=1,4
        if(physical(p)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,lf_tree%nodes(p)%state(:,x,y,z), &
          lf_tree%nodes(p)%temperature(x,y,z),qp,temp_cell,speed,good)
        call require(good,'LF peer primitive')
        qp(6:)=0;qp(n2+5)=1
        if(p==3.and.z==2) then
          qp(2)=-10000._dp;qp(h2+5)=.2_dp;qp(n2+5)=.8_dp
        end if
        call reactive_primitive_to_conserved(species,qp,lf_tree%nodes(p)%state(:,x,y,z), &
          lf_tree%nodes(p)%temperature(x,y,z),speed,good)
        call require(good,'LF peer conserved state')
      end do
      end do
      end do
    end do
    lf_local=local;lf_temperature=lt
    do p=1,count
      if(owners(p)/=rank) cycle
      lf_local%fields(p)%values=lf_tree%nodes(p)%state
      lf_temperature%fields(p)%values(1,:,:,:)=lf_tree%nodes(p)%temperature
    end do
    call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,lf_local,lf_temperature, &
      1e-8_dp,2,rates,flux,MPI_COMM_WORLD,good,lax_friedrichs=.true.)
    call require(good,'LF cut peer boundary coefficient exchange')
    call compute_mol_eb_tree_hydro_3d(species,lf_tree,three,1e-8_dp,2,serial_rates,serial_flux,good, &
      lax_friedrichs=.true.)
    call require(good,'serial LF cut peer boundary coefficient exchange')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(maxval(abs(serial_rates%fields(p)%values-rates%fields(p)%values)/ &
        max(1._dp,abs(rates%fields(p)%values)))<1e-12_dp,'serial and MPI LF peer rates')
      call require(maxval(abs(serial_flux(p)%x-flux(p)%x)/max(1._dp,abs(flux(p)%x)))<1e-12_dp.and. &
        maxval(abs(serial_flux(p)%y-flux(p)%y)/max(1._dp,abs(flux(p)%y)))<1e-12_dp.and. &
        maxval(abs(serial_flux(p)%z-flux(p)%z)/max(1._dp,abs(flux(p)%z)))<1e-12_dp, &
        'serial and MPI LF peer fluxes')
    end do
    left_face=0;right_face=0
    if(owners(2)==rank) left_face=flux(2)%x(:,4,:,:)
    if(owners(3)==rank) right_face=flux(3)%x(:,0,:,:)
    call MPI_Allreduce(MPI_IN_PLACE,left_face,size(left_face),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call MPI_Allreduce(MPI_IN_PLACE,right_face,size(right_face),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
    call require(maxval(abs(left_face-right_face)/max(1._dp,abs(left_face)))<1e-12_dp, &
      'LF cut peer boundary shares one limited face flux')
    changed=0
    do p=1,count
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,lf_tree,p,ghost,gt,good, &
        halo_width=3,target_volume_fraction=three(p)%volume_fraction,ghost_volume_fraction=vf, &
        restrict_snapshot=.false.)
      call require(good,'LF cut peer serial donor')
      call compute_mol_eb_hydro_3d(species,three(p),ghost,gt,1e-8_dp,2,expected,reference_flux,b,w,good, &
        lax_friedrichs=.true.)
      call require(good,'LF cut peer serial spatial rate')
      if(p==2) then
        if(any(flux(p)%x(:,4,:,:)/=reference_flux%x(:,4,:,:))) changed=changed+1
      else if(p==3) then
        if(any(flux(p)%x(:,0,:,:)/=reference_flux%x(:,0,:,:))) changed=changed+1
      end if
      cell_budget=physical(p)%dx*physical(p)%dy*physical(p)%dz* &
        sum(physical(p)%volume_fraction*rates%fields(p)%values(h2+5,:,:,:))*1e-8_dp
      face_budget=physical(p)%dy*physical(p)%dz*sum( &
        physical(p)%x_face_fraction(0,:,:)*flux(p)%x(h2+5,0,:,:)- &
        physical(p)%x_face_fraction(4,:,:)*flux(p)%x(h2+5,4,:,:))+ &
        physical(p)%dx*physical(p)%dz*sum( &
        physical(p)%y_face_fraction(:,0,:)*flux(p)%y(h2+5,:,0,:)- &
        physical(p)%y_face_fraction(:,4,:)*flux(p)%y(h2+5,:,4,:))+ &
        physical(p)%dx*physical(p)%dy*sum( &
        physical(p)%z_face_fraction(:,:,0)*flux(p)%z(h2+5,:,:,0)- &
        physical(p)%z_face_fraction(:,:,4)*flux(p)%z(h2+5,:,:,4))
      call require(abs(cell_budget-face_budget)<1e-13_dp,'LF peer face species budget')
      do z=1,4
      do y=1,4
      do x=1,4
        if(physical(p)%volume_fraction(x,y,z)<=0) cycle
        call require(lf_local%fields(p)%values(h2+5,x,y,z)+ &
          1e-8_dp*rates%fields(p)%values(h2+5,x,y,z)>=0._dp,'LF peer forward species state')
      end do
      end do
      end do
    end do
    call MPI_Allreduce(changed,global_changed,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(global_changed>0,'LF peer face limiter activates on cut boundary')
    if(axis==3.and.count==5) then
      seam_faces='outflow';seam_faces(1:4)='periodic'
      do z=1,4
      do y=1,4
      do x=1,4
        if(physical(3)%volume_fraction(x,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,lf_tree%nodes(3)%state(:,x,y,z), &
          lf_tree%nodes(3)%temperature(x,y,z),qp,temp_cell,speed,good)
        call require(good,'LF periodic reset primitive')
        qp(2)=0._dp
        call reactive_primitive_to_conserved(species,qp,lf_tree%nodes(3)%state(:,x,y,z), &
          lf_tree%nodes(3)%temperature(x,y,z),speed,good)
        call require(good,'LF periodic reset conserved state')
      end do
      end do
      end do
      do z=1,4
      do y=1,4
        if(physical(1)%volume_fraction(4,y,z)<=0) cycle
        call reactive_conserved_to_primitive(species,lf_tree%nodes(1)%state(:,4,y,z), &
          lf_tree%nodes(1)%temperature(4,y,z),qp,temp_cell,speed,good)
        call require(good,'LF periodic primitive')
        qp(6:)=0;qp(n2+5)=1
        if(z==2) then
          qp(h2+5)=.2_dp;qp(n2+5)=.8_dp
        end if
        call reactive_primitive_to_conserved(species,qp,lf_tree%nodes(1)%state(:,4,y,z), &
          lf_tree%nodes(1)%temperature(4,y,z),speed,good)
        call require(good,'LF periodic conserved state')
      end do
      end do
      do p=1,count
        if(owners(p)/=rank) cycle
        lf_local%fields(p)%values=lf_tree%nodes(p)%state
        lf_temperature%fields(p)%values(1,:,:,:)=lf_tree%nodes(p)%temperature
      end do
      call compute_mol_mpi_eb_tree_hydro_3d(species,topology,owners,physical,three,lf_local,lf_temperature, &
        1e-8_dp,2,rates,flux,MPI_COMM_WORLD,good,lax_friedrichs=.true.,boundary_faces=seam_faces)
      call require(good,'LF periodic MPI spatial rate')
      call compute_mol_eb_tree_hydro_3d(species,lf_tree,three,1e-8_dp,2,serial_rates,serial_flux,good, &
        lax_friedrichs=.true.,boundary_faces=seam_faces)
      call require(good,'LF periodic serial spatial rate')
      call fill_reactive_eb_tree_halo_3d(species,lf_tree,1,ghost,gt,good, &
        halo_width=3,target_volume_fraction=three(1)%volume_fraction,ghost_volume_fraction=vf, &
        boundary_faces=seam_faces,restrict_snapshot=.false.)
      call require(good,'LF periodic serial root donors')
      call compute_mol_eb_hydro_3d(species,three(1),ghost,gt,1e-8_dp,2,expected,reference_flux,b,w,good, &
        lax_friedrichs=.true.)
      call require(good,'LF periodic unbounded root face')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(maxval(abs(serial_rates%fields(p)%values-rates%fields(p)%values)/ &
          max(1._dp,abs(rates%fields(p)%values)))<1e-12_dp,'LF periodic serial/MPI rates')
        call require(maxval(abs(serial_flux(p)%x-flux(p)%x)/max(1._dp,abs(flux(p)%x)))<1e-12_dp.and. &
          maxval(abs(serial_flux(p)%y-flux(p)%y)/max(1._dp,abs(flux(p)%y)))<1e-12_dp.and. &
          maxval(abs(serial_flux(p)%z-flux(p)%z)/max(1._dp,abs(flux(p)%z)))<1e-12_dp, &
          'LF periodic serial/MPI fluxes')
        do z=1,4
        do y=1,4
        do x=1,4
          if(physical(p)%volume_fraction(x,y,z)<=0) cycle
          call require(lf_local%fields(p)%values(h2+5,x,y,z)+ &
            1e-8_dp*rates%fields(p)%values(h2+5,x,y,z)>=0._dp,'LF periodic forward fuel')
        end do
        end do
        end do
      end do
      left_face=0;right_face=0
      if(owners(1)==rank) then
        left_face=flux(1)%x(:,0,:,:)
        right_face=flux(1)%x(:,4,:,:)
      end if
      call MPI_Allreduce(MPI_IN_PLACE,left_face,size(left_face),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call MPI_Allreduce(MPI_IN_PLACE,right_face,size(right_face),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,status)
      call require(maxval(abs(left_face-right_face)/max(1._dp,abs(right_face)))<1e-12_dp, &
        'LF periodic root seam shares flux')
      candidate=lf_local;candidate_temperature=lf_temperature
      do p=1,count
        if(owners(p)/=rank) cycle
        candidate%fields(p)%values=candidate%fields(p)%values+1e-8_dp*rates%fields(p)%values
      end do
      call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,physical,candidate, &
        candidate_temperature,flux,MPI_COMM_WORLD,good,peers_only=.true.,boundary_faces=seam_faces)
      call require(good,'LF periodic peer and seam flux synchronization')
      do p=1,count
        if(owners(p)/=rank) cycle
        do z=1,4
        do y=1,4
        do x=1,4
          if(physical(p)%volume_fraction(x,y,z)<=0) cycle
          call require(candidate%fields(p)%values(h2+5,x,y,z)>=0._dp, &
            'LF periodic fuel after peer and seam synchronization')
        end do
        end do
        end do
      end do
      prior=candidate;saved_temperature=candidate_temperature;previous=flux
      if(test_mode=='mol_reflux') then
        call assemble_mpi_eb_tree_reflux_3d(species,topology,owners,physical,prior,saved_temperature,previous, &
          signed_increment,audited_flux,MPI_COMM_WORLD,good,boundary_faces=seam_faces)
        call require(good,'LF signed coarse/fine increment assembly')
        signed_endpoint=prior
        negative_endpoints=0
        do p=1,count
          if(owners(p)/=rank) cycle
          signed_endpoint%fields(p)%values=signed_endpoint%fields(p)%values+signed_increment%fields(p)%values
          do z=1,4
            do y=1,4
              do x=1,4
                if(physical(p)%volume_fraction(x,y,z)<=0.or.covered_amr(topology,p,[x,y,z])) cycle
                if(any(signed_endpoint%fields(p)%values(6:,x,y,z)<0)) negative_endpoints=negative_endpoints+1
              enddo
            enddo
          enddo
        enddo
        call MPI_Allreduce(negative_endpoints,global_negative,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
        call require(global_negative>0,'LF fixture exhausts exposed coarse/fine species budget')
        call integrate_mpi_eb_tree_3d(topology,owners,physical,signed_endpoint,signed_inventory,MPI_COMM_WORLD,good)
        call require(good,'LF unrestricted signed composite inventory')
      endif
      call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,physical,candidate, &
        candidate_temperature,flux,MPI_COMM_WORLD,good,boundary_faces=seam_faces)
      if(test_mode=='mol_reflux') then
        repairs=0
        if(.not.good) then
          do depth=topology%level_count()-2,0,-1
            call repair_mpi_eb_tree_interface_budget_3d(species,topology,owners,physical,candidate, &
              candidate_temperature,flux,MPI_COMM_WORLD,good,depth,boundary_faces=seam_faces)
            call require(good,'LF conservative interface budget repair')
            repairs=repairs+1
            call synchronize_mpi_eb_tree_flux_3d(species,topology,owners,physical,candidate, &
              candidate_temperature,flux,MPI_COMM_WORLD,good,boundary_faces=seam_faces)
            if(good) exit
          enddo
        endif
        call require(good,'LF stressed coarse/fine endpoint must be accepted')
        call require(repairs>0,'LF coarse/fine budget repair activates')
        call integrate_mpi_eb_tree_3d(topology,owners,physical,candidate,accepted_inventory,MPI_COMM_WORLD,good)
        call require(good,'LF repaired composite inventory')
        call require(all(abs(accepted_inventory-signed_inventory)<= &
          128*epsilon(1._dp)*max(1._dp,abs(signed_inventory))),'LF repaired composite vector conservation')
        inventory_residual=maxval(abs(accepted_inventory-signed_inventory)/max(1._dp,abs(signed_inventory)))
        face_residual=0
        do p=1,count
          if(owners(p)/=rank) cycle
          do z=1,4
            do y=1,4
              do x=1,4
                if(physical(p)%volume_fraction(x,y,z)<=0.or.covered_amr(topology,p,[x,y,z])) cycle
                face_change=(physical(p)%x_face_fraction(x-1,y,z)*(flux(p)%x(:,x-1,y,z)-previous(p)%x(:,x-1,y,z))- &
                  physical(p)%x_face_fraction(x,y,z)*(flux(p)%x(:,x,y,z)-previous(p)%x(:,x,y,z)))/physical(p)%dx+ &
                  (physical(p)%y_face_fraction(x,y-1,z)*(flux(p)%y(:,x,y-1,z)-previous(p)%y(:,x,y-1,z))- &
                  physical(p)%y_face_fraction(x,y,z)*(flux(p)%y(:,x,y,z)-previous(p)%y(:,x,y,z)))/physical(p)%dy+ &
                  (physical(p)%z_face_fraction(x,y,z-1)*(flux(p)%z(:,x,y,z-1)-previous(p)%z(:,x,y,z-1))- &
                  physical(p)%z_face_fraction(x,y,z)*(flux(p)%z(:,x,y,z)-previous(p)%z(:,x,y,z)))/physical(p)%dz
                face_change=face_change/physical(p)%volume_fraction(x,y,z)
                observed_change=candidate%fields(p)%values(:,x,y,z)-prior%fields(p)%values(:,x,y,z)
                face_residual=max(face_residual,maxval(abs(observed_change-face_change)/ &
                  max(1._dp,abs(prior%fields(p)%values(:,x,y,z)))))
                call require(all(abs(observed_change-face_change)<=128*epsilon(1._dp)* &
                  max(1._dp,abs(prior%fields(p)%values(:,x,y,z)))),'LF independent repaired face divergence')
              enddo
            enddo
          enddo
        enddo
        call MPI_Allreduce(face_residual,global_face_residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
        if(rank==0) write(*,'(a,i0,a,es12.4,a,es12.4)') 'LF reflux repair interfaces: ',repairs, &
          ' inventory relative: ',inventory_residual,' face divergence relative: ',global_face_residual
        good=.true.
      endif
      do p=1,count
        if(owners(p)/=rank) cycle
        if(good) then
          call require(all(candidate%fields(p)%values(6:,:,:,:)>=0._dp), &
            'LF accepted coarse/fine species remain nonnegative')
        else
          call require(all(candidate%fields(p)%values==prior%fields(p)%values).and. &
            all(candidate_temperature%fields(p)%values==saved_temperature%fields(p)%values), &
            'LF rejected coarse/fine state and temperature rollback')
          call require(all(flux(p)%x==previous(p)%x).and.all(flux(p)%y==previous(p)%y).and. &
            all(flux(p)%z==previous(p)%z),'LF rejected coarse/fine flux rollback')
        end if
      end do
    end if
    if(rank==0) print *, 'PASS: MPI MOL source and flux comparison, collective rollback'
  end subroutine check_mol_sources
  subroutine check_temporal_sources(owners,predictor,local_predictor)
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(in) :: predictor,local_predictor
    type(amr_tree_state_3d) :: first,last,local_first,local_last,halo,previous
    real(dp) :: fractions(2),reference(nv,10,10,10),nan
    integer :: p,c,phase,depth
    logical :: good
    nan=ieee_value(0._dp,ieee_quiet_nan)
    do phase=1,3
      fractions=real(phase-1,dp)/2
      first=predictor;last=predictor;local_first=local_predictor;local_last=local_predictor
      do p=1,count
        do c=1,nv
          where(physical(p)%volume_fraction>0)
            first%fields(p)%values(c,:,:,:)=predictor%fields(p)%values(c,:,:,:)-real(p,dp)
            last%fields(p)%values(c,:,:,:)=predictor%fields(p)%values(c,:,:,:)+3*real(p,dp)
          end where
        end do
        if(phase==1) last%fields(p)%values=nan
        if(phase==3) first%fields(p)%values=nan
        if(owners(p)/=rank) cycle
        local_first%fields(p)%values=first%fields(p)%values
        local_last%fields(p)%values=last%fields(p)%values
      end do
      call fill_sdc_mpi_eb_tree_source_halo_3d(topology,owners,physical,local_predictor,nv,3,halo, &
        MPI_COMM_WORLD,good,boundary_faces=labels,ancestor_start=local_first,ancestor_end=local_last, &
        ancestor_fraction=fractions)
      call require(good,'MPI temporal source halo, unused endpoint NaNs ignored')
      do p=1,count
        call require(allocated(halo%fields(p)%values).eqv.(owners(p)==rank),'owner-local temporal source halo')
        if(owners(p)/=rank) cycle
        depth=topology%nodes(p)%level
        call fill_sdc_eb_tree_source_halo_3d(saved,predictor,p,3,reference,good,boundary_faces=labels, &
          ancestor_start=first,ancestor_end=last,ancestor_fraction=fractions(:depth))
        call require(good,'serial temporal source halo')
        call require(all(halo%fields(p)%values==reference),'exact serial temporal source halo')
      end do
    end do
    previous=halo
    fractions=.5_dp
    if(rank==nranks-1) fractions(1)=.25_dp
    call fill_sdc_mpi_eb_tree_source_halo_3d(topology,owners,physical,local_predictor,nv,3,halo, &
      MPI_COMM_WORLD,good,ancestor_start=local_first,ancestor_end=local_last,ancestor_fraction=fractions)
    call require(.not.good,'rank-dependent source fractions reject collectively')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(halo%fields(p)%values==previous%fields(p)%values),'temporal source halo rollback')
    end do
    fractions=.5_dp;local_first=local_predictor;local_last=local_predictor
    if(rank==owners(1)) local_last%fields(1)%values=nan
    call fill_sdc_mpi_eb_tree_source_halo_3d(topology,owners,physical,local_predictor,nv,3,halo, &
      MPI_COMM_WORLD,good,ancestor_start=local_first,ancestor_end=local_last,ancestor_fraction=fractions)
    call require(.not.good,'used ancestor source NaN rejects collectively')
    do p=1,count
      if(owners(p)/=rank) cycle
      call require(all(halo%fields(p)%values==previous%fields(p)%values),'failed source fetch rollback')
    end do
  end subroutine check_temporal_sources
  subroutine check_sdc_hydro()
    type(eb_reactive_tree_3d) :: stage_tree,end_tree
    type(amr_tree_state_3d) :: local,lt,predictor,local_predictor,rates,reference,prior,provided,source_metrics
    type(amr_tree_state_3d) :: particle_peers,particle_remainder,composed,saved_composed
    type(eb_integrated_flux_3d), allocatable :: flux(:),expected(:),old_flux(:),baseline(:),saved_baseline(:)
    type(eb_integrated_flux_3d), allocatable :: expected_baseline(:)
    type(amr_tree_state_3d) :: end_local,end_source,local_end_source,level_rates,level_reference
    type(eb_integrated_flux_3d), allocatable :: level_flux(:),level_expected(:)
    real(dp) :: dt,nan
    integer :: owners(count),p,mode,failure,c,halo_status,x,y,z
    logical :: good,changed,global_changed
    stage_tree=saved;nan=ieee_value(0._dp,ieee_quiet_nan)
    stage_tree%nodes(1)%state=1.01_dp*stage_tree%nodes(1)%state
    allocate(predictor%fields(count))
    do p=1,count
      allocate(predictor%fields(p)%values(nv,4,4,4),source=0._dp)
      predictor%fields(p)%values(2,:,:,:)=real(p,dp)*2.e4_dp
      predictor%fields(p)%values(3,:,:,:)=-real(p,dp)*3.e4_dp
      predictor%fields(p)%values(5,:,:,:)=real(p,dp)*1.e6_dp
      do c=1,nv
        where(physical(p)%volume_fraction<=0) predictor%fields(p)%values(c,:,:,:)=nan
      end do
    end do
    call compute_sdc_eb_tree_hydro_3d(species,stage_tree,predictor,extended,1.e-8_dp,reference,expected,good, &
      difmag=coefficient,boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t, &
      baseline_flux=expected_baseline)
    call require(good,'serial source-aware SDC hydro')
    do mode=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(mode==2) owners=0
      allocate(local%fields(count),lt%fields(count),local_predictor%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        local%fields(p)%values=stage_tree%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=stage_tree%nodes(p)%temperature
        local_predictor%fields(p)%values=predictor%fields(p)%values
      end do
      call check_mol_sources(stage_tree,owners,local,lt)
      call check_temporal_sources(owners,predictor,local_predictor)
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,baseline_flux=baseline)
      call require(good,'MPI SDC hydro with covered NaNs')
      do p=1,count
        call require(allocated(rates%fields(p)%values).eqv.(owners(p)==rank),'owner-local SDC rates')
        call require(allocated(flux(p)%x).eqv.(owners(p)==rank),'owner-local SDC flux')
        call require(allocated(baseline(p)%x).eqv.(owners(p)==rank),'owner-local SDC baseline')
        if(owners(p)/=rank) cycle
        call require(all(rates%fields(p)%values==reference%fields(p)%values),'exact serial SDC hydro rate')
        call require(all(flux(p)%x==expected(p)%x).and.all(flux(p)%y==expected(p)%y).and. &
          all(flux(p)%z==expected(p)%z),'exact serial SDC hydro flux')
        call require(all(local%fields(p)%values==stage_tree%nodes(p)%state),'SDC preserves input U')
        call require(all(shape(baseline(p)%x)==shape(flux(p)%x)).and. &
          all(shape(baseline(p)%y)==shape(flux(p)%y)).and.all(shape(baseline(p)%z)==shape(flux(p)%z)), &
          'SDC baseline has physical face shapes')
        call require(all(ieee_is_finite(baseline(p)%x)).and.all(ieee_is_finite(baseline(p)%y)).and. &
          all(ieee_is_finite(baseline(p)%z)), 'SDC first-order baseline is finite')
        call require(all(baseline(p)%x==expected_baseline(p)%x).and.all(baseline(p)%y==expected_baseline(p)%y).and. &
          all(baseline(p)%z==expected_baseline(p)%z), 'exact serial/MPI first-order baseline flux')
        call require(all(lt%fields(p)%values(1,:,:,:)==stage_tree%nodes(p)%temperature),'SDC preserves input T')
      end do
      saved_baseline=baseline
      allocate(source_metrics%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        allocate(source_metrics%fields(p)%values(1,10,10,10))
        source_metrics%fields(p)%values(1,:,:,:)=extended(p)%volume_fraction(4:13,4:13,4:13)
      end do
      call fill_sdc_mpi_eb_tree_source_halo_3d(topology,owners,physical,local_predictor,nv,3,provided, &
        MPI_COMM_WORLD,good,boundary_faces=labels,target_volume_fraction=source_metrics)
      call require(good,'prepare explicit predictor halo')
      call fill_sdc_mpi_eb_tree_source_halo_3d(topology,owners,physical,local_predictor,nv,3,particle_peers, &
        MPI_COMM_WORLD,good,boundary_faces=labels,target_volume_fraction=source_metrics,same_level_only=.true.)
      call require(good,'particle source peer-only fill')
      particle_remainder=particle_peers;changed=.false.
      do p=1,count
        if(owners(p)/=rank) cycle
        particle_remainder%fields(p)%values=0
        changed=changed.or.any(particle_peers%fields(p)%values/=provided%fields(p)%values)
        do z=1,10
        do y=1,10
        do x=1,10
          if(all([x,y,z]>=4).and.all([x,y,z]<=7)) cycle
          if(source_metrics%fields(p)%values(1,x,y,z)<=0) cycle
          if(all(particle_peers%fields(p)%values(:,x,y,z)==0)) particle_remainder%fields(p)%values(2,x,y,z)=7
        end do
        end do
        end do
      end do
      call MPI_Allreduce(changed,global_changed,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,halo_status)
      call require(global_changed,'peer-only source excludes physical/ancestor extrapolation')
      call fill_sdc_mpi_eb_particle_predictor_halo_3d(topology,owners,physical,local_predictor,local_predictor, &
        particle_remainder,nv,3,composed,MPI_COMM_WORLD,good,labels,source_metrics, &
        local_predictor,local_predictor,local_predictor,local_predictor,[.5_dp,.5_dp])
      call require(good,'compose complete predictor with deposited particle halo')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(composed%fields(p)%values==particle_peers%fields(p)%values+ &
          particle_remainder%fields(p)%values),'particle-only predictor has no duplicated ancestor contribution')
      end do
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,predictor_halo=composed,baseline_flux=baseline)
      call require(good,'composed particle predictor reaches hydro')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(baseline(p)%x==saved_baseline(p)%x).and.all(baseline(p)%y==saved_baseline(p)%y).and. &
          all(baseline(p)%z==saved_baseline(p)%z),'first-order baseline is independent of transverse predictor sources')
      enddo
      saved_composed=composed
      do p=1,count
        if(owners(p)==rank) particle_remainder%fields(p)%values(2,4:7,4:7,4:7)=1
      end do
      call fill_sdc_mpi_eb_particle_predictor_halo_3d(topology,owners,physical,local_predictor,local_predictor, &
        particle_remainder,nv,3,composed,MPI_COMM_WORLD,good,labels,source_metrics)
      call require(.not.good,'particle remainder cannot replace valid-cell sources')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(composed%fields(p)%values==saved_composed%fields(p)%values),'particle predictor failure rollback')
      end do

      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,predictor_halo=provided)
      call require(good,'explicit predictor halo accepted')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(rates%fields(p)%values==reference%fields(p)%values),'explicit halo preserves baseline rates')
        provided%fields(p)%values(2,1:3,:,:)=provided%fields(p)%values(2,1:3,:,:)+1e5_dp
      end do
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,predictor_halo=provided)
      call require(good,'modified ghost predictor accepted')
      changed=.false.
      do p=1,count
        if(owners(p)/=rank) cycle
        changed=changed.or.any(rates%fields(p)%values/=reference%fields(p)%values)
      end do
      call MPI_Allreduce(changed,global_changed,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,halo_status)
      call require(global_changed,'ghost predictor contributes to hydro flux divergence')
      prior=rates;old_flux=flux
      do p=1,count
        if(owners(p)==rank) provided%fields(p)%values(:,4:7,4:7,4:7)=provided%fields(p)%values(:,4:7,4:7,4:7)+1
      end do
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,predictor_halo=provided)
      call require(.not.good,'inconsistent valid predictor rejected')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(rates%fields(p)%values==prior%fields(p)%values),'explicit predictor rate rollback')
        call require(all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
          all(flux(p)%z==old_flux(p)%z),'explicit predictor flux rollback')
      end do
      deallocate(source_metrics%fields)
      end_tree=stage_tree;end_local=local;end_source=predictor;local_end_source=local_predictor
      do p=1,count
        end_tree%nodes(p)%state=1.01_dp*stage_tree%nodes(p)%state
        do c=1,nv
          where(physical(p)%volume_fraction>0) end_source%fields(p)%values(c,:,:,:)= &
            predictor%fields(p)%values(c,:,:,:)+real(p,dp)
        end do
        if(owners(p)/=rank) cycle
        end_local%fields(p)%values=end_tree%nodes(p)%state
        local_end_source%fields(p)%values=end_source%fields(p)%values
      end do
      call compute_sdc_eb_tree_hydro_3d(species,stage_tree,predictor,extended,1.e-8_dp, &
        level_reference,level_expected,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,target_depth=2,ancestor_start=stage_tree,ancestor_end=end_tree, &
        ancestor_fraction=[.5_dp,.5_dp],source_start=predictor,source_end=end_source,source_fraction=[.25_dp,.25_dp])
      call require(good,'serial selected SDC hydro')
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,level_rates,level_flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t,target_depth=2,ancestor_start=local,ancestor_end=end_local, &
        temperature_start=lt,temperature_end=lt,ancestor_fraction=[.5_dp,.5_dp], &
        source_start=local_predictor,source_end=local_end_source,source_fraction=[.25_dp,.25_dp])
      call require(good,'MPI selected SDC hydro')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(level_rates%fields(p)%values==level_reference%fields(p)%values),'selected MPI hydro rate')
        if(topology%nodes(p)%level/=2) then
          call require(.not.allocated(level_flux(p)%x),'inactive hydro flux absent')
          cycle
        end if
        call require(all(level_flux(p)%x==level_expected(p)%x).and.all(level_flux(p)%y==level_expected(p)%y).and. &
          all(level_flux(p)%z==level_expected(p)%z),'selected MPI hydro flux')
      end do
      prior=rates;old_flux=flux
      do failure=1,2
        dt=1.e-8_dp
        if(failure==1.and.rank==owners(count)) local_predictor%fields(count)%values(2,4,4,4)=nan
        if(failure==2.and.rank==nranks-1) dt=-dt
        call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
          extended,dt,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
          inflow_state=inlet,inflow_temperature=inlet_t)
        call require(.not.good,'collective SDC hydro rejection')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(rates%fields(p)%values==prior%fields(p)%values),'SDC hydro rates rollback')
          call require(all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
            all(flux(p)%z==old_flux(p)%z),'SDC hydro flux rollback')
        end do
        if(rank==owners(count)) local_predictor%fields(count)%values=predictor%fields(count)%values
      end do
      deallocate(local%fields,lt%fields,local_predictor%fields)
    end do
  end subroutine check_sdc_hydro
  subroutine check_mpi_hydro()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    type(eb_geometry_3d), allocatable :: trial_extended(:)
    real(dp), allocatable :: mpi_coefficient
    real(dp) :: b(nv),w(nv),inventory(nv),error,dt,local_coefficient,inlet_density
    integer :: owners(count),p,mode,failure
    logical :: good
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
      initial=local;it=lt;b=17._dp;w=19._dp;local_coefficient=coefficient
      if(test_mode/='default') mpi_coefficient=coefficient
      call advance_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,solver,1.e-8_dp, &
        local,lt,b,w,MPI_COMM_WORLD,good,difmag=mpi_coefficient,extended_geometry=extended, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
      call require(good,'MPI hydro')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local U')
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial hydro U')
        call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial hydro T')
      end do
      call require(all(b==boundary),'exact boundary ledger')
      call require(maxval(abs(w-wall)/max(1._dp,abs(wall)))<3.e-14_dp,'wall ledger')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
      call require(good,'MPI inventory')
      error=maxval(abs(inventory-before-b-w)/max(1._dp,abs(before)))
      call require(error<3.e-14_dp,'MPI hydro conservation')
      if (rank==0) print *, 'MPI axis, nodes, mode, residual: ',axes(axis),count,mode,error
      do failure=1,10
        if(failure==10.and.index(test_mode,'periodic')==0.and.index(test_mode,'wall')==0) cycle
        if(failure==9.and..not.allocated(inlet)) cycle
        if(failure>=7.and.index(test_mode,'ctu')/=1) cycle
        if(allocated(extended)) trial_extended=extended
        if(failure==7.and.rank==0) deallocate(trial_extended)
        if(failure==8.and.rank==0) trial_extended(1)%x_face_centroid_y(0,1,1)=.1_dp
        if(failure==4.and.coefficient==0) cycle
        do p=1,count
          geometry(p)=saved%nodes(p)%geometry
        end do
        if(failure==5.and.rank==0) geometry(1)%halo_volume_fraction(0,0,0)=.5_dp
        if(failure==6.and.rank==0) deallocate(geometry(1)%halo_volume_fraction)
        local_coefficient=coefficient
        if(failure==4.and.rank==0) local_coefficient=.05_dp
        local=initial;lt=it;b=17._dp;w=19._dp;dt=1.e-8_dp
        if(failure==10.and.rank==0) labels(:2)='outflow'
        if(failure==9.and.rank==0) then
          inlet_density=inlet(1,1);inlet(1,1)=inlet_density+.01_dp
        end if
        if (failure==1) then
          call advance_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,'invalid',dt, &
            local,lt,b,w,MPI_COMM_WORLD,good,difmag=local_coefficient,extended_geometry=trial_extended, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
        else
          if (failure==2.and.rank==0) dt=2.e-8_dp
          if (failure==3.and.rank==owners(count)) local%fields(count)%values(1,4,4,4)=-1._dp
          call advance_mpi_eb_tree_hydro_3d(species,topology,owners,geometry,solver,dt, &
            local,lt,b,w,MPI_COMM_WORLD,good,difmag=local_coefficient,extended_geometry=trial_extended, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
        end if
        if(failure==10.and.rank==0) then
          labels(:2)='periodic'
          if(index(test_mode,'wall')>0) labels(:2)='slip_wall'
        end if
        if(failure==9.and.rank==0) inlet(1,1)=inlet_density
        call require(.not.good,'collective hydro rejection')
        if (failure==3.and.rank==owners(count)) then
          call require(local%fields(count)%values(1,4,4,4)==-1._dp,'bad density preserved')
          local%fields(count)%values(1,4,4,4)=initial%fields(count)%values(1,4,4,4)
        end if
        do p=1,count
          if (owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==initial%fields(p)%values),'U rollback')
          call require(all(lt%fields(p)%values==it%fields(p)%values),'T rollback')
        end do
        call require(all(b==17._dp).and.all(w==19._dp),'ledger rollback')
      end do
      do p=1,count
        geometry(p)=saved%nodes(p)%geometry
      end do
      deallocate(local%fields,lt%fields)
    end do
  end subroutine check_mpi_hydro
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_hydro_3d
