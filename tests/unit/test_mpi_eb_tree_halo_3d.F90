program test_mpi_eb_tree_halo_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved,reactive_pelec_flux_x
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    fill_reactive_eb_tree_halo_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: scatter_mpi_amr_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,geometry(5),original(5)
  type(eb_reactive_tree_3d) :: tree
  type(reactive_eb_3d_config) :: config
  type(amr_tree_state_3d) :: global,gt,local,lt,ghost,gh_t,saved,st,saved_ghost,saved_gt
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: reference(:,:,:,:),reference_t(:,:,:)
  integer :: nv,node,count,rank,nranks,status,owners(5),root,mode
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  root=nranks-1;count=5
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
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
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  config%initial_velocity_x=1;config%initial_velocity_y=2;config%initial_velocity_z=3
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')

  do node=2,5
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  do node=1,5
    geometry(node)=tree%nodes(node)%geometry
  end do
  original=geometry
  if (rank==root) then
    allocate(global%fields(5),gt%fields(5))
    do node=1,5
      global%fields(node)%values=tree%nodes(node)%state
      allocate(gt%fields(node)%values(1,4,4,4))
      gt%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
    end do
  end if
  allocate(reference(nv,0:5,0:5,0:5),reference_t(0:5,0:5,0:5))
  do mode=1,2
    owners=[(mod(node-1,nranks),node=1,5)]
    if (mode==2) owners=0
    call scatter_mpi_amr_tree_3d(topology,owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter U')
    call scatter_mpi_amr_tree_3d(topology,owners,gt,lt,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter T')
    saved=local;st=lt
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    call require(ok,'MPI halos')
    call unchanged_inputs()
    do node=1,5
      call require(allocated(ghost%fields(node)%values).eqv.(rank==owners(node)),'owner-local halo U')
      call require(allocated(gh_t%fields(node)%values).eqv.(rank==owners(node)),'owner-local halo T')
      if (rank/=owners(node)) cycle
      call require(all(lbound(ghost%fields(node)%values)==[1,0,0,0]),'halo lower bounds')
      call require(all(lbound(gh_t%fields(node)%values)==[1,0,0,0]),'temperature halo lower bounds')
      call fill_reactive_eb_tree_halo_3d(species,tree,node,reference,reference_t,ok)
      call require(ok,'serial halos')
      call require(all(ghost%fields(node)%values==reference),'exact serial face U and inactive zeros')
      call require(all(gh_t%fields(node)%values(1,:,:,:)==reference_t),'exact serial face T')
    end do
    call temperature_seed_checks()
    call separate_momentum_checks()
    call paired_pressure_checks()
    call wide_checks()
    call wide_checks(.true.)
    call prescribed_checks()
    call periodic_checks()
    call periodic_checks(walls=.true.)
    call wide_checks(walls=.true.)
    call wide_checks(periodic=.true.)
    call temporal_checks()
    ! Root unrefined boundary survives averaging, then fails donor EOS after exchange.
    if (rank==owners(1)) local%fields(1)%values(1,4,4,4)=-1._dp
    saved=local;saved_ghost=ghost;saved_gt=gh_t
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    call require(.not.ok,'bad donor EOS rejected collectively')
    call unchanged_inputs()
    call unchanged_halos()
    if (rank==owners(1)) local%fields(1)%values(1,4,4,4)=tree%nodes(1)%state(1,4,4,4)
    saved=local
    if (rank==0) then
      do node=1,5
        geometry(node)%x_lower=geometry(node)%x_lower+.01_dp
        geometry(node)%x_upper=geometry(node)%x_upper+.01_dp
      end do
    end if
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    call require(.not.ok,'geometry mismatch rejected collectively')
    call unchanged_inputs()
    call unchanged_halos()
    geometry=original
  end do
  if (rank==0) print *, 'PASS: MPI EB face halos, exact serial parity and rollback'
  call MPI_Finalize(status)
contains
  subroutine temperature_seed_checks()
    type(eb_reactive_tree_3d) :: sample
    type(amr_tree_state_3d) :: gas,temp,out_u,out_t,kept_u,kept_t
    real(dp) :: su(nv,0:5,0:5,0:5),stemp(0:5,0:5,0:5),kept_su(nv,0:5,0:5,0:5),kept_st(0:5,0:5,0:5)
    real(dp) :: expected
    integer :: p
    sample=tree
    do p=1,5
      sample%nodes(p)%temperature=sample%nodes(p)%temperature+10._dp
    enddo
    expected=sample%nodes(1)%temperature(1,1,4)
    call fill_reactive_eb_tree_halo_3d(species,sample,2,su,stemp,ok)
    call require(ok,'default seed fixture validates')
    call require(abs(stemp(2,2,5)-expected)>1._dp,'default halo recovers physical temperature')
    call fill_reactive_eb_tree_halo_3d(species,sample,2,su,stemp,ok,preserve_temperature_seed=.true.)
    call require(ok,'serial preserved-seed halo validates')
    call require(stemp(2,2,5)==expected,'serial coarse-fine halo preserves temperature seed')
    kept_su=su;kept_st=stemp
    gas=local;temp=lt
    do p=1,5
      if(owners(p)/=rank) cycle
      temp%fields(p)%values=temp%fields(p)%values+10._dp
    enddo
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      preserve_temperature_seed=.true.)
    call require(ok,'MPI preserved-seed halo validates')
    if(owners(2)==rank) then
      call require(out_t%fields(2)%values(1,2,2,5)==expected,'MPI coarse-fine halo preserves temperature seed')
      call require(all(out_u%fields(2)%values==su),'preserved-seed exact serial MPI gas')
      call require(all(out_t%fields(2)%values(1,:,:,:)==stemp),'preserved-seed exact serial MPI temperature')
    endif
    kept_u=out_u;kept_t=out_t
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      preserve_temperature_seed=rank==0)
    call require(.not.ok,'rank-disagreeing temperature-seed policy rejected')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(out_u%fields(p)%values==kept_u%fields(p)%values).and. &
        all(out_t%fields(p)%values==kept_t%fields(p)%values),'seed-policy refusal rolls back halos')
    enddo
    sample%nodes(1)%state(1,1,1,4)=-1._dp
    call fill_reactive_eb_tree_halo_3d(species,sample,2,su,stemp,ok,preserve_temperature_seed=.true.)
    call require(.not.ok,'preserving seed still rejects invalid serial state')
    call require(all(su==kept_su).and.all(stemp==kept_st),'invalid serial seed-state rollback')
    if(owners(1)==rank) gas%fields(1)%values(1,1,1,4)=-1._dp
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      preserve_temperature_seed=.true.)
    call require(.not.ok,'preserving seed still rejects invalid MPI state')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(out_u%fields(p)%values==kept_u%fields(p)%values).and. &
        all(out_t%fields(p)%values==kept_t%fields(p)%values),'invalid MPI seed-state rollback')
    enddo
  end subroutine temperature_seed_checks
  subroutine separate_momentum_checks()
    type(eb_reactive_tree_3d) :: sample
    type(eb_geometry_3d) :: extended(5)
    type(amr_tree_state_3d) :: gas,temp,held,all_held,metrics,out_u,out_t,out_h,out_v
    type(amr_tree_state_3d) :: legacy_u,legacy_t,legacy_h,legacy_v,keep_u,keep_t,keep_h,keep_v
    real(dp) :: su(nv,-1:6,-1:6,-1:6),stemp(-1:6,-1:6,-1:6),sh(-1:6,-1:6,-1:6),sv(-1:6,-1:6,-1:6)
    real(dp) :: spacing(3),lo(3),hi(3),difference,total_difference,density,scale
    integer :: p,i,j,k
    sample=tree
    allocate(gas%fields(5),temp%fields(5),held%fields(5),all_held%fields(5),metrics%fields(5))
    do p=1,5
      do k=1,4
      do j=1,4
      do i=1,4
        if(sample%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        sample%nodes(p)%state(:,i,j,k)=(1+.01_dp*i)*sample%nodes(p)%state(:,i,j,k)
        density=sample%nodes(p)%state(1,i,j,k)
        sample%nodes(p)%state(5,i,j,k)=sample%nodes(p)%state(5,i,j,k)- &
          .5_dp*sample%nodes(p)%state(3,i,j,k)**2/density
        sample%nodes(p)%state(3,i,j,k)=0
        if(mod(i,2)==0) sample%nodes(p)%state(3,i,j,k)=-real(i,dp)*1.e-26_dp
      enddo
      enddo
      enddo
      allocate(all_held%fields(p)%values(1,4,4,4))
      all_held%fields(p)%values=0
      do k=1,4
      do j=1,4
      do i=1,4
        if(sample%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        all_held%fields(p)%values(1,i,j,k)=sample%nodes(p)%state(5,i,j,k)- &
          .5_dp*sum(sample%nodes(p)%state(2:4,i,j,k)**2)/sample%nodes(p)%state(1,i,j,k)
      enddo
      enddo
      enddo
      spacing=[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
      lo=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]-2*spacing
      hi=[geometry(p)%x_upper,geometry(p)%y_upper,geometry(p)%z_upper]+2*spacing
      call build_axis_plane_eb_geometry_3d(8,8,8,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3), &
        'z',.31_dp,extended(p),ok)
      call require(ok,'candidate halo extended geometry')
      if(owners(p)/=rank) cycle
      gas%fields(p)%values=sample%nodes(p)%state
      allocate(temp%fields(p)%values(1,4,4,4))
      temp%fields(p)%values(1,:,:,:)=sample%nodes(p)%temperature
      held%fields(p)=all_held%fields(p)
      allocate(metrics%fields(p)%values(1,-1:6,-1:6,-1:6))
      metrics%fields(p)%values(1,:,:,:)=extended(p)%volume_fraction
      metrics%fields(p)%values(1,1:4,1:4,1:4)=geometry(p)%volume_fraction
    enddo
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=out_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=out_h,coarse_geometry=extended,separate_momentum=.true.)
    call require(ok,'candidate MPI prepared halo')
    do p=1,5
      call require(allocated(out_u%fields(p)%values).eqv.(owners(p)==rank),'candidate owner-local halo')
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,sample,p,su,stemp,ok,ghost_volume_fraction=sv,halo_width=2, &
        target_volume_fraction=metrics%fields(p)%values(1,:,:,:),internal_energy=all_held, &
        ghost_internal_energy=sh,coarse_geometry=extended,separate_momentum=.true.)
      call require(ok,'candidate serial prepared halo')
      scale=max(1._dp,maxval(abs(su)))
      call require(maxval(abs(out_u%fields(p)%values-su))<=2.e-12_dp*scale,'candidate serial MPI gas parity')
      call require(maxval(abs(out_t%fields(p)%values(1,:,:,:)-stemp))<=2.e-9_dp,'candidate temperature parity')
      call require(maxval(abs(out_h%fields(p)%values(1,:,:,:)-sh))<=2.e-12_dp*scale,'candidate held energy parity')
      call require(all(out_v%fields(p)%values(1,:,:,:)==sv),'candidate exact volume fractions')
    enddo
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,legacy_u,legacy_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=legacy_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=legacy_h,coarse_geometry=extended)
    call require(ok,'legacy prepared halo remains available')
    difference=0
    do p=1,5
      if(owners(p)/=rank) cycle
      difference=max(difference,maxval(abs(legacy_u%fields(p)%values-out_u%fields(p)%values)))
    enddo
    call MPI_Allreduce(difference,total_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.total_difference>1.e-10_dp,'candidate changes noisy-momentum scalar halo')
    keep_u=out_u;keep_t=out_t;keep_h=out_h;keep_v=out_v
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=out_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=out_h,coarse_geometry=extended,separate_momentum=rank==0)
    call require(.not.ok,'mixed momentum policy collectively refused')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(out_u%fields(p)%values==keep_u%fields(p)%values).and. &
        all(out_t%fields(p)%values==keep_t%fields(p)%values).and. &
        all(out_h%fields(p)%values==keep_h%fields(p)%values).and. &
        all(out_v%fields(p)%values==keep_v%fields(p)%values),'mixed policy leaves every halo unchanged')
    enddo
    su=123;stemp=456
    call fill_reactive_eb_tree_halo_3d(species,sample,2,su,stemp,ok,separate_momentum=.true.)
    call require(.not.ok.and.all(su==123).and.all(stemp==456),'candidate requires prepared coarse stencil')
  end subroutine separate_momentum_checks

  subroutine paired_pressure_checks()
    type(eb_reactive_tree_3d) :: sample
    type(eb_geometry_3d) :: extended(5)
    type(amr_tree_state_3d) :: gas,temp,held,all_held,metrics,out_u,out_t,out_h,out_v
    type(amr_tree_state_3d) :: legacy_u,legacy_t,legacy_h,legacy_v,keep_u,keep_t,keep_h,keep_v
    real(dp) :: su(nv,-1:6,-1:6,-1:6),stemp(-1:6,-1:6,-1:6),sh(-1:6,-1:6,-1:6),sv(-1:6,-1:6,-1:6)
    real(dp) :: spacing(3),lo(3),hi(3),difference,total_difference,density,scale
    real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
    real(dp) :: flux(nv),face_density,face_velocity,face_pressure
    integer :: p,i,j,k
    sample=tree
    allocate(gas%fields(5),temp%fields(5),held%fields(5),all_held%fields(5),metrics%fields(5))
    do p=1,5
      do k=1,4
      do j=1,4
      do i=1,4
        if(sample%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        call reactive_conserved_to_primitive(species,sample%nodes(p)%state(:,i,j,k), &
          sample%nodes(p)%temperature(i,j,k),primitive,recovered,speed,ok)
        call require(ok,'pressure fixture primitive')
        primitive(6:)=0;primitive(size(primitive))=1 ! N2, 390--620 K within the pinned EOS range.
        primitive(1)=.25_dp+.02_dp*(p+i+j);primitive(2:4)=[10._dp,-2._dp,3._dp];primitive(5)=60000
        call reactive_primitive_to_conserved(species,primitive,sample%nodes(p)%state(:,i,j,k), &
          sample%nodes(p)%temperature(i,j,k),speed,ok)
        call require(ok,'isobaric variable-temperature fixture')
      enddo
      enddo
      enddo
      allocate(all_held%fields(p)%values(1,4,4,4))
      all_held%fields(p)%values=0
      do k=1,4
      do j=1,4
      do i=1,4
        if(sample%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        all_held%fields(p)%values(1,i,j,k)=sample%nodes(p)%state(5,i,j,k)- &
          .5_dp*sum(sample%nodes(p)%state(2:4,i,j,k)**2)/sample%nodes(p)%state(1,i,j,k)
      enddo
      enddo
      enddo
      spacing=[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
      lo=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]-2*spacing
      hi=[geometry(p)%x_upper,geometry(p)%y_upper,geometry(p)%z_upper]+2*spacing
      call build_axis_plane_eb_geometry_3d(8,8,8,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3), &
        'z',.31_dp,extended(p),ok)
      call require(ok,'candidate halo extended geometry')
      if(owners(p)/=rank) cycle
      gas%fields(p)%values=sample%nodes(p)%state
      allocate(temp%fields(p)%values(1,4,4,4))
      temp%fields(p)%values(1,:,:,:)=sample%nodes(p)%temperature
      held%fields(p)=all_held%fields(p)
      allocate(metrics%fields(p)%values(1,-1:6,-1:6,-1:6))
      metrics%fields(p)%values(1,:,:,:)=extended(p)%volume_fraction
      metrics%fields(p)%values(1,1:4,1:4,1:4)=geometry(p)%volume_fraction
    enddo
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=out_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=out_h,coarse_geometry=extended,separate_momentum=.true.,consistent_velocity=.true., &
      paired_pressure=.true.)
    call require(ok,'candidate MPI prepared halo')
    do p=1,5
      call require(allocated(out_u%fields(p)%values).eqv.(owners(p)==rank),'candidate owner-local halo')
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,sample,p,su,stemp,ok,ghost_volume_fraction=sv,halo_width=2, &
        target_volume_fraction=metrics%fields(p)%values(1,:,:,:),internal_energy=all_held, &
        ghost_internal_energy=sh,coarse_geometry=extended,separate_momentum=.true.,consistent_velocity=.true., &
      paired_pressure=.true.)
      call require(ok,'candidate serial prepared halo')
      scale=max(1._dp,maxval(abs(su)))
      call require(maxval(abs(out_u%fields(p)%values-su))<=2.e-12_dp*scale,'candidate serial MPI gas parity')
      call require(maxval(abs(out_t%fields(p)%values(1,:,:,:)-stemp))<=2.e-9_dp,'candidate temperature parity')
      call require(maxval(abs(out_h%fields(p)%values(1,:,:,:)-sh))<=2.e-12_dp*scale,'candidate held energy parity')
      call require(all(out_v%fields(p)%values(1,:,:,:)==sv),'candidate exact volume fractions')
      do k=-1,6
      do j=-1,6
      do i=-1,6
        if(sv(i,j,k)<=0) cycle
        call reactive_conserved_to_primitive(species,su(:,i,j,k),stemp(i,j,k),primitive,recovered,speed,ok, &
          internal_energy_density=sh(i,j,k))
        call require(ok.and.abs(primitive(5)-60000._dp)<1.e-7_dp,'paired pressure isobaric halo and covered state')
        if(i<6) then
          if(sv(i+1,j,k)>0) then
            call reactive_pelec_flux_x(species,su(:,i,j,k),su(:,i+1,j,k),stemp(i,j,k),stemp(i+1,j,k), &
              flux,ok,face_density,face_velocity,face_pressure)
            call require(ok.and.abs(face_pressure-60000._dp)<1.e-7_dp,'paired pressure Riemann interface pressure')
            call require(abs(face_velocity-10._dp)<1.e-9_dp,'paired pressure Riemann interface velocity')
            call require(abs(flux(2)-10._dp*flux(1)-60000._dp)<1.e-7_dp, &
              'paired pressure separates advected momentum and pressure flux')
          endif
        endif
      enddo
      enddo
      enddo
      call require(all(gas%fields(p)%values==sample%nodes(p)%state),'pressure representation preserves input gas')
      call require(all(held%fields(p)%values==all_held%fields(p)%values),'pressure representation preserves input held')
    enddo
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,legacy_u,legacy_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=legacy_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=legacy_h,coarse_geometry=extended)
    call require(ok,'legacy prepared halo remains available')
    difference=0
    do p=1,5
      if(owners(p)/=rank) cycle
      difference=max(difference,maxval(abs(legacy_u%fields(p)%values-out_u%fields(p)%values)))
    enddo
    call MPI_Allreduce(difference,total_difference,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.total_difference>1.e-10_dp,'pressure representation changes nonlinear EOS restriction')
    keep_u=out_u;keep_t=out_t;keep_h=out_h;keep_v=out_v
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=out_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=out_h,coarse_geometry=extended,separate_momentum=.true.,consistent_velocity=.true., &
      paired_pressure=rank==0)
    call require(.not.ok,'mixed momentum policy collectively refused')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(out_u%fields(p)%values==keep_u%fields(p)%values).and. &
        all(out_t%fields(p)%values==keep_t%fields(p)%values).and. &
        all(out_h%fields(p)%values==keep_h%fields(p)%values).and. &
        all(out_v%fields(p)%values==keep_v%fields(p)%values),'mixed policy leaves every halo unchanged')
    enddo
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=out_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=out_h,coarse_geometry=extended,separate_momentum=.true.,consistent_velocity=.true., &
      paired_pressure=.true.,connected_directions=rank==0)
    call require(.not.ok,'mixed connected-direction policy collectively refused')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(out_u%fields(p)%values==keep_u%fields(p)%values).and. &
        all(out_t%fields(p)%values==keep_t%fields(p)%values).and. &
        all(out_h%fields(p)%values==keep_h%fields(p)%values).and. &
        all(out_v%fields(p)%values==keep_v%fields(p)%values),'mixed connected policy leaves every halo unchanged')
    enddo
    if(owners(1)==rank) held%fields(1)%values(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,gas,temp,out_u,out_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=out_v,halo_width=2,target_volume_fraction=metrics,internal_energy=held, &
      ghost_internal_energy=out_h,coarse_geometry=extended,separate_momentum=.true., &
      consistent_velocity=.true.,paired_pressure=.true.)
    call require(.not.ok,'paired pressure remote held EOS refusal')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(out_u%fields(p)%values==keep_u%fields(p)%values).and. &
        all(out_t%fields(p)%values==keep_t%fields(p)%values).and. &
        all(out_h%fields(p)%values==keep_h%fields(p)%values).and. &
        all(out_v%fields(p)%values==keep_v%fields(p)%values),'pressure EOS refusal rolls back every halo')
    enddo
    su=123;stemp=456
    call fill_reactive_eb_tree_halo_3d(species,sample,2,su,stemp,ok,separate_momentum=.true., &
      consistent_velocity=.true.,paired_pressure=.true.)
    call require(.not.ok.and.all(su==123).and.all(stemp==456),'candidate requires prepared coarse stencil')
  end subroutine paired_pressure_checks

  subroutine periodic_checks(walls)
    logical, optional, intent(in) :: walls
    type(amr_tree_state_3d) :: volume
    real(dp) :: serial_vf(0:5,0:5,0:5)
    character(len=16) :: labels(6)
    integer :: p
    labels='outflow';labels(:4)='periodic'
    if(present(walls)) then
      labels(1:2)='slip_wall';labels(3:4)='no_slip_wall'
    end if
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      boundary_faces=labels)
    call require(ok,'periodic MPI face halos')
    call unchanged_inputs()
    do p=1,5
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,tree,p,reference,reference_t,ok,boundary_faces=labels)
      call require(ok.and.all(ghost%fields(p)%values==reference).and. &
        all(gh_t%fields(p)%values(1,:,:,:)==reference_t),'periodic face exact serial parity')
    end do
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      boundary_faces=labels,ghost_volume_fraction=volume)
    call require(ok,'periodic complete one-layer halo')
    do p=1,5
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,tree,p,reference,reference_t,ok, &
        boundary_faces=labels,ghost_volume_fraction=serial_vf)
      call require(ok.and.all(ghost%fields(p)%values==reference).and. &
        all(gh_t%fields(p)%values(1,:,:,:)==reference_t).and. &
        all(volume%fields(p)%values(1,:,:,:)==serial_vf),'periodic complete halo exact serial parity')
    end do
    saved_ghost=ghost;saved_gt=gh_t
    if(rank==0) labels(:2)='outflow'
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      boundary_faces=labels)
    call require(.not.ok,'rank-disagreeing periodic labels reject')
    call unchanged_halos()
    labels='periodic'
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      boundary_faces=labels)
    call require(.not.ok,'incompatible periodic root faces reject')
    call unchanged_halos()
  end subroutine
  subroutine prescribed_checks()
    real(dp) :: inlet(nv,6),inlet_t(6)
    character(len=16) :: labels(6)
    integer :: p
    labels='outflow';labels(1)='inflow';labels(3)='inflow';inlet_t=t(1,1,3)
    do p=1,6
      inlet(:,p)=real(p+1,dp)*u(:,1,1,3)
    end do
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(ok,'prescribed MPI face halos')
    do p=1,5
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,tree,p,reference,reference_t,ok, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
      call require(ok.and.all(ghost%fields(p)%values==reference).and. &
        all(gh_t%fields(p)%values(1,:,:,:)==reference_t),'prescribed face exact serial parity')
    end do
    saved_ghost=ghost;saved_gt=gh_t
    inlet(1,1)=-1
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(.not.ok,'invalid prescribed EOS collective rejection')
    call unchanged_halos()
  end subroutine
  subroutine wide_checks(prescribed,periodic,walls)
    logical, optional, intent(in) :: prescribed,periodic,walls
    character(len=16), allocatable :: labels(:)
    real(dp), allocatable :: inlet(:,:),inlet_t(:)
    type(amr_tree_state_3d) :: metrics,wide_u,wide_t,wide_vf,old_u,old_t,old_vf,end_u
    type(eb_reactive_tree_3d) :: end_tree
    type(eb_geometry_3d) :: extended
    real(dp) :: serial_u(nv,-5:10,-5:10,-5:10),serial_t(-5:10,-5:10,-5:10),serial_vf(-5:10,-5:10,-5:10)
    real(dp) :: spacing(3),lo(3),hi(3)
    integer :: p,width,depth
    real(dp) :: fractions(2)
    if(present(prescribed)) then
      allocate(labels(6),inlet(nv,6),inlet_t(6))
      labels='outflow';labels(1)='inflow';labels(3)='inflow'
      do p=1,6
        inlet(:,p)=real(p+1,dp)*u(:,1,1,3)
      end do
      inlet_t=t(1,1,3)
    end if
    if(present(periodic)) then
      allocate(labels(6))
      labels='outflow';labels(:4)='periodic'
    end if
    if(present(walls)) then
      allocate(labels(6))
      labels='outflow';labels(1:2)='slip_wall';labels(3:4)='no_slip_wall'
    end if
    allocate(metrics%fields(5))
    do p=1,5
      if(owners(p)/=rank) cycle
      spacing=[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
      lo=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]-6*spacing
      hi=[geometry(p)%x_upper,geometry(p)%y_upper,geometry(p)%z_upper]+6*spacing
      call build_axis_plane_eb_geometry_3d(16,16,16,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),'z',.31_dp,extended,ok)
      call require(ok,'wide target geometry')
      allocate(metrics%fields(p)%values(1,-5:10,-5:10,-5:10))
      metrics%fields(p)%values(1,:,:,:)=extended%volume_fraction
      ! Preserve the exact metrics produced by hierarchy coarsening in the center.
      metrics%fields(p)%values(1,1:4,1:4,1:4)=geometry(p)%volume_fraction
    end do
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,wide_u,wide_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=wide_vf,halo_width=6,target_volume_fraction=metrics, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(ok,'six-layer MPI EB halos')
    do p=1,5
      call require(allocated(wide_u%fields(p)%values).eqv.(owners(p)==rank),'wide halos remain owner local')
      if(owners(p)/=rank) cycle
      call fill_reactive_eb_tree_halo_3d(species,tree,p,serial_u,serial_t,ok,ghost_volume_fraction=serial_vf, &
        halo_width=6,target_volume_fraction=metrics%fields(p)%values(1,:,:,:), &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
      call require(ok,'six-layer serial EB halos')
      call require(all(lbound(wide_u%fields(p)%values)==[1,-5,-5,-5]),'wide halo spatial bounds')
      call require(all(wide_u%fields(p)%values==serial_u),'wide peer/ancestor/outflow U parity')
      call require(all(wide_t%fields(p)%values(1,:,:,:)==serial_t),'wide peer/ancestor/outflow T parity')
      call require(all(wide_vf%fields(p)%values(1,:,:,:)==serial_vf),'wide fine-metric parity')
      call require(all(serial_vf==metrics%fields(p)%values(1,:,:,:)),'target geometry controls every halo layer')
    end do
    if(rank==owners(1).and..not.present(periodic).and..not.present(walls)) then
      call require(all(wide_u%fields(1)%values(:,-5,1,3)==wide_u%fields(1)%values(:,0,1,3)), &
        'all six physical outflow layers use the same boundary donor')
    end if
    end_u=local;end_tree=tree;fractions=[.25_dp,.75_dp]
    do p=1,5
      end_tree%nodes(p)%state=2*tree%nodes(p)%state
      if(owners(p)/=rank) cycle
      end_u%fields(p)%values=2*local%fields(p)%values
    end do
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,wide_u,wide_t,MPI_COMM_WORLD,ok, &
      ancestor_start=local,ancestor_end=end_u,temperature_start=lt,temperature_end=lt,ancestor_fraction=fractions, &
      ghost_volume_fraction=wide_vf,halo_width=6,target_volume_fraction=metrics, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(ok,'wide temporal MPI halos')
    do p=1,5
      if(owners(p)/=rank) cycle
      depth=topology%nodes(p)%level
      call fill_reactive_eb_tree_halo_3d(species,tree,p,serial_u,serial_t,ok,tree,end_tree,fractions(:depth), &
        ghost_volume_fraction=serial_vf,halo_width=6,target_volume_fraction=metrics%fields(p)%values(1,:,:,:), &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
      call require(ok,'wide temporal serial halos')
      call require(all(wide_u%fields(p)%values==serial_u).and.all(wide_t%fields(p)%values(1,:,:,:)==serial_t), &
        'wide ancestor time interpolation matches serial exactly')
    end do
    old_u=wide_u;old_t=wide_t;old_vf=wide_vf
    if(allocated(inlet)) then
      if(rank==0) inlet(1,1)=inlet(1,1)+.01_dp
      call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,wide_u,wide_t,MPI_COMM_WORLD,ok, &
        ghost_volume_fraction=wide_vf,halo_width=6,target_volume_fraction=metrics, &
        boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
      call require(.not.ok,'rank-disagreeing prescribed state rejection')
      if(rank==0) inlet(1,1)=2*u(1,1,1,3)
      do p=1,5
        if(owners(p)/=rank) cycle
        call require(all(wide_u%fields(p)%values==old_u%fields(p)%values).and. &
          all(wide_t%fields(p)%values==old_t%fields(p)%values).and. &
          all(wide_vf%fields(p)%values==old_vf%fields(p)%values),'inlet disagreement rollback')
      end do
    end if
    width=6
    if(rank==0) width=5
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,wide_u,wide_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=wide_vf,halo_width=width,target_volume_fraction=metrics, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(.not.ok,'collective wide halo mismatch rejection')
    do p=1,5
      if(owners(p)/=rank) cycle
      call require(all(wide_u%fields(p)%values==old_u%fields(p)%values).and. &
        all(wide_t%fields(p)%values==old_t%fields(p)%values).and. &
        all(wide_vf%fields(p)%values==old_vf%fields(p)%values),'wide collective rollback')
    end do
    if(rank==owners(1)) metrics%fields(1)%values(1,1,1,1)=.5_dp
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,wide_u,wide_t,MPI_COMM_WORLD,ok, &
      ghost_volume_fraction=wide_vf,halo_width=6,target_volume_fraction=metrics, &
      boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
    call require(.not.ok,'wide central geometry mismatch rejection')
  end subroutine wide_checks

  subroutine temporal_checks()
    type(eb_reactive_tree_3d) :: start_tree,end_tree
    type(amr_tree_state_3d) :: start_u,end_u,start_t,end_t
    real(dp) :: fractions(2)
    integer :: p,trial,depth
    start_tree=tree;end_tree=tree
    start_u=local;end_u=local;start_t=lt;end_t=lt
    do p=1,5
      start_tree%nodes(p)%state=real(p+1,dp)*tree%nodes(p)%state
      end_tree%nodes(p)%state=real(p+3,dp)*tree%nodes(p)%state
      if (owners(p)/=rank) cycle
      start_u%fields(p)%values=start_tree%nodes(p)%state
      end_u%fields(p)%values=end_tree%nodes(p)%state
    end do
    do trial=1,3
      select case(trial)
      case(1);fractions=[0._dp,1._dp]
      case(2);fractions=[1._dp,0._dp]
      case(3);fractions=[.25_dp,.75_dp]
      end select
      call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
        start_u,end_u,start_t,end_t,fractions)
      call require(ok,'MPI temporal halos')
      call unchanged_inputs()
      do p=1,5
        if (owners(p)/=rank) cycle
        depth=topology%nodes(p)%level
        call fill_reactive_eb_tree_halo_3d(species,tree,p,reference,reference_t,ok, &
          start_tree,end_tree,fractions(:depth))
        call require(ok,'serial temporal halos')
        call require(all(ghost%fields(p)%values==reference),'exact temporal U')
        call require(all(gh_t%fields(p)%values(1,:,:,:)==reference_t),'exact temporal T')
        call require(all(start_u%fields(p)%values==start_tree%nodes(p)%state),'start input unchanged')
        call require(all(end_u%fields(p)%values==end_tree%nodes(p)%state),'end input unchanged')
      end do
    end do
    saved_ghost=ghost;saved_gt=gh_t
    if (rank==0) fractions(1)=.5_dp
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      start_u,end_u,start_t,end_t,fractions)
    call require(.not.ok,'rank-disagreeing fractions rejected')
    call unchanged_halos()
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      ancestor_start=start_u)
    call require(.not.ok,'partial temporal context rejected')
    call unchanged_halos()
    fractions=[.25_dp,.75_dp]
    if (rank==owners(1)) end_u%fields(1)%values(1,:,:,:)=-100._dp
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      start_u,end_u,start_t,end_t,fractions)
    call require(.not.ok,'invalid interpolated ancestor EOS rejected')
    call unchanged_halos()
    if (rank==owners(1)) end_u%fields(1)%values=ieee_value(0._dp,ieee_quiet_nan)
    fractions=0._dp
    call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
      start_u,end_u,start_t,end_t,fractions)
    call require(ok,'unused end state ignored at start time')
    saved_ghost=ghost;saved_gt=gh_t
    if (rank==0) then
      call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok)
    else
      call fill_mpi_eb_tree_halo_3d(species,topology,owners,geometry,local,lt,ghost,gh_t,MPI_COMM_WORLD,ok, &
        start_u,end_u,start_t,end_t,fractions)
    end if
    call require(.not.ok,'rank-disagreeing context presence rejected')
    call unchanged_halos()
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (condition) return
    print *, 'FAIL rank ',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
  subroutine unchanged_inputs()
    integer :: n
    do n=1,5
      if (owners(n)/=rank) cycle
      call require(all(local%fields(n)%values==saved%fields(n)%values),'input U unchanged')
      call require(all(lt%fields(n)%values==st%fields(n)%values),'input T unchanged')
    end do
  end subroutine
  subroutine unchanged_halos()
    integer :: n
    do n=1,5
      if (owners(n)/=rank) cycle
      call require(all(ghost%fields(n)%values==saved_ghost%fields(n)%values),'halo U rollback')
      call require(all(gh_t%fields(n)%values==saved_gt%fields(n)%values),'halo T rollback')
    end do
  end subroutine
end program
