program test_mpi_eb_tree_hydro_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_hydro_3d, &
    build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_hydro_3d_mod, only: advance_mpi_eb_tree_hydro_3d,compute_sdc_mpi_eb_tree_hydro_3d
  use eb_reactive_patch_tree_3d_mod, only: compute_sdc_eb_tree_hydro_3d
  use eb_reactive_patch_tree_3d_mod, only: fill_sdc_eb_tree_source_halo_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_sdc_mpi_eb_tree_source_halo_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,extended_finest
  type(eb_geometry_3d), allocatable :: extended(:),physical(:)
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),inlet(:,:),inlet_t(:)
  character(len=16), allocatable :: labels(:)
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),primitive(:)
  real(dp) :: residual,coefficient,et,sound
  character(len=32) :: test_mode,solver
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,rank,nranks,status,i,j,k
  logical :: ok
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
  do axis=1,3
  do count=4,5
  if(index(test_mode,'ctu')==1.and.(axis/=3.or.count/=5)) cycle
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
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  if(index(test_mode,'ctu')==1) then
    allocate(physical(count))
    do node=1,count
      physical(node)=tree%nodes(node)%geometry
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
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
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
    type(amr_tree_state_3d) :: local,lt,predictor,local_predictor,rates,reference,prior
    type(eb_integrated_flux_3d), allocatable :: flux(:),expected(:),old_flux(:)
    type(amr_tree_state_3d) :: end_local,end_source,local_end_source,level_rates,level_reference
    type(eb_integrated_flux_3d), allocatable :: level_flux(:),level_expected(:)
    real(dp) :: dt,nan
    integer :: owners(count),p,mode,failure,c
    logical :: good
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
      difmag=coefficient,boundary_faces=labels,inflow_state=inlet,inflow_temperature=inlet_t)
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
      call check_temporal_sources(owners,predictor,local_predictor)
      call compute_sdc_mpi_eb_tree_hydro_3d(species,topology,owners,physical,local,lt,local_predictor, &
        extended,1.e-8_dp,rates,flux,MPI_COMM_WORLD,good,difmag=coefficient,boundary_faces=labels, &
        inflow_state=inlet,inflow_temperature=inlet_t)
      call require(good,'MPI SDC hydro with covered NaNs')
      do p=1,count
        call require(allocated(rates%fields(p)%values).eqv.(owners(p)==rank),'owner-local SDC rates')
        call require(allocated(flux(p)%x).eqv.(owners(p)==rank),'owner-local SDC flux')
        if(owners(p)/=rank) cycle
        call require(all(rates%fields(p)%values==reference%fields(p)%values),'exact serial SDC hydro rate')
        call require(all(flux(p)%x==expected(p)%x).and.all(flux(p)%y==expected(p)%y).and. &
          all(flux(p)%z==expected(p)%z),'exact serial SDC hydro flux')
        call require(all(local%fields(p)%values==stage_tree%nodes(p)%state),'SDC preserves input U')
        call require(all(lt%fields(p)%values(1,:,:,:)==stage_tree%nodes(p)%temperature),'SDC preserves input T')
      end do
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
