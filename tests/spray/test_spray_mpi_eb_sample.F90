program test_spray_mpi_eb_sample
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use spray_mpi_amr_step_mod, only: kick_mpi_amr_mixture_halo, &
    advance_mpi_eb_ghost_stage,advance_mpi_eb_population_stage, &
    advance_mpi_eb_population_sdc_stage,spray_sdc_population_stage,advance_mpi_eb_recursive_population_stage
  use spray_mpi_amr_film_step_mod, only: kick_mpi_eb_film_halo,advance_mpi_eb_film_sdc_stage,spray_sdc_film_stage, &
    advance_mpi_eb_recursive_film_stage,finish_mpi_eb_film_level
  use spray_mpi_amr_film_mod, only: spray_level_film_population,initialize_mpi_eb_film_levels
  use spray_wall_impact_mod, only: spray_wall_film,spray_wall_options,advance_mixture_wall_film_frozen
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use spray_mpi_amr_tree_mod, only: spray_level_population,initialize_mpi_amr_sdc_populations
  use spray_parcel_mod, only: spray_parcel,spray_liquid,spray_options,initialize_spray_liquid,liquid_energy
  use spray_mixture_mod, only: advance_mixture_frozen_gas
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use spray_mpi_amr_field_mod, only: sample_mpi_amr_tree_spray_field,deposit_mpi_amr_tree_spray_field, &
    sample_mpi_eb_tree_particle_gas
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use spray_eb_tree_stencil_mod, only: eb_tree_spray_stencil
  implicit none
  type(amr_tree_state_3d) :: field,change,initial_field
  real(dp), allocatable :: positions(:,:),samples(:,:),kept(:,:),contribution(:,:)
  integer :: rank,nranks,ierr,owners(2),dims(3)
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: finest,g(2)
  real(dp) :: phi(0:8,0:8,0:8),position(3),zero(3)=0,spacing(3)=1,expected,total
  real(dp), allocatable :: w(:),v(:),saved_w(:),saved_v(:)
  integer, allocatable :: cells(:,:),saved(:,:)
  integer :: j,q,node
  logical :: ok,periodic(3)=.false.
  character(len=16) :: mode
  call get_command_argument(1,mode)
  if(mode=='periodic') periodic=[.true.,.false.,.true.]
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'rank count')
  owners=[0,nranks-1]
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
  call require(ok,'topology')
  do j=0,8
    phi(:,j,:)=real(j,dp)*.5_dp-1.25_dp
  end do
  call build_level_set_eb_geometry_3d(phi,zero,[4._dp,4._dp,4._dp],finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,g(1),ok)
  call require(ok,'root geometry')
  call extract_eb_geometry_3d(finest,[1,1,1],[4,8,8],g(2),ok)
  call require(ok,'child geometry')
  position=[1.9_dp,1.3_dp,1.5_dp]
  if(any(periodic)) position=[3.9_dp,1.3_dp,3.9_dp]
  call run()
  call require(ok.and.all(w>0).and.abs(sum(w)-1)<1e-14_dp,'normalized cut-interface weights')
  call require(any(cells(1,:)==1).and.any(cells(1,:)==2),'coarse and fine support')
  allocate(field%fields(2),positions(3,merge(1,0,rank==0)),samples(1,merge(1,0,rank==0)))
  do node=1,2
    if(owners(node)/=rank) cycle
    dims=tree%nodes(node)%shape
    allocate(field%fields(node)%values(1,dims(1),dims(2),dims(3)))
    field%fields(node)%values=real(10*node,dp)
  end do
  expected=sum(w*real(10*cells(1,:),dp))
  call require(expected>10.and.expected<20,'two-owner reference')
  if(rank==0) positions(:,1)=position
  samples=99;kept=samples
  call sample(0)
  call require(.not.ok.and.all(samples==kept),'request capacity rollback')
  if(rank==nranks-1) phi(0,0,0)=phi(0,0,0)+.01_dp
  call sample(100)
  call require(.not.ok.and.all(samples==kept),'rank geometry mismatch rollback')
  phi(0,0,0)=-1.25_dp
  if(rank==0) positions(2,1)=1._dp
  call sample(100)
  call require(.not.ok.and.all(samples==kept),'solid query rollback')
  if(rank==0) positions(:,1)=position
  call sample(100)
  call require(ok,'owner-routed EB sampling')
  if(rank==0) call require(abs(samples(1,1)-expected)<1e-13_dp,'weighted remote values')
  if(rank/=0) call require(size(samples,2)==0,'empty query ranks')
  call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,field,1,positions,zero,spacing, &
    periodic,samples,ok,100,g,phi,maximum_level=0)
  call require(ok.and.all(abs(samples-10._dp)<1e-13_dp),'covered coarse EB sampling')
  initial_field=field;allocate(contribution(1,size(positions,2)));contribution=2
  call deposit_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,field,1,positions,contribution, &
    zero,spacing,periodic,ok,eb_geometry=g,nodal_phi=phi,gas_change=change,maximum_level=0)
  call require(ok,'coarse EB source deposition')
  total=0
  if(rank==owners(1)) total=sum(change%fields(1)%values(1,:,:,:)*g(1)%volume_fraction)*g(1)%dx*g(1)%dy*g(1)%dz
  call MPI_Allreduce(MPI_IN_PLACE,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(abs(total-2._dp)<1e-13_dp,'coarse EB fluid-volume source inventory')
  if(rank==owners(2)) call require(all(field%fields(2)%values==initial_field%fields(2)%values).and. &
    all(change%fields(2)%values==0),'fine EB field unchanged by coarse source')
  deallocate(positions);allocate(positions(3,0))
  call sample(100)
  call require(ok.and.size(samples,2)==0,'globally empty queries')
  call check_temporal_particle_gas()
  if(any(periodic)) then
    phi(0,8,8)=phi(0,8,8)+.1_dp
    call sample(100)
    call require(.not.ok.and.size(samples,2)==0,'empty periodic query still validates nodal endpoints')
  end if
  if(rank==0) print *, 'PASS: owner-local MPI EB composite particle sampling'
  call MPI_Finalize(ierr)
contains
  subroutine check_temporal_particle_gas()
    type(nasa7_species), allocatable :: species(:)
    type(amr_tree_state_3d) :: gas,temp,ending,metrics,centroids,source,remainder,old_source
    type(gas_transport_species), allocatable :: transport(:)
    type(spray_liquid) :: liquid(1)
    type(spray_options) :: options,child_options
    type(spray_wall_options) :: film_options
    type(spray_level_film_population) :: film_groups(3),old_film_groups(3)
    type(spray_sdc_film_stage) :: film_work,depletion_work,film_clocks(2)
    type(spray_level_film_population), allocatable :: film_physical(:),film_virtuals(:),film_ghosts(:)
    type(spray_level_film_population) :: root_finished_film
    type(amr_tree_state_3d) :: expected_film_source,expected_film_remainder,expected_physical_film
    real(dp) :: old_film_wall(5),expected_film_wall(5)
    type(spray_wall_film), allocatable :: films(:),reference_films(:),kept_films(:)
    type(amr_tree_state_3d) :: film_halo,film_halo_t,film_volume
    real(dp), allocatable :: film_mass(:,:),reference_film_mass(:,:),probes(:,:)
    integer, allocatable :: film_nodes(:)
    real(dp) :: expected_wall,film_wall,film_energy_change
    type(spray_parcel), allocatable :: ghosts(:)
    type(spray_parcel) :: reference
    type(spray_level_population) :: group(3),before_group(3),expected_group(3)
    real(dp) :: accel(3),dv(3),expected_force(5),old_force(5),parcel_weight,birth_position(3)
    integer(int64) :: child_id,impact_id
    type(spray_wall_film), allocatable :: impact_films(:)
    real(dp), allocatable :: impact_mass(:,:)
    real(dp) :: impact_total(3),initial_total(3)
    logical :: rejected
    integer :: trial
    type(spray_sdc_population_stage) :: work,hierarchy_work(2),child_work
    type(spray_level_population), allocatable :: physical(:),virtuals(:),incoming(:)
    type(spray_level_population) :: root_finished
    type(amr_tree_state_3d) :: expected_source,expected_remainder
    type(amr_tree_state_3d) :: group_source,group_remainder,old_group_source,old_group_remainder
    type(amr_tree_state_3d) :: physical_source,old_physical_source
    real(dp), allocatable :: mass(:,:),old_mass(:,:),delta(:),integral(:),ledger(:)
    real(dp) :: wall(5),reference_mass(1)
    real(dp), allocatable :: primitive(:),cell(:),y(:),query(:,:),values(:,:),site_t(:),saved_values(:,:),saved_t(:)
    integer, allocatable :: query_nodes(:)
    real(dp) :: h(3),cell_lower(3),vf
    real(dp) :: gas_constant,t,sound,fraction(1),fractions(3)=[0._dp,.25_dp,1._dp]
    integer :: nv,p,a,b,c,n(3),count,stage,population,substep
    character(len=16) :: faces(6)
    call load_h2o2_full_thermo(species,ok)
    call require(ok,'particle halo thermo')
    nv=reactive_nvar(size(species));allocate(primitive(nv),cell(nv),y(size(species)))
    y=0
    do p=1,size(species)
      if(species(p)%name=='N2') y(p)=1
    end do
    gas_constant=mixture_specific_gas_constant(species,y,ok)
    call require(ok,'particle halo mixture')
    primitive=0;primitive(1)=1;primitive(5)=gas_constant*900;primitive(6:)=y
    call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
    call require(ok,'particle halo initial state')
    allocate(gas%fields(2),temp%fields(2),metrics%fields(2),centroids%fields(2))
    do p=1,2
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      allocate(gas%fields(p)%values(nv,n(1),n(2),n(3)),temp%fields(p)%values(1,n(1),n(2),n(3)))
      gas%fields(p)%values=0;temp%fields(p)%values=900
      do c=1,n(3)
      do b=1,n(2)
      do a=1,n(1)
        if(g(p)%volume_fraction(a,b,c)>0) gas%fields(p)%values(:,a,b,c)=cell
      end do
      end do
      end do
    end do
    do p=1,2
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape;h=[g(p)%dx,g(p)%dy,g(p)%dz]
      allocate(metrics%fields(p)%values(1,0:n(1)+1,0:n(2)+1,0:n(3)+1))
      allocate(centroids%fields(p)%values(3,0:n(1)+1,0:n(2)+1,0:n(3)+1))
      do c=0,n(3)+1
      do b=0,n(2)+1
      do a=0,n(1)+1
        cell_lower=[g(p)%x_lower,g(p)%y_lower,g(p)%z_lower]+real([a,b,c]-1,dp)*h
        vf=max(0._dp,min(1._dp,(cell_lower(2)+h(2)-1.25_dp)/h(2)))
        metrics%fields(p)%values(1,a,b,c)=vf
        centroids%fields(p)%values(:,a,b,c)=cell_lower+.5_dp*h
        if(vf>0) centroids%fields(p)%values(2,a,b,c)=cell_lower(2)+h(2)-.5_dp*vf*h(2)
      end do
      end do
      end do
      metrics%fields(p)%values(1,1:n(1),1:n(2),1:n(3))=g(p)%volume_fraction
    end do
    ending=gas
    if(rank==owners(1)) ending%fields(1)%values=2*gas%fields(1)%values
    count=merge(1,0,rank==owners(2));allocate(query_nodes(count),query(3,count));query_nodes=2
    if(count>0) query(:,1)=[2.25_dp,2.25_dp,2.25_dp]
    faces='outflow'
    if(periodic(1)) faces(1:2)='periodic'
    if(periodic(3)) faces(5:6)='periodic'
    do stage=1,3
      fraction=fractions(stage)
      call sample_mpi_eb_tree_particle_gas(MPI_COMM_WORLD,species,tree,owners,g,gas,temp,1,query_nodes,query, &
        values,site_t,ok,gas,ending,temp,temp,fraction,target_volume_fraction=metrics,boundary_faces=faces, &
        nodal_phi=phi,centroids=centroids)
      call require(ok,'solver-generated temporal particle halo')
      if(count>0) then
        call require(maxval(abs(values(:,1)-(1+fraction(1))*cell)/max(abs(cell),1._dp))<1e-13_dp, &
          'parent time endpoints reach fine ghost particle')
        call require(abs(site_t(1)-900)<1e-8_dp,'particle-site EOS temperature')
      end if
    end do
    if(.not.any(periodic)) then
      call load_h2o2_full_transport(transport,ok)
      call require(ok,'ghost stage transport')
      call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquid(1),ok)
      call require(ok,'ghost stage liquid')
      allocate(ghosts(count),mass(1,count),delta(nv),integral(nv),ledger(nv))
      do stage=1,3
        fraction=fractions(stage);delta=0
        if(count>0) then
          ghosts(1)=spray_parcel(1_int64,query(:,1),[20._dp,0._dp,0._dp],1e-12_dp,350._dp,100._dp)
          mass(:,1)=ghosts(1)%mass;reference=ghosts(1);reference_mass=mass(:,1)
          call reactive_conserved_to_primitive(species,(1+fraction(1))*cell,900._dp,primitive,t,sound,ok)
          call require(ok,'time-interpolated ghost EOS reference')
          call advance_mixture_frozen_gas(species,transport,primitive,t,reference,liquid,reference_mass, &
            options,1e-8_dp,delta,ok)
          call require(ok,'time-interpolated ghost physical reference')
        end if
        do population=1,3
          group(population)%parcels=ghosts;group(population)%component_mass=mass;group(population)%nodes=query_nodes
          if(count>0) group(population)%parcels(1)%id=int(population,int64)
        end do
        call advance_mpi_eb_ghost_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
          1,1,1,100,query_nodes,ghosts,mass,1e-8_dp,0._dp,faces,source,remainder,wall,ledger,ok,phi,centroids,metrics, &
          gas,ending,temp,temp,fraction)
        call require(ok.and.size(ghosts)==count,'stage-time EB halo reaches ghost evolution')
        if(count>0) then
          call require(maxval(abs(ghosts(1)%velocity-reference%velocity))<1e-12_dp.and. &
            abs(ghosts(1)%temperature-reference%temperature)<1e-10_dp.and. &
            abs(mass(1,1)-reference_mass(1))<1e-14_dp*reference_mass(1),'ghost responds to parent endpoint interpolation')
        end if
        integral=0
        do p=1,2
          if(owners(p)/=rank) cycle
          n=tree%nodes(p)%shape
          do c=1,n(3)
          do b=1,n(2)
          do a=1,n(1)
            integral=integral+source%fields(p)%values(:,a,b,c)*g(p)%volume_fraction(a,b,c)*g(p)%dx*g(p)%dy*g(p)%dz
          end do
          end do
          end do
        end do
        ! This query lies at the first fine ghost-cell center. Its reaction
        ! stays in the halo remainder rather than entering the coarse gas.
        do p=1,2
          if(owners(p)/=rank) cycle
          call require(all(source%fields(p)%values==0),'ghost-only reaction never feeds coarse valid cells')
          n=tree%nodes(p)%shape
          do c=0,n(3)+1
          do b=0,n(2)+1
          do a=0,n(1)+1
            integral=integral+remainder%fields(p)%values(:,a,b,c)*metrics%fields(p)%values(1,a,b,c)* &
              g(p)%dx*g(p)%dy*g(p)%dz
          end do
          end do
          end do
        end do
        call MPI_Allreduce(MPI_IN_PLACE,integral,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call require(all(abs(integral-delta)<1e-12_dp*max(abs(delta),1e-30_dp)), &
          'stage-time EB valid plus retained halo reaction inventory')
        call advance_mpi_eb_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
          1,1,1,100,group,1e-8_dp,0._dp,faces,group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics, &
          gas,ending,temp,temp,fraction,physical_source=physical_source)
        call require(ok,'combined physical virtual ghost stage')
        call integrate_mpi_eb_tree_3d(tree,owners,g,physical_source,integral,MPI_COMM_WORLD,ok)
        call require(ok,'integrate physical-only source on composite leaves')
        call require(all(abs(integral-delta)<1e-12_dp*max(abs(delta),1e-30_dp)), &
          'physical-only export excludes virtual and ghost copies and retains off-level reaction')
        do population=1,3
          call require(all(group(population)%component_mass==mass),'all populations sample identical stage-time gas')
          if(count>0) call require(all(group(population)%parcels(1)%velocity==ghosts(1)%velocity).and. &
            group(population)%parcels(1)%temperature==ghosts(1)%temperature,'combined population physical update')
        end do
        do p=1,2
          if(owners(p)/=rank) cycle
          call require(all(group_source%fields(p)%values==3*source%fields(p)%values),'combined valid reaction')
          call require(all(abs(group_remainder%fields(p)%values-3*remainder%fields(p)%values)< &
            1e-12_dp*max(abs(remainder%fields(p)%values),1e-30_dp)),'combined retained halo reaction')
        end do
      end do
      ! Reject the third population after the first two have advanced privately.
      if(count>0) group(3)%nodes=0
      before_group=group;old_group_source=group_source;old_group_remainder=group_remainder
      old_physical_source=physical_source
      wall=17;ledger=19
      call advance_mpi_eb_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,group,1e-8_dp,0._dp,faces,group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics, &
        gas,ending,temp,temp,fraction,physical_source=physical_source)
      call require(.not.ok.and.all(wall==17).and.all(ledger==19),'late combined population failure')
      do population=1,3
        call require(all(group(population)%component_mass==before_group(population)%component_mass).and. &
          all(group(population)%nodes==before_group(population)%nodes),'all population payloads roll back')
        if(count>0) call require(all(group(population)%parcels(1)%velocity== &
          before_group(population)%parcels(1)%velocity).and.group(population)%parcels(1)%temperature== &
          before_group(population)%parcels(1)%temperature,'earlier population kicks roll back')
      end do
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==old_group_source%fields(p)%values).and. &
          all(group_remainder%fields(p)%values==old_group_remainder%fields(p)%values),'combined source rollback')
        call require(all(physical_source%fields(p)%values==old_physical_source%fields(p)%values), &
          'physical-only export rollback after late copy failure')
      end do
      group(3)%nodes=query_nodes
      do population=1,3
        if(count>0) group(population)%parcels(1)%velocity(1)=-20._dp
      end do
      before_group=group;fraction=.25_dp
      accel=[1000._dp,-2000._dp,3000._dp];dv=1e-8_dp*accel
      expected_force=0
      if(count>0) then
        parcel_weight=group(1)%parcels(1)%mass*group(1)%parcels(1)%multiplicity
        expected_force(2:4)=parcel_weight*dv
        expected_force(5)=parcel_weight*dot_product(group(1)%parcels(1)%velocity+.5_dp*dv,dv)
      end if
      call MPI_Allreduce(MPI_IN_PLACE,expected_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,group,work,0._dp,2e-8_dp,0,2,faces,group_source,group_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction,acceleration=accel)
      call require(ok.and.work%next_stage==1,'old SDC kick and drift cached')
      call require(all(work%force_change==expected_force),'old external force counts physical particles only')
      old_force=work%force_change
      do population=1,3
        call require(all(group(population)%component_mass==before_group(population)%component_mass), &
          'old stage does not publish physical population')
        if(count>0) call require(maxval(abs(work%midpoint(population)%parcels(1)%position- &
          group(population)%parcels(1)%position-2e-8_dp*work%midpoint(population)%parcels(1)%velocity))<1e-14_dp, &
          'old stage performs one full drift')
      end do
      old_group_source=group_source;old_group_remainder=group_remainder
      expected_group=work%midpoint;fraction=1
      call advance_mpi_eb_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,expected_group,1e-8_dp,0._dp,faces,expected_source,expected_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction)
      call require(ok,'single new half-kick reference')
      expected_force=0
      if(count>0) then
        parcel_weight=expected_group(1)%parcels(1)%mass*expected_group(1)%parcels(1)%multiplicity
        expected_force(2:4)=parcel_weight*dv
        expected_force(5)=parcel_weight*dot_product(expected_group(1)%parcels(1)%velocity+.5_dp*dv,dv)
        do population=1,3
          expected_group(population)%parcels(1)%velocity=expected_group(population)%parcels(1)%velocity+dv
        end do
      end if
      call MPI_Allreduce(MPI_IN_PLACE,expected_force,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      expected_force=expected_force+old_force
      if(rank==0) accel(1)=2*accel(1)
      call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,group,work,0._dp,2e-8_dp,1,2,faces,group_source,group_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction,acceleration=accel)
      call require(.not.ok.and.work%next_stage==1.and.all(work%force_change==old_force), &
        'changed rank-local acceleration retains stage and external ledger')
      accel=[1000._dp,-2000._dp,3000._dp]
      if(count>0) work%midpoint(3)%nodes=0
      call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,group,work,0._dp,2e-8_dp,1,2,faces,group_source,group_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction,acceleration=accel)
      call require(ok.and.work%next_stage==2,'first new SDC population iteration')
      call require(all(work%force_change==old_force),'intermediate iteration does not repeat external force')
      do population=1,3
        call require(all(group(population)%component_mass==before_group(population)%component_mass), &
          'intermediate iteration leaves population uncommitted')
      end do
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==old_group_source%fields(p)%values).and. &
          all(group_remainder%fields(p)%values==old_group_remainder%fields(p)%values), &
          'intermediate SDC iteration retains old source rates')
      end do
      old_group_source=group_source;old_group_remainder=group_remainder;wall=17;ledger=19
      if(count>0) work%midpoint(3)%nodes=0
      call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,group,work,0._dp,2e-8_dp,2,2,faces,group_source,group_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction,acceleration=accel)
      call require(.not.ok.and.work%next_stage==2.and.all(wall==17).and.all(ledger==19),'SDC late failure retains stage')
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==old_group_source%fields(p)%values).and. &
          all(group_remainder%fields(p)%values==old_group_remainder%fields(p)%values),'SDC late failure retains rates')
      end do
      work%midpoint(3)%nodes=query_nodes
      call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,group,work,0._dp,2e-8_dp,2,2,faces,group_source,group_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction,acceleration=accel)
      call require(ok.and.work%next_stage==0,'final SDC particle iteration commits')
      call require(all(work%force_change==expected_force),'final physical external momentum and work')
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==expected_source%fields(p)%values/1e-8_dp).and. &
          all(group_remainder%fields(p)%values==expected_remainder%fields(p)%values/1e-8_dp), &
          'only final SDC iteration computes the new half-kick rate')
      end do
      do population=1,3
        call require(all(group(population)%component_mass==expected_group(population)%component_mass), &
          'final iteration applies one new half kick')
        if(count>0) call require(all(group(population)%parcels(1)%velocity==expected_group(population)%parcels(1)%velocity) &
          .and.all(group(population)%parcels(1)%position==expected_group(population)%parcels(1)%position), &
          'final SDC particles match one drift and one new kick')
      end do
      call initialize_mpi_amr_sdc_populations(MPI_COMM_WORLD,tree,owners,group(1)%parcels, &
        group(1)%component_mass,zero,spacing,periodic,1,100,physical,virtuals,incoming,ok,phi)
      call require(ok,'recursive population initialization')
      fraction=0
      call advance_mpi_eb_recursive_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,0,1,100,physical,virtuals,incoming,hierarchy_work,0._dp,2e-8_dp,0,2,faces, &
        group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction)
      call require(ok,'recursive root old stage creates child ghosts')
      wall=17;ledger=19
      call advance_mpi_eb_recursive_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,physical,virtuals,incoming,hierarchy_work,0._dp,1e-8_dp,0,2,faces, &
        group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction)
      call require(.not.ok.and.hierarchy_work(2)%next_stage==0.and.all(wall==17).and.all(ledger==19), &
        'child rejects before parent interval completes')
      do stage=1,2
        call advance_mpi_eb_recursive_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
          1,0,1,100,physical,virtuals,incoming,hierarchy_work,0._dp,2e-8_dp,stage,2,faces, &
          group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction)
        call require(ok,'recursive root new stage')
      end do
      root_finished=physical(1)
      do substep=0,1
        do stage=0,2
          fraction=.5_dp*substep
          if(stage>0) fraction=fraction+.5_dp
          call advance_mpi_eb_recursive_population_stage(MPI_COMM_WORLD,species,transport,liquid,options, &
            tree,owners,g,gas,temp,1,1,1,100,physical,virtuals,incoming,hierarchy_work,substep*1e-8_dp,1e-8_dp,stage,2,faces, &
            group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction)
          call require(ok,'recursive child interval and stage')
        end do
        call require(hierarchy_work(2)%time+hierarchy_work(2)%interval==(substep+1)*1e-8_dp, &
          'child population clock advances once per interval')
        call require(size(incoming(2)%parcels)==count,'incoming ghosts retained across child substeps')
      end do
      call require(all(physical(1)%component_mass==root_finished%component_mass), &
        'child ghost feedback does not advance parent physical population again')
      wall=17;ledger=19;fraction=1
      call advance_mpi_eb_recursive_population_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,physical,virtuals,incoming,hierarchy_work,2e-8_dp,1e-8_dp,0,2,faces, &
        group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction)
      call require(.not.ok.and.all(wall==17).and.all(ledger==19),'child cannot advance beyond parent interval')
      ! Pinned SpraySB.cpp places each child at release location + full dt *
      ! its own velocity, although allocation occurs after the parent loop.
      child_options=options;child_options%breakup_model='khrt';child_options%khrt_children=.true.
      child_options%khrt_reference_diameter_cubed=1e-25_dp
      child_options%surface_tension=.072_dp;child_options%liquid_viscosity=.001_dp
      child_options%drag=.false.;child_options%heat=.false.;child_options%evaporation=.false.
      group=before_group;child_id=1;fraction=.25_dp
      birth_position=[1.75_dp,2.25_dp,2.25_dp]
      do population=1,3
        if(count==0) cycle
        group(population)%parcels(1)=spray_parcel(1_int64,birth_position, &
          [1000._dp,10._dp,5._dp],1e-12_dp,350._dp,100._dp)
        group(population)%component_mass=1e-12_dp;group(population)%nodes=query_nodes
      end do
      call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,child_options,tree,owners,g,gas,temp, &
        1,1,1,100,group,child_work,0._dp,2e-8_dp,0,2,faces,group_source,group_remainder,wall,ledger,ok, &
        phi,centroids,metrics,gas,ending,temp,temp,fraction,physical_source=physical_source,last_id=child_id)
      call require(ok.and.child_id==5,'old physical stage reserves exactly four child IDs')
      call require(size(group(1)%parcels)==count,'old stage keeps physical publication deferred')
      do population=2,3
        call require(size(child_work%midpoint(population)%parcels)==count,'copies create no shedding children')
      end do
      do p=1,size(child_work%midpoint(1)%parcels)
        if(child_work%midpoint(1)%parcels(p)%id==1) cycle
        call require(all(child_work%midpoint(1)%parcels(p)%position==birth_position+ &
          2e-8_dp*child_work%midpoint(1)%parcels(p)%velocity), &
          'child starts at release position plus its own full-step displacement')
      end do
      expected_group=child_work%midpoint
      do stage=1,2
        call advance_mpi_eb_population_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,child_options,tree,owners, &
          g,gas,temp,1,1,1,100,group,child_work,0._dp,2e-8_dp,stage,2,faces, &
          group_source,group_remainder,wall,ledger,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction, &
          physical_source=physical_source,last_id=child_id)
        call require(ok.and.child_id==5,'intermediate and new-time stages allocate no child IDs')
      end do
      do population=1,3
        call require(size(group(population)%parcels)==size(expected_group(population)%parcels), &
          'final particle populations contain no extra children')
        do p=1,size(group(population)%parcels)
          call require(all(group(population)%parcels(p)%position==expected_group(population)%parcels(p)%position), &
            'new-time half kick does not drift children twice')
        end do
      end do
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,g,gas,temp,film_halo,film_halo_t,MPI_COMM_WORLD,ok, &
        ghost_volume_fraction=film_volume,halo_width=1,target_volume_fraction=metrics, &
        boundary_faces=faces,restrict_snapshot=.false.)
      call require(ok,'frozen film stage halo')
      ! Exercise the same halo exchange followed by deposition/splash flight.
      child_options=options;child_options%surface_tension=.072_dp;child_options%liquid_viscosity=.001_dp
      child_options%drag=.false.;child_options%heat=.false.;child_options%evaporation=.false.
      film_options%temperature=300;film_options%maximum_films=100;film_options%enabled=.true.
      do trial=1,2
        allocate(impact_films(0),impact_mass(1,0));ghosts=before_group(1)%parcels;mass=before_group(1)%component_mass
        if(count>0) then
          ghosts(1)=spray_parcel(1_int64,[1.5_dp,1.250001_dp,1.5_dp], &
            [0._dp,-real(merge(1,100,trial==1),dp),0._dp],1e-9_dp,350._dp,100._dp)
          mass(:,1)=ghosts(1)%mass
        end if
        initial_total=0
        if(count>0) then
          initial_total(1)=ghosts(1)%mass*ghosts(1)%multiplicity
          initial_total(2)=initial_total(1)*ghosts(1)%velocity(2)
          initial_total(3)=initial_total(1)*(liquid_energy(liquid(1),350._dp)+.5_dp*sum(ghosts(1)%velocity**2))
        end if
        call MPI_Allreduce(MPI_IN_PLACE,initial_total,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        impact_id=1;wall=17;ledger=19;old_mass=mass;old_source=source
        ! Fail source routing after private impacts, then retry the same input.
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,child_options,tree,owners, &
          film_halo,film_halo_t,1,query_nodes,ghosts,mass,zero,spacing,periodic,1e-8_dp,source,ok, &
          film_volume,phi,centroids,2e-6_dp,faces,wall,ledger, &
          target_level=1,level_source=group_source,source_remainder=group_remainder, &
          eb_geometry=g,maximum_received=-1,last_id=impact_id,wall_films=impact_films,film_mass=impact_mass, &
          film_settings=film_options,model_rejected=rejected)
        call require(.not.ok.and..not.rejected.and.impact_id==1.and.size(impact_films)==0.and. &
          all(mass==old_mass).and.all(wall==17).and.all(ledger==19),'late source failure rolls back impact payloads')
        if(count>0) call require(all(ghosts(1)%position==[1.5_dp,1.250001_dp,1.5_dp]), &
          'late source failure retains pre-impact position')
        do p=1,2
          if(owners(p)/=rank) cycle
          call require(all(source%fields(p)%values==old_source%fields(p)%values),'late impact source rollback')
        end do
        call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquid,child_options,tree,owners, &
          film_halo,film_halo_t,1,query_nodes,ghosts,mass,zero,spacing,periodic,1e-8_dp,source,ok, &
          film_volume,phi,centroids,2e-6_dp,faces,wall,ledger, &
          last_id=impact_id,wall_films=impact_films,film_mass=impact_mass,film_settings=film_options,model_rejected=rejected)
        call require(ok.and..not.rejected,'halo impact flight')
        impact_total=0
        do p=1,size(ghosts)
          parcel_weight=ghosts(p)%mass*ghosts(p)%multiplicity
          impact_total(1)=impact_total(1)+parcel_weight
          impact_total(2)=impact_total(2)+parcel_weight*ghosts(p)%velocity(2)
          impact_total(3)=impact_total(3)+parcel_weight*(liquid_energy(liquid(1),ghosts(p)%temperature)+ &
            .5_dp*sum(ghosts(p)%velocity**2))
        end do
        do p=1,size(impact_films)
          impact_total(1)=impact_total(1)+impact_films(p)%mass
          impact_total(3)=impact_total(3)+impact_films(p)%mass*liquid_energy(liquid(1),impact_films(p)%temperature)
        end do
        call MPI_Allreduce(MPI_IN_PLACE,impact_total,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call require(all(abs(impact_total-initial_total-wall([1,3,5]))<= &
          1e-12_dp*max(abs(initial_total),1e-30_dp)),'halo impact gas particle film wall inventory')
        call require(all(ledger==0),'interior impact has no outflow')
        if(count>0) then
          call require(size(impact_films)>0,'halo contact deposits a physical film')
          if(trial==2) call require(size(ghosts)>1.and.impact_id>1,'halo splash creates children and IDs')
        end if
        deallocate(impact_films,impact_mass)
      end do

      call reactive_conserved_to_primitive(species,cell,900._dp,primitive,t,sound,ok)
      call require(ok,'frozen film reference gas')
      allocate(films(2*count),film_mass(1,2*count),probes(3,2*count),film_nodes(2*count))
      film_nodes=2;film_mass=1e-8_dp
      do p=1,2*count
        films(p)=spray_wall_film([1.5_dp,1.25_dp,1.5_dp],[0._dp,1._dp,0._dp], &
          1e-8_dp,350._dp,1e-3_dp,1e-8_dp/(1000*.25_dp*acos(-1._dp)*1e-6_dp),p)
        if(p==2) films(p)%position(1)=2.25_dp
        probes(:,p)=films(p)%position+[0._dp,1e-8_dp,0._dp]
      end do
      film_options%enabled=.true.;film_options%evolve=.true.;film_options%temperature=400
      film_options%liquid_conductivity=.6_dp;film_options%minimum_height=1e-15_dp
      reference_films=films;reference_film_mass=film_mass;expected_wall=0;delta=0
      do p=1,size(films)
        call advance_mixture_wall_film_frozen(reference_films(p),liquid,reference_film_mass(:,p),species,transport, &
          primitive,900._dp,options,.6_dp,400._dp,1e-10_dp,1e-15_dp,ledger,film_wall,ok)
        call require(ok,'independent film half-step reference')
        delta=delta+ledger;expected_wall=expected_wall-film_wall
      end do
      call MPI_Allreduce(MPI_IN_PLACE,delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,expected_wall,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call kick_mpi_eb_film_halo(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_nodes,probes,films,film_mass,zero,spacing,periodic,1e-10_dp, &
        film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
      call require(ok,'film half step with frozen gas halo')
      call require(maxval(abs(film_mass-reference_film_mass),mask=.true.)<1e-22_dp,'film component exchange reference')
      call require(abs(wall(5)-expected_wall)<1e-12_dp*abs(expected_wall),'film wall heat reference')
      do p=1,size(films)
        call require(all(films(p)%position==reference_films(p)%position),'film stays at wall')
        call require(abs(films(p)%temperature-reference_films(p)%temperature)<1e-10_dp,'film temperature reference')
      end do
      call integrate_mpi_eb_tree_3d(tree,owners,g,physical_source,integral,MPI_COMM_WORLD,ok)
      call require(ok.and.all(abs(integral-delta)<=1e-12_dp*max(abs(delta),1e-30_dp)), &
        'physical film exchange deposited on finest containing leaves')
      film_energy_change=0
      do p=1,size(films)
        film_energy_change=film_energy_change+films(p)%mass*liquid_energy(liquid(1),films(p)%temperature)- &
          1e-8_dp*liquid_energy(liquid(1),350._dp)
      end do
      call MPI_Allreduce(MPI_IN_PLACE,film_energy_change,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call require(abs(film_energy_change+integral(5)-wall(5))<1e-14_dp,'film gas and wall energy closure')
      ! The first film can update privately before material 2 leaves its range.
      kept_films=films;reference_film_mass=film_mass
      old_group_source=group_source;old_group_remainder=group_remainder;old_physical_source=physical_source
      film_options%port_conductivity_models(2)='polynomial'
      film_options%port_conductivity_coefficients(1,2)=.6_dp
      film_options%port_conductivity_temperature_ranges(:,2)=[300._dp,340._dp];wall=17
      call kick_mpi_eb_film_halo(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_nodes,probes,films,film_mass,zero,spacing,periodic,1e-10_dp, &
        film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
      call require(.not.ok.and.all(wall==17).and.all(film_mass==reference_film_mass),'late material failure rolls back film')
      do p=1,size(films)
        call require(films(p)%temperature==kept_films(p)%temperature.and.films(p)%mass==kept_films(p)%mass, &
          'film payload rollback')
      end do
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==old_group_source%fields(p)%values).and. &
          all(group_remainder%fields(p)%values==old_group_remainder%fields(p)%values).and. &
          all(physical_source%fields(p)%values==old_physical_source%fields(p)%values),'all film sources roll back')
      end do
      film_options%port_conductivity_models(2)=''
      film_options%port_conductivity_coefficients(:,2)=0
      film_options%port_conductivity_temperature_ranges(:,2)=0
      do population=1,3
        film_groups(population)%films=films;film_groups(population)%component_mass=film_mass
        film_groups(population)%nodes=film_nodes;film_groups(population)%probes=probes
      end do
      old_film_groups=film_groups
      reference_films=films;reference_film_mass=film_mass
      call kick_mpi_eb_film_halo(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_nodes,probes,reference_films,reference_film_mass,zero,spacing,periodic,1e-10_dp, &
        film_volume,phi,centroids,expected_film_source,expected_film_remainder,old_film_wall,ok,expected_physical_film)
      call require(ok,'single-population old film half-step reference')
      call advance_mpi_eb_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_groups,film_work,zero,spacing,periodic,0._dp,2e-10_dp,0,2, &
        film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
      call require(ok.and.film_work%next_stage==1.and.all(wall==old_film_wall),'old film stage physical-only wall ledger')
      do population=1,3
        call require(all(film_groups(population)%component_mass==old_film_groups(population)%component_mass), &
          'old film stage defers publication')
        call require(all(film_work%midpoint(population)%component_mass==reference_film_mass),'all film copies advance once')
      end do
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(abs(group_source%fields(p)%values-3*expected_film_source%fields(p)%values/1e-10_dp)<= &
          1e-12_dp*max(abs(group_source%fields(p)%values),1e-30_dp)),'combined film rates include all copies')
        call require(all(physical_source%fields(p)%values==expected_physical_film%fields(p)%values/1e-10_dp), &
          'physical film rates exclude virtual and ghost copies')
      end do
      old_group_source=group_source;old_group_remainder=group_remainder;old_physical_source=physical_source
      if(count>0) film_work%midpoint(3)%films(2)%material_id=17
      call advance_mpi_eb_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_groups,film_work,zero,spacing,periodic,0._dp,2e-10_dp,1,2, &
        film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
      call require(ok.and.film_work%next_stage==2,'intermediate film stage skips exchange')
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==old_group_source%fields(p)%values).and. &
          all(physical_source%fields(p)%values==old_physical_source%fields(p)%values),'intermediate film rates stay cached')
      end do
      wall=17
      call advance_mpi_eb_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_groups,film_work,zero,spacing,periodic,0._dp,2e-10_dp,2,2, &
        film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
      call require(.not.ok.and.film_work%next_stage==2.and.all(wall==17),'late film copy failure retains stage')
      do population=1,3
        call require(all(film_groups(population)%component_mass==old_film_groups(population)%component_mass), &
          'late film copy failure rolls back every population')
      end do
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(group_source%fields(p)%values==old_group_source%fields(p)%values).and. &
          all(group_remainder%fields(p)%values==old_group_remainder%fields(p)%values).and. &
          all(physical_source%fields(p)%values==old_physical_source%fields(p)%values),'late film copy source rollback')
      end do
      if(count>0) film_work%midpoint(3)%films(2)%material_id=2
      call kick_mpi_eb_film_halo(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_nodes,probes,reference_films,reference_film_mass,zero,spacing,periodic,1e-10_dp, &
        film_volume,phi,centroids,expected_film_source,expected_film_remainder,expected_film_wall,ok,expected_physical_film)
      call require(ok,'single-population new film half-step reference')
      call advance_mpi_eb_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        film_halo,film_halo_t,1,1,film_groups,film_work,zero,spacing,periodic,0._dp,2e-10_dp,2,2, &
        film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
      call require(ok.and.film_work%next_stage==0.and.all(wall==old_film_wall+expected_film_wall), &
        'final film stage commits both physical wall half steps once')
      do population=1,3
        call require(all(film_groups(population)%component_mass==reference_film_mass),'film final half-step state')
      end do
      ! Depletion removes film payloads, nodes and probes together before the new half step.
      film_groups=old_film_groups;film_options%minimum_height=1e-3_dp
      do stage=0,2
        call advance_mpi_eb_film_sdc_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
          film_halo,film_halo_t,1,1,film_groups,depletion_work,zero,spacing,periodic,0._dp,2e-10_dp,stage,2, &
          film_volume,phi,centroids,group_source,group_remainder,wall,ok,physical_source)
        call require(ok,'depleted film continuation')
      end do
      do population=1,3
        call require(size(film_groups(population)%films)==0.and.size(film_groups(population)%nodes)==0.and. &
          size(film_groups(population)%probes,2)==0.and.size(film_groups(population)%component_mass,2)==0, &
          'depleted film population and metadata compact together')
      end do
      film_options%minimum_height=1e-15_dp
      call initialize_mpi_eb_film_levels(MPI_COMM_WORLD,tree,owners,films,film_mass,zero,spacing,100,phi,film_physical,ok)
      call require(ok,'recursive film hierarchy initialization')
      allocate(film_virtuals(2),film_ghosts(2))
      do p=1,2
        film_virtuals(p)%films=film_physical(p)%films(:0)
        film_virtuals(p)%component_mass=film_physical(p)%component_mass(:,:0)
        film_virtuals(p)%nodes=film_physical(p)%nodes(:0);film_virtuals(p)%probes=film_physical(p)%probes(:,:0)
        film_ghosts(p)=film_virtuals(p)
      end do
      fraction=0
      call advance_mpi_eb_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        gas,temp,1,0,1,100,film_physical,film_virtuals,film_ghosts,film_clocks,0._dp,2e-10_dp,0,2,faces, &
        group_source,group_remainder,wall,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction,physical_source=physical_source)
      call require(ok.and.size(film_ghosts(2)%films)==count,'root old film stage creates child ghosts')
      wall=17
      call advance_mpi_eb_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        gas,temp,1,1,1,100,film_physical,film_virtuals,film_ghosts,film_clocks,0._dp,1e-10_dp,0,2,faces, &
        group_source,group_remainder,wall,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction,physical_source=physical_source)
      call require(.not.ok.and.film_clocks(2)%next_stage==0.and.all(wall==17),'film child waits for parent interval')
      do stage=1,2
        call advance_mpi_eb_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
          gas,temp,1,0,1,100,film_physical,film_virtuals,film_ghosts,film_clocks,0._dp,2e-10_dp,stage,2,faces, &
          group_source,group_remainder,wall,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction, &
          physical_source=physical_source)
        call require(ok,'root new film stages')
      end do
      root_finished_film=film_physical(1)
      call finish_mpi_eb_film_level(MPI_COMM_WORLD,0,0._dp,2e-10_dp,1,1,film_virtuals,film_ghosts,film_clocks,ok)
      call require(.not.ok,'root film cleanup waits for descendant synchronization')
      do substep=0,1
        do stage=0,2
          fraction=.5_dp*substep
          if(stage>0) fraction=fraction+.5_dp
          call advance_mpi_eb_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options, &
            tree,owners,g,gas,temp,1,1,1,100,film_physical,film_virtuals,film_ghosts,film_clocks, &
            substep*1e-10_dp,1e-10_dp,stage,2,faces,group_source,group_remainder,wall,ok, &
            phi,centroids,metrics,gas,ending,temp,temp,fraction,physical_source=physical_source)
          call require(ok,'recursive film child interval')
        end do
        call finish_mpi_eb_film_level(MPI_COMM_WORLD,1,substep*1e-10_dp,1e-10_dp,substep+1,2, &
          film_virtuals,film_ghosts,film_clocks,ok)
        call require(ok.and.size(film_virtuals(2)%films)==0,'film virtuals expire after every child interval')
        call require(size(film_ghosts(2)%films)==merge(count,0,substep==0),'film ghosts expire only after final child cycle')
      end do
      call require(all(film_physical(1)%component_mass==root_finished_film%component_mass), &
        'child film copies do not advance parent physical films again')
      wall=17
      call advance_mpi_eb_recursive_film_stage(MPI_COMM_WORLD,species,transport,liquid,options,film_options,tree,owners,g, &
        gas,temp,1,1,1,100,film_physical,film_virtuals,film_ghosts,film_clocks,2e-10_dp,1e-10_dp,0,2,faces, &
        group_source,group_remainder,wall,ok,phi,centroids,metrics,gas,ending,temp,temp,fraction,physical_source=physical_source)
      call require(.not.ok.and.all(wall==17),'film child cannot exceed parent clock')
      call finish_mpi_eb_film_level(MPI_COMM_WORLD,0,0._dp,2e-10_dp,1,1,film_virtuals,film_ghosts,film_clocks,ok)
      call require(ok.and.size(film_virtuals(1)%films)==0.and.size(film_ghosts(1)%films)==0,'root film cleanup')
      old_mass=mass;old_source=source;wall=17;ledger=19;fraction=1.1_dp
      call advance_mpi_eb_ghost_stage(MPI_COMM_WORLD,species,transport,liquid,options,tree,owners,g,gas,temp, &
        1,1,1,100,query_nodes,ghosts,mass,1e-8_dp,0._dp,faces,source,remainder,wall,ledger,ok,phi,centroids,metrics, &
        gas,ending,temp,temp,fraction)
      call require(.not.ok.and.all(mass==old_mass).and.all(wall==17).and.all(ledger==19),'bad stage time rollback')
      do p=1,2
        if(owners(p)==rank) call require(all(source%fields(p)%values==old_source%fields(p)%values),'stage-time source rollback')
      end do
    end if
    saved_values=values;saved_t=site_t;fraction=1.1_dp
    call sample_mpi_eb_tree_particle_gas(MPI_COMM_WORLD,species,tree,owners,g,gas,temp,1,query_nodes,query, &
      values,site_t,ok,gas,ending,temp,temp,fraction,target_volume_fraction=metrics,boundary_faces=faces, &
        nodal_phi=phi,centroids=centroids)
    call require(.not.ok.and.all(values==saved_values).and.all(site_t==saved_t),'temporal halo failure rollback')
    fraction=.25_dp
    if(count>0) query(2,1)=10
    call sample_mpi_eb_tree_particle_gas(MPI_COMM_WORLD,species,tree,owners,g,gas,temp,1,query_nodes,query, &
      values,site_t,ok,gas,ending,temp,temp,fraction,target_volume_fraction=metrics,boundary_faces=faces, &
        nodal_phi=phi,centroids=centroids)
    call require(.not.ok.and.all(values==saved_values).and.all(site_t==saved_t),'late sampling rollback')
  end subroutine check_temporal_particle_gas
  subroutine sample(limit)
    integer, intent(in) :: limit
    call sample_mpi_amr_tree_spray_field(MPI_COMM_WORLD,tree,owners,field,1,positions,zero,spacing, &
      periodic,samples,ok,limit,g,phi)
  end subroutine
  subroutine run()
    call eb_tree_spray_stencil(tree,g,phi,position,zero,spacing,cells,w,v,ok,periodic=periodic)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
