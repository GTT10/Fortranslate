program test_spray_mpi_amr_step
  use mpi_f08
  use spray_mpi_amr_step_mod, only: advance_mpi_amr_mixture_coupling,kick_mpi_amr_mixture_halo, &
    advance_mpi_amr_ghost_halo
  use spray_amr_tree_partition_mod, only: locate_amr_tree_spray_cells
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties,advance_mixture_frozen_gas
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patches(2)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: gas,temp,saved,saved_t,expected,expected_t,local,local_t
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel), allocatable :: drops(:),saved_drops(:)
  real(dp), allocatable :: masses(:,:),saved_mass(:,:)
  real(dp), allocatable :: root(:,:,:,:),primitive(:),cell(:),y(:),before(:),after(:),scale(:)
  real(dp), allocatable :: positions(:,:)
  integer, allocatable :: locations(:,:)
  real(dp) :: spacing(3)=.001_dp,zero(3)=0._dp,t,sound,r,interval
  integer :: nv,node,n(3),i,j,k,p,water,nitrogen,steps
  integer :: rank,nranks,ierr,owners(3),q
  logical :: ok
  integer(int64) :: last_id
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks==2.or.nranks==4,'ranks')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermodynamics')
  nv=reactive_nvar(size(species))
  allocate(root(nv,4,4,4),primitive(nv),cell(nv),y(size(species)),before(nv),after(nv),scale(nv))
  water=0;nitrogen=0
  do p=1,size(species)
    if (species(p)%name=='H2O') water=p
    if (species(p)%name=='N2') nitrogen=p
  end do
  call require(water>0.and.nitrogen>0,'species identity')
  y=0._dp;y(nitrogen)=1._dp
  r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0._dp;primitive(1)=1._dp;primitive(5)=r*1000._dp;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  do k=1,4
  do j=1,4
  do i=1,4
    root(:,i,j,k)=cell
  end do
  end do
  end do
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'fine patch')
  call initialize_amr_patch_3d(4,8,8,1,2,1,2,1,2,2,patches(2),ok)
  call require(ok,'nested patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  call prolong_amr_patch_tree_3d(tree,root,gas,ok)
  call require(ok,'gas fields')
  allocate(temp%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(temp%fields(node)%values(1,n(1),n(2),n(3)))
    temp%fields(node)%values=1000._dp
  end do
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water liquid')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide liquid')
  owners=[0,1,nranks-1]
  allocate(local%fields(3),local_t%fields(3))
  do node=1,3
    if (owners(node)/=rank) cycle
    local%fields(node)%values=gas%fields(node)%values
    local_t%fields(node)%values=temp%fields(node)%values
  end do
  ! All parcels begin on the last rank, including a coarse/fine crossing and wrap.
  allocate(drops(merge(3,0,rank==nranks-1)),masses(2,merge(3,0,rank==nranks-1)))
  do p=1,size(drops)
    drops(p)=spray_parcel(int(p,int64),spacing*[1.99_dp,.8_dp,.8_dp], &
      [20._dp,10._dp,5._dp],1e-12_dp,350._dp,100._dp)
    if (p==2) drops(p)%position=spacing*[.99_dp,.9_dp,.9_dp]
    if (p==3) drops(p)%position=spacing*[3.99_dp,1.1_dp,1.1_dp]
    masses(:,p)=drops(p)%mass*[.25_dp,.75_dp]
  end do
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses
  interval=2.e-5_dp;last_id=3
  options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0,'late substep exhaustion')
  call unchanged()
  options%maximum_substeps=100000
  interval=merge(2.e-5_dp,1.e-6_dp,rank==0)
  call advance()
  call require(.not.ok,'interval disagreement');call unchanged()
  interval=2.e-5_dp
  options%maximum_parcels=1
  call advance()
  call require(.not.ok,'late migration capacity failure');call unchanged()
  options%maximum_parcels=100
  call totals(before)
  call advance()
  call require(ok.and.steps>1,'coupled substeps')
  call totals(after)
  scale=max(abs(before),abs(after),1.e-30_dp)
  call require(all(abs(after-before)<=3.e-12_dp*scale),'composite gas plus liquid conservation')
  allocate(positions(3,size(drops)))
  do p=1,size(drops)
    call require(drops(p)%temperature>350._dp.and.drops(p)%mass<1e-12_dp,'heating and evaporation active')
    call require(drops(p)%velocity(1)<20._dp,'drag active')
    positions(:,p)=drops(p)%position
    call require(all(positions(:,p)>=0).and.all(positions(:,p)<4*spacing),'periodic wrapped trajectory')
  end do
  call locate_amr_tree_spray_cells(tree,positions,zero,spacing,[.true.,.true.,.true.],locations,ok)
  call require(ok,'final locations')
  do p=1,size(drops)
    call require(owners(locations(1,p))==rank,'finest patch migration')
  end do
  ! Empty ranks participate in both field queries and final migration.
  q=size(drops)
  call MPI_Allreduce(MPI_IN_PLACE,q,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(q==3,'global parcel count')
  call check_level_view()
  call check_halo_kick()
  ! Exercise collective child IDs and the transverse-energy gas debit.
  local=saved;local_t=saved_t;drops=saved_drops;masses=saved_mass
  do p=1,size(drops)
    drops(p)%velocity=[1000._dp,10._dp,5._dp]
  end do
  options%breakup_model='khrt';options%khrt_children=.true.
  options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
  options%khrt_reference_diameter_cubed=1e-25_dp
  options%drag=.false.;options%heat=.false.;options%evaporation=.false.
  interval=1.e-8_dp;last_id=3
  call totals(before)
  call advance()
  call require(ok.and.last_id>3,'KH children and collective IDs')
  call totals(after)
  scale=max(abs(before),abs(after),1.e-30_dp)
  call require(all(abs(after-before)<=3.e-12_dp*scale),'KH coupled conservation')
  q=size(drops)
  call MPI_Allreduce(MPI_IN_PLACE,q,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call require(q>3,'KH child migration')
  deallocate(drops,masses);allocate(drops(0),masses(2,0))
  options=spray_options();interval=1.e-6_dp
  saved=local;saved_t=local_t;saved_drops=drops;saved_mass=masses;last_id=3
  call advance()
  call require(ok,'globally empty population');call unchanged()
  if (rank==0) print *, 'PASS: periodic MPI AMR mixed-particle coupling'
  call MPI_Finalize(ierr)
contains
  subroutine check_halo_kick()
    type(amr_tree_state_3d) :: halo,ht,source,kept_source,routed,surplus,kept_routed,kept_surplus
    type(spray_parcel), allocatable :: hd(:),reference(:)
    real(dp), allocatable :: hm(:,:),rm(:,:)
    real(dp) :: delta(nv),expected_delta(nv),integral(nv),hp(nv),guess,cs,dt
    integer, allocatable :: nodes(:)
    real(dp) :: wall(5),exterior(nv),position_before(3),drift_dt
    character(len=16) :: faces(6)
    integer :: a,b,c,d,sz(3),count
    count=merge(1,0,rank==owners(2))
    allocate(hd(count),hm(2,count),nodes(count),halo%fields(3),ht%fields(3))
    nodes=2;dt=1e-8_dp
    do a=1,count
      ! This point lies just beyond the fine patch's valid x extent.
      hd(a)=spray_parcel(1_int64,spacing*[2.1_dp,1.5_dp,1.5_dp], &
        [20._dp,10._dp,5._dp],1e-12_dp,350._dp,100._dp)
      hm(:,a)=hd(a)%mass*[.25_dp,.75_dp]
    end do
    do a=1,3
      if(owners(a)/=rank) cycle
      sz=tree%nodes(a)%shape
      allocate(halo%fields(a)%values(nv,0:sz(1)+1,0:sz(2)+1,0:sz(3)+1))
      allocate(ht%fields(a)%values(1,0:sz(1)+1,0:sz(2)+1,0:sz(3)+1))
      do d=0,sz(3)+1
      do c=0,sz(2)+1
      do b=0,sz(1)+1
        halo%fields(a)%values(:,b,c,d)=cell
      end do
      end do
      end do
      ht%fields(a)%values=1000._dp
    end do
    reference=hd;rm=hm;expected_delta=0
    do a=1,count
      call reactive_conserved_to_primitive(species,cell,1000._dp,hp,guess,cs,ok)
      call require(ok,'halo kick EOS reference')
      call advance_mixture_frozen_gas(species,transport,hp,guess,reference(a),liquids,rm(:,a),options,dt,delta,ok)
      call require(ok,'halo kick reference physics')
      reference(a)%position=hd(a)%position
      expected_delta=expected_delta+delta
    end do
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],dt,source,ok, &
      target_level=1,level_source=routed,source_remainder=surplus,route_to_coarse=.true.)
    call require(ok,'halo kick accepts fine ghost particle')
    call require(all(hm==rm),'halo kick components')
    do a=1,count
      call require(all(hd(a)%position==reference(a)%position).and. &
        maxval(abs(hd(a)%velocity-reference(a)%velocity))<1e-13_dp.and. &
        abs(hd(a)%temperature-reference(a)%temperature)<1e-10_dp.and. &
        hd(a)%mass==reference(a)%mass,'halo kick physics and fixed position')
    end do
    integral=0
    do a=1,3
      if(owners(a)/=rank) cycle
      do d=lbound(source%fields(a)%values,4),ubound(source%fields(a)%values,4)
      do c=lbound(source%fields(a)%values,3),ubound(source%fields(a)%values,3)
      do b=lbound(source%fields(a)%values,2),ubound(source%fields(a)%values,2)
        integral=integral+source%fields(a)%values(:,b,c,d)*product(spacing/real(2**tree%nodes(a)%level,dp))
      end do
      end do
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,integral,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,expected_delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(integral-expected_delta)<=1e-12_dp*max(abs(expected_delta),1e-30_dp)), &
      'halo kick conserves extensive reaction')
    integral=0
    do a=1,3
      if(owners(a)/=rank) cycle
      sz=tree%nodes(a)%shape
      do d=1,sz(3)
      do c=1,sz(2)
      do b=1,sz(1)
        integral=integral+routed%fields(a)%values(:,b,c,d)*product(spacing/real(2**tree%nodes(a)%level,dp))
      end do
      end do
      end do
      call require(all(surplus%fields(a)%values==0),'routed halo kick leaves no in-domain remainder')
    end do
    call MPI_Allreduce(MPI_IN_PLACE,integral,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call require(all(abs(integral-expected_delta)<=1e-12_dp*max(abs(expected_delta),1e-30_dp)), &
      'routed coarse and fine kick conserves reaction')
    if(rank==owners(1)) call require(any(routed%fields(1)%values/=0),'fine ghost reaction reaches coarse owner')
    kept_source=source;kept_routed=routed;kept_surplus=surplus;reference=hd;rm=hm
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],dt,source,ok, &
      target_level=1,level_source=routed,source_remainder=surplus,route_to_coarse=.true.,maximum_received=0)
    call require(.not.ok.and.all(hm==rm),'transfer capacity failure rolls back kick components')
    do a=1,count
      call require(all(hd(a)%velocity==reference(a)%velocity).and.hd(a)%temperature==reference(a)%temperature, &
        'transfer capacity failure rolls back particle kick')
    end do
    do a=1,3
      if(owners(a)/=rank) cycle
      call require(all(source%fields(a)%values==kept_source%fields(a)%values),'transfer failure preserves raw source')
      call require(all(routed%fields(a)%values==kept_routed%fields(a)%values),'transfer failure preserves routed source')
      call require(all(surplus%fields(a)%values==kept_surplus%fields(a)%values),'transfer failure preserves remainder')
    end do
    ! The drift is independently timed; zero kick isolates periodic motion.
    drift_dt=2e-4_dp;faces='periodic';reference=hd;rm=hm
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],0._dp,source,ok, &
      drift_interval=drift_dt,boundary_faces=faces,wall_change=wall,boundary_change=exterior)
    call require(ok.and.all(wall==0).and.all(exterior==0),'halo periodic drift has no external exchange')
    do a=1,count
      position_before=modulo(reference(a)%position+drift_dt*reference(a)%velocity,4*spacing)
      call require(maxval(abs(hd(a)%position-position_before))<1e-15_dp,'halo kick-drift periodic position')
      call require(all(hd(a)%velocity==reference(a)%velocity).and.all(hm==rm),'pure drift preserves inventory')
    end do
    kept_source=source;reference=hd;rm=hm;wall=17;exterior=19
    ! More than 32 crossings rejects after the private kick has been computed.
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],dt,source,ok, &
      drift_interval=1._dp,boundary_faces=faces,wall_change=wall,boundary_change=exterior)
    call require(.not.ok.and.all(wall==17).and.all(exterior==19),'late drift rejection preserves ledgers')
    call require(all(hm==rm),'late drift rejection preserves components')
    do a=1,count
      call require(all(hd(a)%position==reference(a)%position).and.all(hd(a)%velocity==reference(a)%velocity), &
        'late drift rejection preserves particle position and velocity')
    end do
    do a=1,3
      if(owners(a)==rank) call require(all(source%fields(a)%values==kept_source%fields(a)%values),'late source rollback')
    end do
    expected_delta=0
    do a=1,count
      expected_delta(1)=expected_delta(1)-hd(a)%mass*hd(a)%multiplicity
      do b=1,2
        expected_delta(5+liquids(b)%vapor_index)=expected_delta(5+liquids(b)%vapor_index)-hm(b,a)*hd(a)%multiplicity
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,expected_delta,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    faces='outflow'
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,nodes,hd,hm,zero,spacing,[.false.,.false.,.false.],0._dp,source,ok, &
      drift_interval=drift_dt,boundary_faces=faces,wall_change=wall,boundary_change=exterior)
    call require(ok.and.all(hm==0).and.all(hd%mass==0),'halo drift retires outflow components')
    call require(exterior(1)==expected_delta(1).and.all(exterior(6:)==expected_delta(6:)), &
      'halo drift exports mass and each vapor species')
    ! Complete stage must retain a ghost which moves into the fine valid box.
    faces='periodic';nodes=2;drift_dt=1e-5_dp
    do a=1,count
      hd(a)=spray_parcel(1_int64,spacing*[2.1_dp,1.5_dp,1.5_dp], &
        [-20._dp,0._dp,0._dp],1e-12_dp,350._dp,100._dp)
      hm(:,a)=hd(a)%mass*[.25_dp,.75_dp]
    end do
    reference=hd;rm=hm;kept_routed=routed;kept_surplus=surplus;wall=17;exterior=19
    call advance_mpi_amr_ghost_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,1,1,0,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],dt,drift_dt,faces,routed,surplus,wall,exterior,ok)
    call require(.not.ok.and.all(wall==17).and.all(exterior==19),'ghost stage final capacity failure preserves ledgers')
    call require(size(hd)==count.and.all(hm==rm).and.all(nodes==2),'ghost stage capacity rollback')
    do a=1,count
      call require(all(hd(a)%position==reference(a)%position).and.all(hd(a)%velocity==reference(a)%velocity), &
        'ghost stage rolls back kick and drift')
    end do
    do a=1,3
      if(owners(a)/=rank) cycle
      call require(all(routed%fields(a)%values==kept_routed%fields(a)%values),'ghost stage rolls back reaction transfer')
      call require(all(surplus%fields(a)%values==kept_surplus%fields(a)%values),'ghost stage rolls back remainder')
    end do
    call advance_mpi_amr_ghost_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,1,1,100,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],dt,drift_dt,faces,routed,surplus,wall,exterior,ok)
    call require(ok.and.size(hd)==count.and.all(nodes==2),'ghost stage keeps a copy entering valid fine cells')
    call require(all(wall==0).and.all(exterior==0),'interior ghost stage has no boundary exchange')
    if(rank==owners(1)) call require(all(routed%fields(1)%values==0),'ghost stage excludes coarse reaction')
    if(rank==owners(2)) then
      call require(any(routed%fields(2)%values/=0),'ghost stage keeps valid fine reaction')
      call require(any(surplus%fields(2)%values/=0),'ghost stage retains outside reaction separately')
    end if
    do a=1,count
      call require(hd(a)%position(1)<2*spacing(1),'ghost enters fine valid region')
      call require(maxval(abs(hd(a)%position-reference(a)%position-drift_dt*hd(a)%velocity))<1e-15_dp, &
        'ghost stage uses kicked velocity for full drift')
      call require(hd(a)%mass<reference(a)%mass.and.hd(a)%temperature>reference(a)%temperature, &
        'ghost stage retains evaporation and heating')
    end do
    kept_source=source;reference=hd;rm=hm
    ! A failed query happens before publication of either population or source.
    if(count>0) nodes=0
    call kick_mpi_amr_mixture_halo(MPI_COMM_WORLD,species,transport,liquids,options,tree,owners,halo,ht, &
      1,nodes,hd,hm,zero,spacing,[.true.,.true.,.true.],dt,source,ok)
    call require(.not.ok.and.all(hm==rm),'halo rejection preserves components')
    do a=1,count
      call require(all(hd(a)%velocity==reference(a)%velocity).and.hd(a)%mass==reference(a)%mass, &
        'halo rejection preserves parcels')
    end do
    do a=1,3
      if(owners(a)==rank) call require(all(source%fields(a)%values==kept_source%fields(a)%values),'halo source rollback')
    end do
  end subroutine check_halo_kick
  subroutine check_level_view()
    type(amr_patch_tree_3d) :: coarse
    type(amr_patch_3d) :: no_patches(0)
    type(amr_tree_state_3d) :: capped,ct,flat,ft,cs,fs
    type(spray_parcel), allocatable :: cd(:),fd(:)
    real(dp), allocatable :: cm(:,:),fm(:,:)
    integer :: a,ns,flat_steps
    integer(int64) :: cid,fid
    call initialize_amr_patch_tree_3d([4,4,4],[integer::],no_patches,coarse,ok)
    call require(ok,'single-level reference tree')
    capped=saved;ct=saved_t;cd=saved_drops;cm=saved_mass
    fd=cd;fm=cm;cid=3;fid=3
    allocate(flat%fields(1),ft%fields(1))
    if(rank==owners(1)) then
      flat%fields(1)%values=saved%fields(1)%values
      ft%fields(1)%values=saved_t%fields(1)%values
    end if
    ! Fine velocities differ from the coarse state: any accidental fine donor
    ! changes particle drag and invalidates the single-level comparison.
    do a=2,3
      if(owners(a)/=rank) cycle
      capped%fields(a)%values(2,:,:,:)=5._dp
      capped%fields(a)%values(5,:,:,:)=capped%fields(a)%values(5,:,:,:)+12.5_dp
    end do
    cs=capped
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,capped,ct, &
      tree,owners,cd,cm,zero,spacing,interval,ns,ok,cid,gas_change=fs,maximum_level=0)
    call require(ok,'coarse-view particle evolution')
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,flat,ft, &
      coarse,owners(:1),fd,fm,zero,spacing,interval,flat_steps,ok,fid)
    call require(ok.and.ns==flat_steps,'coarse-view substeps match single-level evolution')
    call require(size(cd)==size(fd).and.cid==fid,'coarse ownership and IDs')
    call require(all(cm==fm),'coarse component evolution')
    do a=1,size(cd)
      call require(all(cd(a)%position==fd(a)%position).and.all(cd(a)%velocity==fd(a)%velocity).and. &
        cd(a)%temperature==fd(a)%temperature.and.cd(a)%mass==fd(a)%mass,'coarse particle evolution')
    end do
    if(rank==owners(1)) then
      call require(all(capped%fields(1)%values==flat%fields(1)%values),'coarse gas matches single-level update')
      call require(all(ct%fields(1)%values==ft%fields(1)%values),'coarse temperature')
    end if
    do a=1,3
      if(owners(a)/=rank) cycle
      call require(all(abs(capped%fields(a)%values-cs%fields(a)%values-fs%fields(a)%values)<= &
        3.e-13_dp*max(1._dp,abs(capped%fields(a)%values))),'capped source reconstruction')
      if(a==1) cycle
      call require(all(capped%fields(a)%values==cs%fields(a)%values),'finer gas untouched')
      call require(all(ct%fields(a)%values==saved_t%fields(a)%values),'finer temperature untouched')
      call require(all(fs%fields(a)%values==0),'finer source zero')
    end do
    cs=capped
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,capped,ct, &
      tree,owners,cd,cm,zero,spacing,interval,ns,ok,cid,maximum_level=merge(0,1,rank==0))
    call require(.not.ok,'rank-disagreeing level rejected')
    do a=1,3
      if(owners(a)==rank) call require(all(capped%fields(a)%values==cs%fields(a)%values),'level rejection rollback')
    end do
  end subroutine check_level_view
  subroutine advance()
    type(amr_tree_state_3d) :: previous,source
    integer :: a
    previous=local;source=local
    do a=1,size(owners)
      if(owners(a)==rank) source%fields(a)%values=17
    end do
    call advance_mpi_amr_mixture_coupling(MPI_COMM_WORLD,species,transport,liquids,options,local,local_t, &
      tree,owners,drops,masses,zero,spacing,interval,steps,ok,last_id,gas_change=source)
    do a=1,size(owners)
      call require(allocated(source%fields(a)%values).eqv.(owners(a)==rank),'owner-local particle source')
      if(owners(a)/=rank) cycle
      if(ok) then
        call require(all(abs(local%fields(a)%values-previous%fields(a)%values-source%fields(a)%values)<= &
          3.e-13_dp*max(1._dp,abs(local%fields(a)%values))),'source reconstructs accepted particle feedback')
      else
        call require(all(source%fields(a)%values==17),'particle source rollback')
      end if
    end do
  end subroutine
  subroutine unchanged()
    integer :: a
    call require(size(drops)==size(saved_drops).and.last_id==3,'particle count and IDs rollback')
    call require(all(masses==saved_mass),'component rollback')
    do a=1,size(drops)
      call require(all(drops(a)%position==saved_drops(a)%position).and. &
        all(drops(a)%velocity==saved_drops(a)%velocity).and.drops(a)%mass==saved_drops(a)%mass.and. &
        drops(a)%temperature==saved_drops(a)%temperature,'particle rollback')
    end do
    do a=1,3
      if (owners(a)/=rank) cycle
      call require(all(local%fields(a)%values==saved%fields(a)%values),'gas rollback')
      call require(all(local_t%fields(a)%values==saved_t%fields(a)%values),'temperature rollback')
    end do
  end subroutine
  subroutine totals(total)
    real(dp), intent(out) :: total(:)
    type(amr_tree_state_3d) :: assembled
    real(dp) :: liquid_total(nv),rho,cp,e,fractions(2),amount
    integer :: a,b,c
    assembled=gas
    do a=1,3
      assembled%fields(a)%values=0
      if (owners(a)==rank) assembled%fields(a)%values=local%fields(a)%values
      call MPI_Allreduce(MPI_IN_PLACE,assembled%fields(a)%values,size(assembled%fields(a)%values), &
        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    end do
    call composite_integrals_amr_patch_tree_3d(tree,assembled,spacing,total,ok)
    call require(ok,'gas composite integral')
    liquid_total=0
    do a=1,size(drops)
      fractions=masses(:,a)/drops(a)%mass
      call liquid_mixture_properties(liquids,fractions,drops(a)%temperature,rho,cp,e,ok)
      call require(ok,'liquid energy')
      amount=drops(a)%mass*drops(a)%multiplicity
      liquid_total(1)=liquid_total(1)+amount
      liquid_total(2:4)=liquid_total(2:4)+amount*drops(a)%velocity
      liquid_total(5)=liquid_total(5)+amount*(e+.5_dp*sum(drops(a)%velocity**2))
      do b=1,2
        c=5+liquids(b)%vapor_index
        liquid_total(c)=liquid_total(c)+masses(b,a)*drops(a)%multiplicity
      end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE,liquid_total,nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    total=total+liquid_total
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL rank ',rank,': ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
