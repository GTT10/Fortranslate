program test_spray_eb_step
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_eb_step_mod, only: advance_eb_mixture_coupling,advance_eb_mixture_coupling_allocated
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(eb_geometry_3d) :: g
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(spray_liquid) :: liquids(2)
  type(spray_options) :: options
  type(spray_parcel) :: drops(1),saved_drop
  real(dp) :: phi(0:4,0:4,0:4),t,sound,r,temperature(4,4,4),saved_t(4,4,4),masses(2,1),saved_mass(2,1)
  real(dp) :: ledger(5),interval=2e-5_dp
  real(dp), allocatable :: primitive(:),cell(:),y(:),state(:,:,:,:),saved(:,:,:,:),before(:),after(:),expected(:)
  integer :: i,j,k,nv,nitrogen,steps
  logical :: ok
  do i=0,4
    phi(i,:,:)=real(i,dp)*.001_dp-.00125_dp
  end do
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],g,ok)
  call require(ok,'geometry')
  call load_h2o2_full_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_full_transport(transport,ok)
  call require(ok,'transport')
  call initialize_spray_liquid(species,'H2O',1000._dp,4184._dp,373.15_dp,101325._dp,2.257e6_dp,liquids(1),ok)
  call require(ok,'water')
  call initialize_spray_liquid(species,'H2O2',800._dp,2000._dp,350._dp,101325._dp,1e6_dp,liquids(2),ok)
  call require(ok,'peroxide')
  nv=reactive_nvar(size(species))
  allocate(primitive(nv),cell(nv),y(size(species)),state(nv,4,4,4),before(nv),after(nv),expected(nv))
  nitrogen=0
  do i=1,size(species)
    if(species(i)%name=='N2') nitrogen=i
  end do
  call require(nitrogen>0,'nitrogen identity')
  y=0;y(nitrogen)=1;r=mixture_specific_gas_constant(species,y,ok)
  call require(ok,'gas constant')
  primitive=0;primitive(1)=1;primitive(2:4)=[1._dp,2._dp,3._dp];primitive(5)=r*1000;primitive(6:)=y
  call reactive_primitive_to_conserved(species,primitive,cell,t,sound,ok)
  call require(ok,'initial gas')
  do k=1,4
  do j=1,4
  do i=1,4
    state(:,i,j,k)=cell
  end do
  end do
  end do
  state(:,1,:,:)=0;temperature=t;temperature(1,:,:)=42
  drops(1)%id=1;drops(1)%position=[.001251_dp,.0015_dp,.0015_dp]
  drops(1)%velocity=[-20._dp,10._dp,5._dp];drops(1)%mass=1e-12_dp
  drops(1)%temperature=350;drops(1)%multiplicity=100
  masses(:,1)=drops(1)%mass*[.25_dp,.75_dp]
  saved=state;saved_t=temperature;saved_drop=drops(1);saved_mass=masses
  ledger=42;options%maximum_substeps=1
  call advance()
  call require(.not.ok.and.steps==0,'whole-interval failure')
  call require(all(state==saved).and.all(temperature==saved_t).and.all(masses==saved_mass),'gas and composition rollback')
  call require(all(drops(1)%position==saved_drop%position).and.all(ledger==42),'particle and wall ledger rollback')
  options%maximum_substeps=100000
  call inventory(before)
  call advance()
  call require(ok.and.steps>1,'coupled substeps')
  call inventory(after)
  expected=0;expected(1:5)=ledger
  call require(all(abs(after-before-expected)<=3e-12_dp*max(abs(before),abs(after),1e-30_dp)), &
    'gas plus particles minus wall ledger conservation')
  call require(drops(1)%temperature>350.and.drops(1)%mass<1e-12_dp,'heating and evaporation')
  call require(drops(1)%velocity(1)>0.and.drops(1)%velocity(2)<10,'reflection and drag')
  call require(ledger(2)>0.and.drops(1)%position(1)>.00125_dp,'wall impulse and fluid endpoint')
  call require(all(state(:,1,:,:)==0).and.all(temperature(1,:,:)==42),'covered cells unchanged')
  call check_fixed()
  call check_reference_etab()
  print *, 'PASS: EB mixture coupling with heat, evaporation, drag and reflection'
contains
  subroutine check_fixed()
    real(dp) :: constraint(5),initial(nv),final(nv),external(nv)
    type(spray_options) :: old_options
    integer :: attempt
    old_options=options;options%fixed_parts=.true.;options%mixture_time_scheme='euler'
    state=saved;temperature=saved_t;drops(1)=saved_drop;masses=saved_mass
    call inventory(initial)
    do attempt=1,2
      constraint=37;ledger=42
      options%maximum_substeps=merge(1,100000,attempt==1)
      call advance_eb_mixture_coupling(species,transport,liquids,options,g,phi,state,temperature, &
        drops,masses,interval,steps,ledger,ok,force_change=constraint)
      if(attempt==1) then
        call require(.not.ok.and.all(constraint==37).and.all(ledger==42),'fixed EB rejected ledger rollback')
        call require(all(state==saved).and.all(temperature==saved_t).and.all(masses==saved_mass), &
          'fixed EB rejected gas and mass rollback')
        call require(all(drops(1)%position==saved_drop%position).and.all(drops(1)%velocity==saved_drop%velocity), &
          'fixed EB rejected motion rollback')
      else
        call require(ok.and.steps>1,'fixed EB accepted substeps')
      end if
    end do
    call inventory(final);external=0;external(:5)=constraint+ledger
    call require(all(abs(final-initial-external)<=3e-12_dp*max(abs(initial),abs(final),1e-30_dp)), &
      'fixed EB gas parcel constraint closure')
    call require(all(drops(1)%position==saved_drop%position).and.all(drops(1)%velocity==0), &
      'fixed EB stationary particle')
    call require(drops(1)%temperature>350.and.drops(1)%mass<saved_drop%mass,'fixed EB thermal and mass exchange')
    options=old_options
  end subroutine check_fixed
  subroutine check_reference_etab()
    type(spray_parcel), allocatable :: particles(:),original(:)
    real(dp), allocatable :: composition(:,:),original_mass(:,:),grid(:,:,:,:),grid_t(:,:,:)
    real(dp) :: wall_change(5),initial(nv),final(nv)
    integer(int64) :: watermark
    integer :: count,q
    logical :: accepted,dispersed
    particles=[spray_parcel(1_int64,[.0025_dp,.0015_dp,.0015_dp], &
      [50._dp,0._dp,0._dp],1e-12_dp,350._dp,16._dp)]
    allocate(composition(2,1));composition(:,1)=particles(1)%mass*[.25_dp,.75_dp]
    original=particles;original_mass=composition;grid=saved;grid_t=saved_t
    options%breakup_model='etab';options%surface_tension=.072_dp;options%liquid_viscosity=.001_dp
    options%maximum_multiplicity=1;options%etab_parcel_factor=.5_dp;options%maximum_parcels=3
    options%drag=.false.;options%heat=.false.;options%evaporation=.false.
    watermark=1;wall_change=42
    call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,grid,grid_t, &
      particles,composition,1e-9_dp,count,wall_change,accepted,last_id=watermark)
    call require(.not.accepted.and.watermark==1.and.all(grid==saved).and. &
      all(grid_t==saved_t).and.all(particles%id==original%id).and.all(composition==original_mass), &
      'EB reference ETAB capacity rollback')
    options%maximum_parcels=10;wall_change=0
    call eb_inventory(grid,particles,composition,initial)
    call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,grid,grid_t, &
      particles,composition,1e-9_dp,count,wall_change,accepted,last_id=watermark)
    call require(accepted.and.size(particles)==4.and.watermark==5, &
      'EB reference ETAB child replacement')
    call require(all(particles%id==[2_int64,3_int64,4_int64,5_int64]).and. &
      all(particles%parent_id==1_int64),'EB reference ETAB child lineage')
    call eb_inventory(grid,particles,composition,final)
    call require(all(abs(final-initial-wall_change)<=3e-12_dp*max(abs(initial),abs(final),1e-30_dp)), &
      'EB reference ETAB gas and child conservation')
    particles=original;particles(1)%deformation=1.1_dp
    composition=original_mass;grid=saved;grid_t=saved_t
    options%maximum_parcels=10000;watermark=1;wall_change=0
    call eb_inventory(grid,particles,composition,initial)
    call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,grid,grid_t, &
      particles,composition,1e-9_dp,count,wall_change,accepted,last_id=watermark)
    call require(accepted.and.size(particles)>1.and.watermark>int(size(particles),int64), &
      'EB active ETAB parent replacement')
    dispersed=.false.
    do q=2,size(particles)
      dispersed=dispersed.or.any(abs(particles(q)%velocity-particles(1)%velocity)>1e-6_dp)
    enddo
    call require(dispersed,'EB active ETAB transverse dispersion')
    call eb_inventory(grid,particles,composition,final)
    call require(all(abs(final-initial-wall_change)<=3e-12_dp*max(abs(initial),abs(final),1e-30_dp)), &
      'EB active ETAB compensation conservation')
    particles=original;particles(1)%deformation=1.1_dp
    composition=original_mass;grid=saved;grid_t=saved_t
    options%etab_parcel_factor=0;watermark=1;wall_change=0
    call eb_inventory(grid,particles,composition,initial)
    call advance_eb_mixture_coupling_allocated(species,transport,liquids,options,g,phi,grid,grid_t, &
      particles,composition,1e-9_dp,count,wall_change,accepted,last_id=watermark)
    call require(accepted.and.size(particles)>1.and.watermark==1+int(size(particles),int64), &
      'EB active ETAB factor-zero replacement')
    dispersed=.false.
    do q=2,size(particles)
      dispersed=dispersed.or.any(abs(particles(q)%velocity-particles(1)%velocity)>1e-6_dp)
    enddo
    call require(dispersed,'EB active ETAB factor-zero transverse dispersion')
    call eb_inventory(grid,particles,composition,final)
    call require(all(abs(final-initial-wall_change)<=3e-12_dp*max(abs(initial),abs(final),1e-30_dp)), &
      'EB active ETAB factor-zero compensation conservation')
  end subroutine
  subroutine eb_inventory(grid,particles,composition,total)
    real(dp), intent(in) :: grid(:,:,:,:),composition(:,:)
    type(spray_parcel), intent(in) :: particles(:)
    real(dp), intent(out) :: total(:)
    real(dp) :: rho,cp,e,amount
    integer :: ii,jj,kk,q,c
    total=0
    do kk=1,g%nz
    do jj=1,g%ny
    do ii=1,g%nx
      total=total+grid(:,ii,jj,kk)*g%volume_fraction(ii,jj,kk)*g%dx*g%dy*g%dz
    enddo
    enddo
    enddo
    do q=1,size(particles)
      if(particles(q)%mass<=0) cycle
      call liquid_mixture_properties(liquids,composition(:,q)/particles(q)%mass, &
        particles(q)%temperature,rho,cp,e,ok)
      call require(ok,'EB ETAB liquid inventory')
      amount=particles(q)%mass*particles(q)%multiplicity
      total(1)=total(1)+amount;total(2:4)=total(2:4)+amount*particles(q)%velocity
      total(5)=total(5)+amount*(e+.5_dp*sum(particles(q)%velocity**2))
      do c=1,size(liquids)
        total(5+liquids(c)%vapor_index)=total(5+liquids(c)%vapor_index)+ &
          composition(c,q)*particles(q)%multiplicity
      enddo
    enddo
  end subroutine
  subroutine advance()
    call advance_eb_mixture_coupling(species,transport,liquids,options,g,phi,state,temperature, &
      drops,masses,interval,steps,ledger,ok)
  end subroutine
  subroutine inventory(total)
    real(dp), intent(out) :: total(:)
    real(dp) :: rho,cp,e,amount
    integer :: ii,jj,kk,n
    total=0
    do kk=1,4
    do jj=1,4
    do ii=1,4
      total=total+state(:,ii,jj,kk)*g%volume_fraction(ii,jj,kk)*g%dx*g%dy*g%dz
    end do
    end do
    end do
    call liquid_mixture_properties(liquids,masses(:,1)/drops(1)%mass,drops(1)%temperature,rho,cp,e,ok)
    call require(ok,'liquid inventory')
    amount=drops(1)%mass*drops(1)%multiplicity
    total(1)=total(1)+amount;total(2:4)=total(2:4)+amount*drops(1)%velocity
    total(5)=total(5)+amount*(e+.5_dp*sum(drops(1)%velocity**2))
    do n=1,2
      total(5+liquids(n)%vapor_index)=total(5+liquids(n)%vapor_index)+masses(n,1)*drops(1)%multiplicity
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
