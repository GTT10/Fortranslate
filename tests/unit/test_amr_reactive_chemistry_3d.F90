program test_amr_reactive_chemistry_3d
  use precision_mod, only: dp
  use state_indices_mod, only: ncons
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use h2o2_elementary_mechanism_mod, only: &
    load_h2o2_elementary_mechanism
  use mesh_3d_mod, only: uniform_cell_centers_3d
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_entropy_wave_3d_problem_mod, only: &
    initialize_reactive_problem_3d
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,advance_reactive_boundary_hydro_3d
  use amr_hierarchy_3d_mod, only: &
    amr_patch_3d, initialize_amr_patch_3d, average_down_3d, &
    restrict_average_3d, composite_integrals_amr_3d
  use amr_reactive_3d_mod, only: &
    composite_element_integrals_amr_3d, element_integrals_from_reactive_integrals_3d, &
    advance_amr_reactive_chemistry_3d, &
    advance_amr_reactive_strang_3d, advance_amr_reactive_hydro_3d
  implicit none

  integer, parameter :: nx = 8, ny = 4, nz = 4, ratio = 2
  real(dp), parameter :: chemistry_interval = 1.0e-6_dp
  real(dp), parameter :: dx = 1.0e-3_dp
  real(dp), parameter :: dy = 1.0e-3_dp
  real(dp), parameter :: dz = 1.0e-3_dp
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(elementary_reaction), allocatable :: invalid_reactions(:)
  type(amr_patch_3d) :: patch
  real(dp), allocatable :: coarse_state(:, :, :, :)
  real(dp), allocatable :: fine_state(:, :, :, :)
  real(dp), allocatable :: coarse_temperature(:, :, :)
  real(dp), allocatable :: fine_temperature(:, :, :)
  real(dp), allocatable :: initial_coarse(:, :, :, :)
  real(dp), allocatable :: initial_fine(:, :, :, :)
  real(dp), allocatable :: initial_coarse_temperature(:, :, :)
  real(dp), allocatable :: initial_fine_temperature(:, :, :)
  real(dp) :: initial_elements(3), final_elements(3)
  real(dp) :: chemistry_change, element_error
  real(dp) :: reflux, reference_reflux
  logical :: ok

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "elementary thermodynamics load")
  call load_h2o2_elementary_mechanism(reactions, ok)
  call require(ok, "elementary mechanism load")
  call initialize_amr_patch_3d( &
    nx, ny, nz, 3, 6, 2, 3, 2, 3, ratio, patch, ok)
  call require(ok .and. patch%is_strictly_interior(), &
    "strictly interior chemistry patch")
  call initialize_uniform_hierarchy( &
    species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, ok)
  call require(ok, "uniform chemistry hierarchy initialization")
  initial_coarse = coarse_state
  initial_coarse_temperature = coarse_temperature
  initial_fine = fine_state
  initial_fine_temperature = fine_temperature

  call composite_element_integrals_amr_3d( &
    species, patch, coarse_state, fine_state, dx, dy, dz, &
    initial_elements, ok)
  call require(ok, "initial composite element integrals")
  call advance_amr_reactive_chemistry_3d( &
    species, reactions, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, chemistry_interval, &
    2.0e-7_dp, 1.0e-12_dp, ok)
  call require(ok, "transactional hierarchy chemistry source")
  chemistry_change = max(maxval(abs( &
    coarse_state(ncons + 1:, :, :, :) - &
    initial_coarse(ncons + 1:, :, :, :))), maxval(abs( &
    fine_state(ncons + 1:, :, :, :) - &
    initial_fine(ncons + 1:, :, :, :))))
  call require(chemistry_change > 1.0e-12_dp, "chemistry activity")
  call require(maxval(abs( &
    coarse_state(1:ncons, :, :, :) - &
    initial_coarse(1:ncons, :, :, :))) == 0.0_dp, &
    "coarse Euler invariants under chemistry")
  call require(maxval(abs( &
    fine_state(1:ncons, :, :, :) - &
    initial_fine(1:ncons, :, :, :))) == 0.0_dp, &
    "fine Euler invariants under chemistry")
  call require(levels_are_synchronized(coarse_state, fine_state, patch), &
    "chemistry average-down synchronization")
  call composite_element_integrals_amr_3d( &
    species, patch, coarse_state, fine_state, dx, dy, dz, &
    final_elements, ok)
  call require(ok, "final composite element integrals")
  element_error = maxval(abs(final_elements - initial_elements) / &
    max(1.0e-30_dp, abs(initial_elements)))
  call require(element_error <= 2.0e-10_dp, &
    "composite H/O/N conservation")

  invalid_reactions = reactions
  invalid_reactions(1)%reactant_stoich(1) = -1.0_dp
  call assert_source_rejection_is_transactional( &
    species, invalid_reactions, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, "explicit", ok)
  call require(ok, "invalid reaction rollback")
  call assert_source_rejection_is_transactional( &
    species, reactions, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, "invalid", ok)
  call require(ok, "invalid chemistry integrator rollback")

  coarse_state = initial_coarse
  coarse_temperature = initial_coarse_temperature
  fine_state = initial_fine
  fine_temperature = initial_fine_temperature
  call advance_amr_reactive_chemistry_3d( &
    species, reactions, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, 0.5_dp * chemistry_interval, &
    2.0e-7_dp, 1.0e-12_dp, ok)
  call require(ok, "first reference half chemistry phase")
  call advance_amr_reactive_chemistry_3d( &
    species, reactions, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, 0.5_dp * chemistry_interval, &
    2.0e-7_dp, 1.0e-12_dp, ok)
  call require(ok, "second reference half chemistry phase")

  call advance_amr_reactive_strang_3d( &
    species, reactions, patch, initial_coarse, initial_coarse_temperature, &
    initial_fine, initial_fine_temperature, dx, dy, dz, &
    chemistry_interval, "rusanov", .true., 2.0e-7_dp, 1.0e-12_dp, &
    reflux, ok, "characteristic_plm", "mc")
  call require(ok, "AMR reaction-hydro-reaction split")
  call require(maxval(abs(initial_coarse - coarse_state)) <= 5.0e-12_dp .and. &
    maxval(abs(initial_fine - fine_state)) <= 5.0e-12_dp .and. &
    maxval(abs(initial_coarse_temperature - coarse_temperature)) <= &
      5.0e-8_dp .and. &
    maxval(abs(initial_fine_temperature - fine_temperature)) <= 5.0e-8_dp, &
    "uniform split reduces to two source half phases")

  initial_coarse = coarse_state
  initial_coarse_temperature = coarse_temperature
  initial_fine = fine_state
  initial_fine_temperature = fine_temperature
  call advance_amr_reactive_strang_3d( &
    species, reactions, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, dx, dy, dz, chemistry_interval, &
    "invalid", .true., 2.0e-7_dp, 1.0e-12_dp, reflux, ok, &
    "characteristic_plm", "mc")
  call require(.not. ok .and. reflux == 0.0_dp .and. &
    hierarchy_matches(coarse_state, coarse_temperature, fine_state, &
      fine_temperature, initial_coarse, initial_coarse_temperature, &
      initial_fine, initial_fine_temperature), &
    "failed middle hydro rolls back the full split")

  call initialize_uniform_hierarchy( &
    species, patch, coarse_state, coarse_temperature, &
    fine_state, fine_temperature, ok)
  call require(ok, "disabled-path hierarchy reset")
  initial_coarse = coarse_state
  initial_coarse_temperature = coarse_temperature
  initial_fine = fine_state
  initial_fine_temperature = fine_temperature
  call advance_amr_reactive_hydro_3d( &
    species, patch, coarse_state, coarse_temperature, fine_state, &
    fine_temperature, dx, dy, dz, 1.0e-7_dp, "rusanov", &
    reference_reflux, ok, "characteristic_plm", "mc")
  call require(ok, "disabled-path hydro reference")
  call advance_amr_reactive_strang_3d( &
    species, reactions, patch, initial_coarse, initial_coarse_temperature, &
    initial_fine, initial_fine_temperature, dx, dy, dz, 1.0e-7_dp, &
    "rusanov", .false., 2.0e-7_dp, 1.0e-12_dp, reflux, ok, &
    "characteristic_plm", "mc")
  call require(ok .and. reflux == reference_reflux .and. &
    hierarchy_matches(initial_coarse, initial_coarse_temperature, &
      initial_fine, initial_fine_temperature, coarse_state, &
      coarse_temperature, fine_state, fine_temperature), &
    "chemistry-disabled exact hydro dispatch")

  call check_physical_split()
  call check_boundary_touching_hydro()

  write(*, '(2(a,es24.16))') &
    "chemistry change=", chemistry_change, ", element error=", element_error
  write(*, '(a)') "test_amr_reactive_chemistry_3d: PASS"

contains

  subroutine check_boundary_touching_hydro()
    type(amr_patch_3d) :: p
    real(dp), allocatable :: c0(:,:,:,:),t0(:,:,:),f0(:,:,:,:),ft0(:,:,:)
    real(dp), allocatable :: c(:,:,:,:),t(:,:,:),f(:,:,:,:),ft(:,:,:),ref(:,:,:,:),rt(:,:,:)
    real(dp), allocatable :: inlet(:,:),it(:),before(:),after(:),ledger(:),rl(:),total_ledger(:),scale(:)
    real(dp) :: r,interval,velocity(3),factor,inlet_density
    integer :: mode,scheme,lo(3),hi(3),a,i,j,k,nv,rr,substep
    character(len=32) :: faces(6),method,bad_faces(6)
    interval=1e-8_dp;velocity=[30._dp,-20._dp,10._dp]
    do mode=1,8
      lo=2;hi=3;rr=2
      if (mode<=6) then
        a=(mode+1)/2
        if (modulo(mode,2)==1) lo(a)=1
        if (modulo(mode,2)==0) hi(a)=4
      else if (mode==7) then
        lo=1;hi=4
      else
        lo=1;hi=2;rr=3
      end if
      call initialize_amr_patch_3d(4,4,4,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3),rr,p,ok)
      call require(ok,'boundary touching patch')
      call initialize_uniform_hierarchy(species,p,c0,t0,f0,ft0,ok)
      call require(ok,'boundary touching initial state')
      nv=size(c0,1)
      allocate(inlet(nv,6),it(6),before(nv),after(nv),ledger(nv),rl(nv),total_ledger(nv),scale(nv))
      do k=1,4
        do j=1,4
          do i=1,4
            c0(2:4,i,j,k)=c0(1,i,j,k)*velocity
            c0(5,i,j,k)=c0(5,i,j,k)+.5_dp*c0(1,i,j,k)*sum(velocity**2)
          end do
        end do
      end do
      inlet=spread(c0(:,1,1,1),2,6);it=t0(1,1,1)
      inlet(:,1)=1.1_dp*inlet(:,1)
      inlet_density=inlet(1,1)
      do k=1,size(f0,4)
        do j=1,size(f0,3)
          do i=1,size(f0,2)
            factor=1+.05_dp*sin(real(i+2*j+3*k,dp))
            f0(:,i,j,k)=factor*f0(:,i,j,k)
            f0(2:4,i,j,k)=f0(1,i,j,k)*velocity
            f0(5,i,j,k)=f0(5,i,j,k)+.5_dp*f0(1,i,j,k)*sum(velocity**2)
          end do
        end do
      end do
      call average_down_3d(c0,f0,p,ok)
      call require(ok,'boundary touching synchronization')
      faces=['inflow                          ','outflow                         ', &
        'slip_wall                       ','no_slip_wall                    ', &
        'outflow                         ','slip_wall                       ']
      call composite_integrals_amr_3d(c0,f0,p,dx,dy,dz,before,ok)
      call require(ok,'boundary touching initial inventory')
      do scheme=1,2
        method='pcm'
        if (scheme==2) method='characteristic_plm'
        c=c0;t=t0;f=f0;ft=ft0
        call advance_amr_reactive_hydro_3d(species,p,c,t,f,ft,dx,dy,dz,interval,'rusanov',r,ok, &
          method,'mc',faces,ledger,inlet,it)
        if (.not.ok) write(*,*) 'boundary hydro case:',mode,trim(method)
        call require(ok,'boundary touching hydro accepted')
        call composite_integrals_amr_3d(c,f,p,dx,dy,dz,after,ok)
        scale=max(abs(before),abs(ledger),1e-20_dp);scale(2:4)=maxval(scale(2:4))
        call require(ok.and.all(abs(after-before-ledger)<=2e-11_dp*scale),'boundary touching conserved balance')
        if (mode==7) then
          ref=f0;rt=ft0;total_ledger=0
          do substep=1,rr
            call advance_reactive_boundary_hydro_3d(species,ref,rt,[dx,dy,dz]/rr,interval/rr, &
              faces,method,'mc','rusanov',rl,ok,inlet,it)
            call require(ok,'full refinement uniform-grid reference')
            total_ledger=total_ledger+rl
          end do
          call require(all(ref==f).and.all(rt==ft).and.r==0,'full refinement exact fine-grid equality')
          call require(all(abs(total_ledger-ledger)<=2e-11_dp*scale),'full refinement boundary transfer')
        end if
        c=c0;t=t0;f=f0;ft=ft0;inlet(1,1)=-1
        call advance_amr_reactive_hydro_3d(species,p,c,t,f,ft,dx,dy,dz,interval,'rusanov',r,ok, &
          method,'mc',faces,ledger,inlet,it)
        call require(.not.ok.and.hierarchy_matches(c,t,f,ft,c0,t0,f0,ft0).and.r==0.and.all(ledger==0), &
          'boundary touching invalid inlet rollback')
        inlet(1,1)=inlet_density
        bad_faces='periodic'
        call advance_amr_reactive_hydro_3d(species,p,c,t,f,ft,dx,dy,dz,interval,'rusanov',r,ok, &
          method,'mc',bad_faces,ledger)
        call require(.not.ok.and.hierarchy_matches(c,t,f,ft,c0,t0,f0,ft0).and.r==0.and.all(ledger==0), &
          'boundary touching periodic seam rejected')
      end do
      deallocate(inlet,it,before,after,ledger,rl,total_ledger,scale)
    end do
  end subroutine check_boundary_touching_hydro

  subroutine check_physical_split()
    real(dp), allocatable :: c(:,:,:,:),ct(:,:,:),f(:,:,:,:),ft(:,:,:)
    real(dp), allocatable :: c0(:,:,:,:),ct0(:,:,:),f0(:,:,:,:),ft0(:,:,:)
    real(dp), allocatable :: cr(:,:,:,:),ctr(:,:,:),fr(:,:,:,:),ftr(:,:,:)
    real(dp), allocatable :: inlet(:,:),it(:),ledger(:),reference_ledger(:),before(:),after(:)
    real(dp) :: e0(3),e1(3),eb(3),r,rr,scale(5)
    character(len=32) :: faces(6),method
    integer :: mode,nv
    call initialize_uniform_hierarchy(species,patch,c0,ct0,f0,ft0,ok)
    call require(ok,'physical split initial hierarchy')
    nv=size(c0,1)
    allocate(inlet(nv,6),it(6),ledger(nv),reference_ledger(nv),before(nv),after(nv))
    inlet=spread(c0(:,1,1,1),2,6);it=ct0(1,1,1)
    faces='periodic';faces(1)='inflow';faces(2)='outflow'
    ! Different reservoir density makes the boundary budget nonzero.
    inlet(:,1)=1.1_dp*inlet(:,1)
    call composite_integrals_amr_3d(c0,f0,patch,dx,dy,dz,before,ok)
    call require(ok,'physical split initial conserved inventory')
    call element_integrals_from_reactive_integrals_3d(species,before,e0,ok)
    call require(ok,'physical split initial elements')
    do mode=1,2
      method='pcm'
      if (mode==2) method='characteristic_plm'
      c=c0;ct=ct0;f=f0;ft=ft0
      call advance_amr_reactive_strang_3d(species,reactions,patch,c,ct,f,ft,dx,dy,dz, &
        chemistry_interval,'rusanov',.true.,2e-7_dp,1e-12_dp,r,ok,method,'mc', &
        boundary_faces=faces,boundary_change=ledger,inflow_state=inlet,inflow_temperature=it)
      call require(ok,'physical reaction-hydro-reaction split')
      call composite_integrals_amr_3d(c,f,patch,dx,dy,dz,after,ok)
      scale=max(abs(before(1:5)),abs(ledger(1:5)),1e-20_dp)
      scale(2:4)=maxval(scale(2:4))
      call require(ok.and.all(abs(after(1:5)-before(1:5)-ledger(1:5)) &
        <=2e-11_dp*scale),'physical split Euler boundary balance')
      call require(abs(ledger(1))>1e-20_dp,'physical split nonzero mass boundary flux')
      call element_integrals_from_reactive_integrals_3d(species,after,e1,ok)
      call require(ok,'physical split final elements')
      call element_integrals_from_reactive_integrals_3d(species,ledger,eb,ok)
      call require(ok.and.all(abs(e1-e0-eb)<=2e-10_dp*max(abs(e0),1e-30_dp)), &
        'physical split elemental boundary balance')
      cr=c0;ctr=ct0;fr=f0;ftr=ft0
      call advance_amr_reactive_chemistry_3d(species,reactions,patch,cr,ctr,fr,ftr, &
        .5_dp*chemistry_interval,2e-7_dp,1e-12_dp,ok)
      call require(ok,'physical reference first reaction')
      call advance_amr_reactive_hydro_3d(species,patch,cr,ctr,fr,ftr,dx,dy,dz,chemistry_interval, &
        'rusanov',rr,ok,method,'mc',faces,reference_ledger,inlet,it)
      call require(ok,'physical reference hydro')
      call advance_amr_reactive_chemistry_3d(species,reactions,patch,cr,ctr,fr,ftr, &
        .5_dp*chemistry_interval,2e-7_dp,1e-12_dp,ok)
      call require(ok.and.hierarchy_matches(c,ct,f,ft,cr,ctr,fr,ftr).and.r==rr.and. &
        all(ledger==reference_ledger),'physical split matches explicit phases exactly')
      call require(maxval(abs(after(6:)-before(6:)-ledger(6:)))>1e-20_dp, &
        'physical split exercises chemical conversion')
      c=c0;ct=ct0;f=f0;ft=ft0
      inlet(1,1)=-1
      call advance_amr_reactive_strang_3d(species,reactions,patch,c,ct,f,ft,dx,dy,dz, &
        chemistry_interval,'rusanov',.true.,2e-7_dp,1e-12_dp,r,ok,method,'mc', &
        boundary_faces=faces,boundary_change=ledger,inflow_state=inlet,inflow_temperature=it)
      call require(.not.ok.and.hierarchy_matches(c,ct,f,ft,c0,ct0,f0,ft0).and.r==0.and.all(ledger==0), &
        'invalid physical inlet rolls back first reaction and all outputs')
      inlet(1,1)=1.1_dp*c0(1,1,1,1)
      cr=c0;ctr=ct0;fr=f0;ftr=ft0
      call advance_amr_reactive_hydro_3d(species,patch,cr,ctr,fr,ftr,dx,dy,dz,chemistry_interval, &
        'rusanov',rr,ok,method,'mc',faces,reference_ledger,inlet,it)
      call require(ok,'physical disabled chemistry reference')
      call advance_amr_reactive_strang_3d(species,reactions,patch,c,ct,f,ft,dx,dy,dz, &
        chemistry_interval,'rusanov',.false.,2e-7_dp,1e-12_dp,r,ok,method,'mc', &
        boundary_faces=faces,boundary_change=ledger,inflow_state=inlet,inflow_temperature=it)
      call require(ok.and.hierarchy_matches(c,ct,f,ft,cr,ctr,fr,ftr).and.r==rr.and. &
        all(ledger==reference_ledger),'physical chemistry-disabled exact dispatch')
    end do
  end subroutine check_physical_split

  subroutine initialize_uniform_hierarchy( &
      species, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), allocatable, intent(out) :: coarse_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: coarse_temperature(:, :, :)
    real(dp), allocatable, intent(out) :: fine_state(:, :, :, :)
    real(dp), allocatable, intent(out) :: fine_temperature(:, :, :)
    logical, intent(out) :: ok

    type(reactive_3d_config) :: config, fine_config
    real(dp), allocatable :: x(:), y(:), z(:), xf(:), yf(:), zf(:)
    real(dp), allocatable :: mass_fractions(:), recovered(:, :, :)
    real(dp) :: local_dx, local_dy, local_dz, density, fine_density
    logical :: local_ok
    integer :: nvar

    ok = .false.
    config = reactive_3d_config()
    config%nx = patch%coarse_nx
    config%ny = patch%coarse_ny
    config%nz = patch%coarse_nz
    config%x_upper = real(config%nx, dp) * dx
    config%y_upper = real(config%ny, dp) * dy
    config%z_upper = real(config%nz, dp) * dz
    config%problem = "uniform_reactor"
    config%thermo_model = "elementary"
    config%initial_temperature = 1200.0_dp
    config%initial_velocity_x = 0.0_dp
    config%initial_velocity_y = 0.0_dp
    config%initial_velocity_z = 0.0_dp
    nvar = reactive_nvar(size(species))
    allocate(coarse_state(nvar, config%nx, config%ny, config%nz))
    allocate(coarse_temperature(config%nx, config%ny, config%nz))
    allocate(fine_state(nvar, patch%fine_nx(), patch%fine_ny(), &
      patch%fine_nz()))
    allocate(fine_temperature( &
      patch%fine_nx(), patch%fine_ny(), patch%fine_nz()))
    allocate(x(config%nx), y(config%ny), z(config%nz))
    allocate(xf(patch%fine_nx()), yf(patch%fine_ny()), zf(patch%fine_nz()))
    allocate(mass_fractions(size(species)))
    allocate(recovered(config%nx, config%ny, config%nz))
    call uniform_cell_centers_3d( &
      config%nx, config%ny, config%nz, config%x_lower, config%x_upper, &
      config%y_lower, config%y_upper, config%z_lower, config%z_upper, &
      x, y, z, local_dx, local_dy, local_dz)
    call initialize_reactive_problem_3d( &
      species, config, x, y, z, coarse_state, coarse_temperature, &
      density, mass_fractions, local_ok)
    if (.not. local_ok) return
    fine_config = config
    fine_config%nx = patch%fine_nx()
    fine_config%ny = patch%fine_ny()
    fine_config%nz = patch%fine_nz()
    fine_config%x_lower = config%x_lower + &
      real(patch%coarse_i_lower - 1, dp) * dx
    fine_config%x_upper = config%x_lower + &
      real(patch%coarse_i_upper, dp) * dx
    fine_config%y_lower = config%y_lower + &
      real(patch%coarse_j_lower - 1, dp) * dy
    fine_config%y_upper = config%y_lower + &
      real(patch%coarse_j_upper, dp) * dy
    fine_config%z_lower = config%z_lower + &
      real(patch%coarse_k_lower - 1, dp) * dz
    fine_config%z_upper = config%z_lower + &
      real(patch%coarse_k_upper, dp) * dz
    call uniform_cell_centers_3d( &
      fine_config%nx, fine_config%ny, fine_config%nz, &
      fine_config%x_lower, fine_config%x_upper, fine_config%y_lower, &
      fine_config%y_upper, fine_config%z_lower, fine_config%z_upper, &
      xf, yf, zf, local_dx, local_dy, local_dz)
    call initialize_reactive_problem_3d( &
      species, fine_config, xf, yf, zf, fine_state, fine_temperature, &
      fine_density, mass_fractions, local_ok)
    if (.not. local_ok) return
    call average_down_3d(coarse_state, fine_state, patch, local_ok)
    if (.not. local_ok) return
    call recover_reactive_temperatures_3d( &
      species, coarse_state, coarse_temperature, config%nx, config%ny, &
      config%nz, recovered, local_ok)
    if (.not. local_ok) return
    coarse_temperature = recovered
    ok = .true.
  end subroutine initialize_uniform_hierarchy

  subroutine assert_source_rejection_is_transactional( &
      species, reactions, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, integrator, ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_patch_3d), intent(in) :: patch
    real(dp), intent(inout) :: coarse_state(:, :, :, :)
    real(dp), intent(inout) :: coarse_temperature(:, :, :)
    real(dp), intent(inout) :: fine_state(:, :, :, :)
    real(dp), intent(inout) :: fine_temperature(:, :, :)
    character(len=*), intent(in) :: integrator
    logical, intent(out) :: ok

    real(dp), allocatable :: saved_coarse(:, :, :, :)
    real(dp), allocatable :: saved_fine(:, :, :, :)
    real(dp), allocatable :: saved_coarse_temperature(:, :, :)
    real(dp), allocatable :: saved_fine_temperature(:, :, :)
    logical :: local_ok

    saved_coarse = coarse_state
    saved_coarse_temperature = coarse_temperature
    saved_fine = fine_state
    saved_fine_temperature = fine_temperature
    call advance_amr_reactive_chemistry_3d( &
      species, reactions, patch, coarse_state, coarse_temperature, &
      fine_state, fine_temperature, chemistry_interval, &
      2.0e-7_dp, 1.0e-12_dp, local_ok, integrator)
    ok = .not. local_ok .and. hierarchy_matches( &
      coarse_state, coarse_temperature, fine_state, fine_temperature, &
      saved_coarse, saved_coarse_temperature, saved_fine, &
      saved_fine_temperature)
  end subroutine assert_source_rejection_is_transactional

  logical function levels_are_synchronized(coarse, fine, patch) result(match)
    real(dp), intent(in) :: coarse(:, :, :, :), fine(:, :, :, :)
    type(amr_patch_3d), intent(in) :: patch

    real(dp), allocatable :: restricted(:, :, :, :)
    logical :: local_ok

    allocate(restricted(size(coarse, 1), &
      patch%coarse_i_upper - patch%coarse_i_lower + 1, &
      patch%coarse_j_upper - patch%coarse_j_lower + 1, &
      patch%coarse_k_upper - patch%coarse_k_lower + 1))
    call restrict_average_3d(fine, patch, restricted, local_ok)
    match = local_ok
    if (match) match = maxval(abs(restricted - coarse(:, &
      patch%coarse_i_lower:patch%coarse_i_upper, &
      patch%coarse_j_lower:patch%coarse_j_upper, &
      patch%coarse_k_lower:patch%coarse_k_upper))) == 0.0_dp
  end function levels_are_synchronized

  logical function hierarchy_matches( &
      coarse_left, coarse_temperature_left, fine_left, &
      fine_temperature_left, coarse_right, coarse_temperature_right, &
      fine_right, fine_temperature_right) result(match)
    real(dp), intent(in) :: coarse_left(:, :, :, :), coarse_right(:, :, :, :)
    real(dp), intent(in) :: fine_left(:, :, :, :), fine_right(:, :, :, :)
    real(dp), intent(in) :: coarse_temperature_left(:, :, :)
    real(dp), intent(in) :: coarse_temperature_right(:, :, :)
    real(dp), intent(in) :: fine_temperature_left(:, :, :)
    real(dp), intent(in) :: fine_temperature_right(:, :, :)

    match = all(shape(coarse_left) == shape(coarse_right)) .and. &
      all(shape(fine_left) == shape(fine_right)) .and. &
      all(shape(coarse_temperature_left) == &
        shape(coarse_temperature_right)) .and. &
      all(shape(fine_temperature_left) == shape(fine_temperature_right))
    if (.not. match) return
    match = maxval(abs(coarse_left - coarse_right)) == 0.0_dp .and. &
      maxval(abs(fine_left - fine_right)) == 0.0_dp .and. &
      maxval(abs(coarse_temperature_left - &
        coarse_temperature_right)) == 0.0_dp .and. &
      maxval(abs(fine_temperature_left - fine_temperature_right)) == 0.0_dp
  end function hierarchy_matches

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') "FAILED:", trim(label)
      error stop 1
    end if
  end subroutine require

end program test_amr_reactive_chemistry_3d
