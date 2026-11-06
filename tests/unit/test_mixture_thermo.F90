program test_mixture_thermo
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_is_finite, ieee_positive_inf, ieee_quiet_nan, ieee_value
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant, &
    configure_nasa7_eos_policy,nasa7_eos_bounds,nasa7_eos_policy_identity
  use thermo_database_mod, only: &
    load_gri30_thermo_subset, gri_h2_index, gri_o2_index, &
    gri_h2o_index, gri_n2_index
  use mixture_thermo_mod, only: &
    valid_mixture_composition, mixture_mass_properties, mixture_pressure, &
    mixture_density, mixture_sound_speed, temperature_from_internal_energy, &
    mixture_temperature_bounds, mixture_molecular_weight, &
    mixture_specific_gas_constant, mass_fractions_from_mole_fractions, &
    mole_fractions_from_mass_fractions
  implicit none

  type(nasa7_species), allocatable :: species(:)
  real(dp) :: mass_fractions(4), invalid_mass_fractions(4)
  real(dp) :: molecular_weight, gas_constant, cp, cv, gamma
  real(dp) :: enthalpy, internal_energy, entropy, pressure, density
  real(dp) :: sound_speed, recovered_temperature, target_temperature
  logical :: ok
  integer :: iterations, i
  real(dp), parameter :: temperatures(3) = [300.0_dp, 1200.0_dp, 2500.0_dp]

  call load_gri30_thermo_subset(species, ok)
  if (.not. ok) error stop "Failed to load mixture species"

  mass_fractions = 0.0_dp
  mass_fractions(gri_o2_index) = 0.23291751145757963_dp
  mass_fractions(gri_n2_index) = 0.7670824885424203_dp
  if (.not. valid_mixture_composition(species, mass_fractions)) then
    error stop "Valid air composition was rejected"
  end if

  call mixture_mass_properties( &
    species, mass_fractions, 1200.0_dp, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
  if (.not. ok) error stop "Mixture property evaluation failed"

  call assert_close( &
    molecular_weight, 28.85067068107281_dp, 2.0e-12_dp, "Wmix")
  call assert_close( &
    gas_constant, 288.18957833128849_dp, 2.0e-12_dp, "Rmix")
  call assert_close(cp, 1182.4886092588142_dp, 2.0e-12_dp, "cp")
  call assert_close(cv, 894.29903092752579_dp, 2.0e-12_dp, "cv")
  call assert_close(gamma, 1.3222519183906434_dp, 2.0e-12_dp, "gamma")
  call assert_close( &
    enthalpy, 986630.08716612426_dp, 2.0e-12_dp, "enthalpy")
  call assert_close( &
    internal_energy, 640802.59316857834_dp, 2.0e-12_dp, &
    "internal energy")

  pressure = mixture_pressure( &
    species, mass_fractions, 1.2_dp, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture pressure evaluation failed"
  call assert_close( &
    pressure, 414992.99279705540_dp, 2.0e-12_dp, "pressure")

  density = mixture_density( &
    species, mass_fractions, pressure, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture density evaluation failed"
  call assert_close(density, 1.2_dp, 2.0e-13_dp, "density round trip")

  sound_speed = mixture_sound_speed( &
    species, mass_fractions, 1200.0_dp, ok)
  if (.not. ok) error stop "Mixture sound-speed evaluation failed"
  call assert_close( &
    sound_speed, 676.21820987790932_dp, 2.0e-12_dp, "sound speed")

  do i = 1, size(temperatures)
    target_temperature = temperatures(i)
    call mixture_mass_properties( &
      species, mass_fractions, target_temperature, molecular_weight, &
      gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
    if (.not. ok) error stop "Mixture energy evaluation failed"
    call temperature_from_internal_energy( &
      species, mass_fractions, internal_energy, 900.0_dp, &
      recovered_temperature, ok, iterations)
    if (.not. ok) error stop "e-to-T inversion failed"
    call assert_close( &
      recovered_temperature, target_temperature, 2.0e-11_dp, &
      "e-to-T round trip")
    if (iterations < 1 .or. iterations > 100) then
      error stop "Invalid temperature inversion iteration count"
    end if
  end do

  call mixture_mass_properties( &
    species, mass_fractions, 300.0_dp, molecular_weight, gas_constant, &
    cp, cv, gamma, enthalpy, internal_energy, entropy, ok)
  call temperature_from_internal_energy( &
    species, mass_fractions, internal_energy - 1.0e6_dp, 300.0_dp, &
    recovered_temperature, ok)
  if (ok) error stop "Out-of-range internal energy was accepted"

  invalid_mass_fractions = mass_fractions
  invalid_mass_fractions(gri_h2_index) = 0.01_dp
  if (valid_mixture_composition(species, invalid_mass_fractions)) then
    error stop "Non-unit mixture composition was accepted"
  end if

  invalid_mass_fractions = mass_fractions
  invalid_mass_fractions(gri_h2o_index) = -1.0e-4_dp
  invalid_mass_fractions(gri_n2_index) = &
    invalid_mass_fractions(gri_n2_index) + 1.0e-4_dp
  if (valid_mixture_composition(species, invalid_mass_fractions)) then
    error stop "Negative mixture composition was accepted"
  end if

  call check_nonfinite_contract(species, mass_fractions, ok)
  if (.not. ok) error stop "Mixture thermodynamics finite contract failed"

  block
    real(dp) :: near_zero(4),target,below,above,lower_energy,upper_energy
    integer :: sample
    near_zero=[0._dp,1.0351682798876758_dp,-1._dp,1e-6_dp]
    mass_fractions=0;mass_fractions(gri_n2_index)=1
    do sample=1,size(near_zero)
      target=near_zero(sample)
      call temperature_from_internal_energy(species,mass_fractions,target,440._dp,recovered_temperature,ok,iterations)
      if(.not.ok.or.iterations>=100) error stop 'Near-zero energy inversion stagnates'
      below=nearest(recovered_temperature,-1._dp);above=nearest(recovered_temperature,1._dp)
      call mixture_mass_properties(species,mass_fractions,below,molecular_weight,gas_constant, &
        cp,cv,gamma,enthalpy,lower_energy,entropy,ok)
      if(.not.ok) error stop 'Lower neighbor thermodynamics'
      call mixture_mass_properties(species,mass_fractions,above,molecular_weight,gas_constant, &
        cp,cv,gamma,enthalpy,upper_energy,entropy,ok)
      if(.not.ok) error stop 'Upper neighbor thermodynamics'
      if(target<lower_energy.or.target>upper_energy) error stop 'Recovered temperature does not bracket target energy'
    end do
  end block

  block
    type(nasa7_species), allocatable :: full_species(:)
    real(dp), allocatable :: pure_n2(:)
    real(dp) :: target,guesses(2)
    integer :: n,guess
    call load_h2o2_full_thermo(full_species,ok)
    if(.not.ok) error stop 'Midpoint inversion species load'
    allocate(pure_n2(size(full_species)));pure_n2=0
    do n=1,size(full_species)
      if(full_species(n)%name=='N2') pure_n2(n)=1
    enddo
    call mixture_mass_properties(full_species,pure_n2,1000._dp,molecular_weight,gas_constant, &
      cp,cv,gamma,enthalpy,target,entropy,ok)
    if(.not.ok) error stop 'Midpoint energy evaluation'
    guesses=[1000._dp,nearest(1000._dp,1._dp)]
    do guess=1,size(guesses)
      call temperature_from_internal_energy(full_species,pure_n2,target,guesses(guess), &
        recovered_temperature,ok,iterations)
      if(.not.ok.or.abs(recovered_temperature-1000._dp)>1e-9_dp) &
        error stop 'Midpoint high-fit inversion depends on exact or upper-neighbor warm start'
    enddo
  end block

  call check_pinned_extrapolation()
  call check_pinned_convergence()

  write(*, '(a)') "test_mixture_thermo: PASS"

contains

  subroutine check_pinned_convergence()
    type(nasa7_species) :: table(2)
    character(len=128) :: default_identity,converged_identity
    real(dp) :: guess,recovered
    logical :: good
    table%name='constant';table%molecular_weight=universal_gas_constant
    table%temperature_min=300;table%temperature_mid=1000;table%temperature_max=3500
    table(1)%low_coefficients=[2._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp]
    table(1)%high_coefficients=table(1)%low_coefficients;table(2)=table(1)
    call configure_nasa7_eos_policy(table,.true.,good)
    if(.not.good) error stop 'Default caloric policy setup'
    default_identity=nasa7_eos_policy_identity(table)
    guess=2000._dp+5e-7_dp
    ! cv=1 and e=T: the native default deliberately retains this warm start.
    call temperature_from_internal_energy(table,[1._dp,0._dp],2000._dp,guess,recovered,good)
    if(.not.good.or.recovered/=guess) error stop 'Default caloric stopping rule changed'
    call configure_nasa7_eos_policy(table,.true.,good,converged_inverse=.true.)
    if(.not.good) error stop 'Converged caloric policy setup'
    converged_identity=nasa7_eos_policy_identity(table)
    if(converged_identity==default_identity.or.index(converged_identity,'PELEF_FUEGO_EOS_CONVERGENCE_V1')/=1) &
      error stop 'Caloric precision missing from restart identity'
    call temperature_from_internal_energy(table,[1._dp,0._dp],2000._dp,guess,recovered,good)
    if(.not.good.or.abs(recovered-2000)>1e-12_dp) error stop 'Converged caloric analytic solution'
    call configure_nasa7_eos_policy(table,.false.,good,converged_inverse=.true.)
    if(good.or.nasa7_eos_policy_identity(table)/=converged_identity) error stop 'Invalid precision setup changed policy'
    table(1)%pinned_eos_convergence=.false.
    if(valid_mixture_composition(table,[1._dp,0._dp])) error stop 'Mixed caloric precision accepted'
    if(nasa7_eos_policy_identity(table)/='PELEF_INVALID_EOS_POLICY') error stop 'Mixed caloric identity accepted'
    call configure_nasa7_eos_policy(table,.true.,good)
    if(.not.good.or.nasa7_eos_policy_identity(table)/=default_identity) error stop 'Default caloric identity not restored'
  end subroutine check_pinned_convergence

  subroutine check_pinned_extrapolation()
    type(nasa7_species), allocatable :: table(:)
    type(nasa7_species) :: constant(1)
    real(dp), allocatable :: y(:)
    real(dp) :: w,r,cp,cv,g,h,e,s,target,recovered,cs,p,boundaries(2),offsets(2)
    integer :: i,j,count
    logical :: good
    call load_h2o2_full_thermo(table,good)
    if(.not.good) error stop 'Pinned EOS table load'
    allocate(y(size(table)));y=0
    do i=1,size(table)
      if(table(i)%name=='N2') y(i)=1
    enddo
    boundaries=[300._dp,3500._dp];offsets=[-10._dp,25._dp]
    do j=1,2
      call mixture_mass_properties(table,y,boundaries(j),w,r,cp,cv,g,h,e,s,good)
      if(.not.good) error stop 'Pinned boundary properties'
      target=e+cv*offsets(j)
      call temperature_from_internal_energy(table,y,target,1000._dp,recovered,good)
      if(good.or.recovered/=0) error stop 'Strict mode still rejects boundary energy outside fits'
      table%pinned_eos_extrapolation=.true.
      call temperature_from_internal_energy(table,y,target,1000._dp,recovered,good,count)
      if(.not.good.or.count/=0) error stop 'Pinned boundary extrapolation refused'
      call assert_close(recovered,boundaries(j)+offsets(j),2e-13_dp,'pinned extrapolated temperature')
      p=mixture_pressure(table,y,1._dp,recovered,good)
      if(.not.good.or.p<=0) error stop 'Extrapolated pressure'
      cs=mixture_sound_speed(table,y,recovered,good)
      if(.not.good.or.cs<=0) error stop 'Extrapolated sound speed'
      table%pinned_eos_extrapolation=.false.
    enddo
    table%pinned_eos_extrapolation=.true.
    call mixture_mass_properties(table,y,300._dp,w,r,cp,cv,g,h,e,s,good)
    target=e-301._dp*cv
    call temperature_from_internal_energy(table,y,target,1000._dp,recovered,good)
    if(good.or.recovered/=0) error stop 'Negative extrapolated temperature refused'
    ! R=1, cp/R=2 gives exactly cv=1 and e(300)=300: no rounded
    ! cancellation can turn the intended zero into a small positive result.
    constant=table(:1);constant(1)%molecular_weight=universal_gas_constant
    constant(1)%low_coefficients=[2._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp]
    constant(1)%high_coefficients=constant(1)%low_coefficients
    call mixture_mass_properties(constant,[1._dp],300._dp,w,r,cp,cv,g,h,target,s,good)
    if(.not.good.or.cv/=1.or.target/=300) error stop 'Exact zero-temperature fixture'
    call temperature_from_internal_energy(constant,[1._dp],0._dp,1000._dp,recovered,good)
    if(good.or.recovered/=0) error stop 'Zero extrapolated temperature refused'
    target=ieee_value(1._dp,ieee_quiet_nan)
    call temperature_from_internal_energy(table,y,target,1000._dp,recovered,good)
    if(good.or.recovered/=0) error stop 'Nonfinite extrapolated energy refused'
    call temperature_from_internal_energy(table,y,huge(1._dp)/4,1000._dp,recovered,good)
    if(good.or.recovered/=0) error stop 'Overflowing extrapolated properties refused'
    table(1)%pinned_eos_extrapolation=.false.
    if(valid_mixture_composition(table,y)) error stop 'Mixed EOS policies refused'
    ! Generated headers use mechanism-specific common fit limits, with a
    ! few explicit exceptions. Both selections retain coefficient metadata.
    table%pinned_eos_extrapolation=.true.;table%temperature_min=200;table%temperature_max=5000
    call mixture_mass_properties(table,y,200._dp,w,r,cp,cv,g,h,e,s,good)
    target=e-10._dp*cv
    call temperature_from_internal_energy(table,y,target,1000._dp,recovered,good)
    if(.not.good) error stop 'Mechanism-specific inversion bounds refused'
    call assert_close(recovered,190._dp,2e-13_dp,'mechanism-specific lower inversion boundary')
    call configure_nasa7_eos_policy(table,.true.,good,300._dp,3500._dp)
    if(.not.good) error stop 'Explicit inversion boundary setup'
    if(any(table%temperature_min/=200).or.any(table%temperature_max/=5000)) &
      error stop 'EOS policy changed fitted metadata'
    call nasa7_eos_bounds(table,target,recovered,good)
    if(.not.good.or.target/=300.or.recovered/=3500) error stop 'Explicit inversion bounds'
    call mixture_mass_properties(table,y,300._dp,w,r,cp,cv,g,h,e,s,good)
    target=e-10._dp*cv
    call temperature_from_internal_energy(table,y,target,1000._dp,recovered,good)
    if(.not.good) error stop 'Explicit reference inversion bounds refused'
    call assert_close(recovered,290._dp,2e-13_dp,'explicit reference inversion boundary')
    call configure_nasa7_eos_policy(table,.true.,good,3500._dp,300._dp)
    if(good.or.any(table%pinned_eos_temperature_min/=300).or.any(table%pinned_eos_temperature_max/=3500)) &
      error stop 'Invalid boundary setup was not transactional'
    call configure_nasa7_eos_policy(table,.false.,good,300._dp,3500._dp)
    if(good.or..not.all(table%pinned_eos_extrapolation)) error stop 'Inactive nonzero bounds accepted'
    table(1)%pinned_eos_temperature_max=4000
    if(valid_mixture_composition(table,y)) error stop 'Mixed inversion bounds accepted'
    table(1)%pinned_eos_temperature_max=3500
    if(index(nasa7_eos_policy_identity(table),'PELEF_FUEGO_EOS_V2')/=1) error stop 'EOS policy identity'
    ! PelePhysics 6e8026fc, Mechanisms/LiDryer/mechanism.H, species 8.
    ! Exact header coefficients; the existing Cantera table is a different dataset.
    table=table(:1);y=[1._dp]
    table(1)%name='N2';table(1)%molecular_weight=28.014_dp
    table(1)%temperature_min=300;table(1)%temperature_mid=1000;table(1)%temperature_max=3500
    table(1)%low_coefficients=[3.298677_dp,1.408240e-3_dp,-3.963222e-6_dp, &
      5.641515e-9_dp,-2.444855e-12_dp,-1020.900_dp,3.950372_dp]
    table(1)%high_coefficients=[2.926640_dp,1.487977e-3_dp,-5.684761e-7_dp, &
      1.009704e-10_dp,-6.753351e-15_dp,-922.7977_dp,5.980528_dp]
    ! Energies produced by compiling the byte-verified, unmodified C++ header.
    call temperature_from_internal_energy(table,y,-94478.993881851493_dp,700._dp,recovered,good)
    if(.not.good) error stop 'Pinned C++ low energy refused'
    call assert_close(recovered,290._dp,2e-12_dp,'pinned C++ low extrapolation')
    call temperature_from_internal_energy(table,y,2960713.5535413623_dp,700._dp,recovered,good)
    if(.not.good) error stop 'Pinned C++ high energy refused'
    call assert_close(recovered,3525._dp,2e-12_dp,'pinned C++ high extrapolation')
    call temperature_from_internal_energy(table,y,469600.84168867877_dp,700._dp,recovered,good)
    if(.not.good) error stop 'Pinned C++ midpoint energy refused'
    call assert_close(recovered,1000._dp,2e-12_dp,'pinned C++ midpoint inversion')
  end subroutine check_pinned_extrapolation


  subroutine check_nonfinite_contract(active_species, valid_mass_fractions, ok_out)
    type(nasa7_species), intent(in) :: active_species(:)
    real(dp), intent(in) :: valid_mass_fractions(:)
    logical, intent(out) :: ok_out

    type(nasa7_species), allocatable :: bad_species(:)
    real(dp) :: bad_values(2), bad_value
    real(dp) :: bad_mass_fractions(4), bad_mole_fractions(4)
    real(dp) :: output_mass_fractions(4), output_mole_fractions(4)
    real(dp) :: lower, upper, molecular_weight, gas_constant
    real(dp) :: cp, cv, gamma, enthalpy, internal_energy, entropy
    real(dp) :: pressure, density, sound_speed, recovered_temperature
    logical :: local_ok
    integer :: bad_index

    ok_out = .false.
    if (size(active_species) /= 4 .or. size(valid_mass_fractions) /= 4) return
    allocate(bad_species, source=active_species)
    bad_values(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    bad_values(2) = ieee_value(0.0_dp, ieee_positive_inf)

    do bad_index = 1, 2
      bad_value = bad_values(bad_index)

      bad_mass_fractions = valid_mass_fractions
      bad_mass_fractions(1) = bad_value
      if (valid_mixture_composition(active_species, bad_mass_fractions)) return

      bad_species = active_species
      bad_species(1)%molecular_weight = bad_value
      call mixture_temperature_bounds(bad_species, lower, upper, local_ok)
      if (local_ok .or. lower /= 0.0_dp .or. upper /= 0.0_dp) return
      molecular_weight = mixture_molecular_weight( &
        bad_species, valid_mass_fractions, local_ok)
      if (local_ok .or. molecular_weight /= 0.0_dp) return
      gas_constant = mixture_specific_gas_constant( &
        bad_species, valid_mass_fractions, local_ok)
      if (local_ok .or. gas_constant /= 0.0_dp) return

      call mixture_mass_properties( &
        active_species, valid_mass_fractions, bad_value, molecular_weight, &
        gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, local_ok)
      if (local_ok .or. .not. all([molecular_weight, gas_constant, cp, cv, &
          gamma, enthalpy, internal_energy, entropy] == 0.0_dp)) return

      pressure = mixture_pressure( &
        active_species, valid_mass_fractions, bad_value, 1200.0_dp, local_ok)
      if (local_ok .or. pressure /= 0.0_dp) return
      pressure = mixture_pressure( &
        active_species, valid_mass_fractions, 1.2_dp, bad_value, local_ok)
      if (local_ok .or. pressure /= 0.0_dp) return
      pressure = mixture_pressure( &
        active_species, valid_mass_fractions, huge(1.0_dp), &
        huge(1.0_dp), local_ok)
      if (local_ok .or. pressure /= 0.0_dp) return

      density = mixture_density( &
        active_species, valid_mass_fractions, bad_value, 1200.0_dp, local_ok)
      if (local_ok .or. density /= 0.0_dp) return
      density = mixture_density( &
        active_species, valid_mass_fractions, 101325.0_dp, bad_value, local_ok)
      if (local_ok .or. density /= 0.0_dp) return
      density = mixture_density( &
        active_species, valid_mass_fractions, huge(1.0_dp), &
        tiny(1.0_dp), local_ok)
      if (local_ok .or. density /= 0.0_dp) return

      sound_speed = mixture_sound_speed( &
        active_species, valid_mass_fractions, bad_value, local_ok)
      if (local_ok .or. sound_speed /= 0.0_dp) return

      call mixture_mass_properties( &
        active_species, valid_mass_fractions, 1200.0_dp, molecular_weight, &
        gas_constant, cp, cv, gamma, enthalpy, internal_energy, entropy, local_ok)
      if (.not. local_ok .or. .not. ieee_is_finite(internal_energy)) return
      recovered_temperature = 77.0_dp
      call temperature_from_internal_energy( &
        active_species, valid_mass_fractions, bad_value, 900.0_dp, &
        recovered_temperature, local_ok)
      if (local_ok .or. recovered_temperature /= 0.0_dp) return
      recovered_temperature = 77.0_dp
      call temperature_from_internal_energy( &
        active_species, valid_mass_fractions, internal_energy, bad_value, &
        recovered_temperature, local_ok)
      if (local_ok .or. recovered_temperature /= 0.0_dp) return

      bad_mole_fractions = valid_mass_fractions
      bad_mole_fractions(1) = bad_value
      output_mass_fractions = 1.0_dp
      call mass_fractions_from_mole_fractions( &
        active_species, bad_mole_fractions, output_mass_fractions, local_ok)
      if (local_ok .or. .not. all(output_mass_fractions == 0.0_dp)) return
      output_mole_fractions = 1.0_dp
      call mole_fractions_from_mass_fractions( &
        active_species, bad_mass_fractions, output_mole_fractions, local_ok)
      if (local_ok .or. .not. all(output_mole_fractions == 0.0_dp)) return
    end do

    bad_species = active_species
    bad_species(1)%molecular_weight = 0.125_dp * tiny(1.0_dp)
    bad_mass_fractions = 0.0_dp
    bad_mass_fractions(1) = 1.0_dp
    molecular_weight = mixture_molecular_weight( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. molecular_weight /= 0.0_dp) return
    gas_constant = mixture_specific_gas_constant( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. gas_constant /= 0.0_dp) return
    output_mole_fractions = 1.0_dp
    call mole_fractions_from_mass_fractions( &
      bad_species, bad_mass_fractions, output_mole_fractions, local_ok)
    if (local_ok .or. .not. all(output_mole_fractions == 0.0_dp)) return

    bad_species = active_species
    bad_species(1)%molecular_weight = huge(1.0_dp)
    molecular_weight = mixture_molecular_weight( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. molecular_weight /= 0.0_dp) return
    gas_constant = mixture_specific_gas_constant( &
      bad_species, bad_mass_fractions, local_ok)
    if (local_ok .or. gas_constant /= 0.0_dp) return

    ok_out = .true.
  end subroutine check_nonfinite_contract

  subroutine assert_close(actual, expected, relative_tolerance, label)
    real(dp), intent(in) :: actual, expected, relative_tolerance
    character(len=*), intent(in) :: label
    real(dp) :: error

    error = abs(actual - expected) / max(1.0_dp, abs(expected))
    if (error > relative_tolerance) then
      write(*, '(a,2(1x,es24.16),1x,es12.4)') &
        trim(label), actual, expected, error
      error stop "Mixture thermodynamics reference mismatch"
    end if
  end subroutine assert_close

end program test_mixture_thermo
