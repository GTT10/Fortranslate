! Shared chemistry identity serializer for regular and AMR checkpoints.
module reactive_checkpoint_context_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_eos_policy_identity
  use elementary_kinetics_mod, only: elementary_reaction, reaction_kind_sri, reaction_kind_sri_activated, &
    reaction_kind_tsang, reaction_kind_tsang_activated, reaction_kind_chebyshev, reaction_kind_plog
  implicit none
  private
  public :: write_selected_checkpoint_context,eos_checkpoint_header
contains
  pure function eos_checkpoint_header(base,species) result(header)
    character(len=*), intent(in) :: base
    type(nasa7_species), intent(in) :: species(:)
    character(len=:), allocatable :: header
    character(len=128) :: policy
    policy=nasa7_eos_policy_identity(species)
    header=''
    if(trim(policy)=='PELEF_INVALID_EOS_POLICY') return
    header=trim(base)
    if(len_trim(policy)>0) header=header//' '//trim(policy)
  end function eos_checkpoint_header
  subroutine write_selected_checkpoint_context( &
      unit, reactions, bundle_sha256, chemistry_integrator, &
      base_mole_fractions, relative_tolerance, absolute_tolerance, status, fixed_model)
    integer, intent(in) :: unit
    type(elementary_reaction), intent(in) :: reactions(:)
    character(len=*), intent(in) :: bundle_sha256, chemistry_integrator
    real(dp), intent(in) :: base_mole_fractions(:)
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    integer, intent(out) :: status

    logical, optional, intent(in) :: fixed_model
    character(len=16) :: context_marker
    integer :: reaction_index, efficiency_size,plog_index

    context_marker="SELECTED_CONTEXT"
    if (present(fixed_model)) then
      if (fixed_model) context_marker="FIXED_CONTEXT"
    end if
    status = 0
    write(unit, '(a)', iostat=status) trim(context_marker)
    if (status /= 0) return
    write(unit, '(a)', iostat=status) trim(bundle_sha256)
    if (status /= 0) return
    write(unit, '(a)', iostat=status) trim(chemistry_integrator)
    if (status /= 0) return
    write(unit, '(i0)', iostat=status) size(base_mole_fractions)
    if (status /= 0) return
    write(unit, '(*(es27.18e3,1x))', iostat=status) base_mole_fractions
    if (status /= 0) return
    write(unit, '(a)', iostat=status) "REACTIONS"
    if (status /= 0) return
    write(unit, '(i0)', iostat=status) size(reactions)
    if (status /= 0) return
    do reaction_index = 1, size(reactions)
      efficiency_size = 0
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) then
        efficiency_size = &
          size(reactions(reaction_index)%third_body_efficiencies)
      end if
      write(unit, '(a)', iostat=status) &
        trim(reactions(reaction_index)%equation)
      if (status /= 0) return
      write(unit, '(*(i0,1x))', iostat=status) &
        reactions(reaction_index)%kind, &
        merge(1, 0, reactions(reaction_index)%reversible), &
        merge(1, 0, &
          allocated(reactions(reaction_index)%third_body_efficiencies)), &
        efficiency_size, &
        size(reactions(reaction_index)%reactant_stoich), &
        size(reactions(reaction_index)%product_stoich)
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%reactant_stoich
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%product_stoich
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%forward_rate%pre_exponential, &
        reactions(reaction_index)%forward_rate%temperature_exponent, &
        reactions(reaction_index)%forward_rate%activation_energy, &
        reactions(reaction_index)%low_pressure_rate%pre_exponential, &
        reactions(reaction_index)%low_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%low_pressure_rate%activation_energy, &
        reactions(reaction_index)%high_pressure_rate%pre_exponential, &
        reactions(reaction_index)%high_pressure_rate%temperature_exponent, &
        reactions(reaction_index)%high_pressure_rate%activation_energy
      if (status /= 0) return
      write(unit, '(i0)', iostat=status) &
        merge(1, 0, reactions(reaction_index)%troe%enabled)
      if (status /= 0) return
      write(unit, '(*(es27.18e3,1x))', iostat=status) &
        reactions(reaction_index)%troe%alpha, &
        reactions(reaction_index)%troe%temperature_3, &
        reactions(reaction_index)%troe%temperature_1, &
        reactions(reaction_index)%troe%temperature_2
      if (status /= 0) return
      if (reactions(reaction_index)%kind==reaction_kind_sri.or. &
          reactions(reaction_index)%kind==reaction_kind_sri_activated) then
        write(unit,'(*(es27.18e3,1x))',iostat=status) reactions(reaction_index)%sri
        if (status/=0) return
      end if
      if (reactions(reaction_index)%kind==reaction_kind_tsang.or. &
          reactions(reaction_index)%kind==reaction_kind_tsang_activated) then
        write(unit,'(*(es27.18e3,1x))',iostat=status) reactions(reaction_index)%tsang
        if (status/=0) return
      end if
      if (reactions(reaction_index)%kind==reaction_kind_chebyshev) then
        write(unit,'(2(i0,1x))',iostat=status) reactions(reaction_index)%chebyshev_shape
        if (status/=0) return
        write(unit,'(*(es27.18e3,1x))',iostat=status) reactions(reaction_index)%chebyshev_bounds
        if (status/=0) return
        write(unit,'(*(es27.18e3,1x))',iostat=status) reactions(reaction_index)%chebyshev
        if (status/=0) return
      end if
      if (reactions(reaction_index)%kind==reaction_kind_plog) then
        write(unit,'(i0)',iostat=status) reactions(reaction_index)%plog_count
        if (status/=0) return
        do plog_index=1,reactions(reaction_index)%plog_count
          write(unit,'(*(es27.18e3,1x))',iostat=status) reactions(reaction_index)%plog(:,plog_index)
          if (status/=0) return
        end do
      end if
      if (allocated(reactions(reaction_index)%third_body_efficiencies)) then
        write(unit, '(*(es27.18e3,1x))', iostat=status) &
          reactions(reaction_index)%third_body_efficiencies
        if (status /= 0) return
      end if
    end do
    write(unit, '(a)', iostat=status) "CHEMISTRY_CONTROLS"
    if (status /= 0) return
    write(unit, '(*(es27.18e3,1x))', iostat=status) &
      relative_tolerance, absolute_tolerance
    if (status /= 0) return
    write(unit, '(a)', iostat=status) "END_"//trim(context_marker)
  end subroutine write_selected_checkpoint_context
end module reactive_checkpoint_context_mod
