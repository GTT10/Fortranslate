module h2o2_full_thermo_mod
  use h2o2_full_mechanism_mod, only: &
    full_nspecies => h2o2_full_nspecies, &
    full_h2_index => h2o2_full_h2_index, &
    full_h_index => h2o2_full_h_index, &
    full_o_index => h2o2_full_o_index, &
    full_o2_index => h2o2_full_o2_index, &
    full_oh_index => h2o2_full_oh_index, &
    full_h2o_index => h2o2_full_h2o_index, &
    full_ho2_index => h2o2_full_ho2_index, &
    full_h2o2_index => h2o2_full_h2o2_index, &
    full_ar_index => h2o2_full_ar_index, &
    full_n2_index => h2o2_full_n2_index, &
    load_h2o2_full_thermo => load_h2o2_full_thermo_data
  implicit none
  private

  public :: full_nspecies
  public :: full_h2_index, full_h_index, full_o_index, full_o2_index
  public :: full_oh_index, full_h2o_index, full_ho2_index
  public :: full_h2o2_index, full_ar_index, full_n2_index
  public :: load_h2o2_full_thermo

end module h2o2_full_thermo_mod
