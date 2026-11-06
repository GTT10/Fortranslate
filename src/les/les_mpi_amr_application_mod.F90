module les_mpi_amr_application_mod
  use mpi_f08
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_application_3d_mod, only: run_mpi_amr_patch_tree_application_3d
  use les_mpi_amr_tree_mod, only: compute_mpi_amr_tree_les_rhs_3d
  implicit none
  private
  public :: run_mpi_amr_les_application_3d
contains
  subroutine run_mpi_amr_les_application_3d(topology_path,input_path,config,species,reactions,transport, &
      mole_fractions,comm,chemistry_integrator,mechanism_identity,schedule_path,evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    character(len=*), intent(in) :: topology_path,input_path
    character(len=*), optional, intent(in) :: chemistry_integrator,mechanism_identity,schedule_path
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(nasa7_species), allocatable :: eos_species(:)
    logical :: eos_policy_ok
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(MPI_Comm), intent(in) :: comm
    logical :: enabled,harmonic_mean
    real(dp) :: cs,prt,sct,ci,pressure(6)
    character(len=32) :: model,test_model
    character(len=512) :: identity
    integer :: test_ratio,grid_ratio,unit,ios,status
    character(len=:), allocatable :: dynamic_filter
    integer, allocatable :: dynamic_ratio,dynamic_grid
    logical, allocatable :: dynamic_harmonic
    namelist /amr_tree_les/ enabled,cs,prt,sct,ci,model,test_model,test_ratio,grid_ratio,harmonic_mean
    eos_species=species
    call configure_nasa7_eos_policy(eos_species,config%pinned_eos_extrapolation,eos_policy_ok, &
      config%pinned_eos_temperature_min,config%pinned_eos_temperature_max)
    if(.not.eos_policy_ok) error stop 'Invalid application EOS policy'
    enabled=.false.;cs=.17_dp;prt=.7_dp;sct=.7_dp;ci=0
    model='smagorinsky';test_model='box';test_ratio=2;grid_ratio=1;harmonic_mean=.false.
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) call fatal('Cannot open LES input')
    read(unit,nml=amr_tree_les,iostat=ios)
    close(unit)
    if(ios>0) call fatal('Invalid amr_tree_les namelist')
    if(.not.enabled) then
      call run_mpi_amr_patch_tree_application_3d(topology_path,input_path,config,eos_species,reactions,transport, &
        mole_fractions,comm,chemistry_integrator,mechanism_identity,schedule_path,handled_input_group="amr_tree_les", &
        evaluate_reaction=evaluate_reaction)
      return
    end if
    pressure=0
    where(config%inflow_mach>=0) pressure=config%inflow_pressure
    if(model=='dynamic') then
      dynamic_filter=trim(test_model);dynamic_ratio=test_ratio;dynamic_grid=grid_ratio;dynamic_harmonic=harmonic_mean
    else
      if(test_model/='box'.or.test_ratio/=2.or.grid_ratio/=1.or.harmonic_mean) call fatal('Filter controls require dynamic LES')
    end if
    write(identity,'(a,4es25.17,2a,2i6,l2)') 'les-sdc-v1:',cs,prt,sct,ci,trim(model)//':', &
      trim(test_model)//':',test_ratio,grid_ratio,harmonic_mean
    call run_mpi_amr_patch_tree_application_3d(topology_path,input_path,config,eos_species,reactions,transport, &
      mole_fractions,comm,chemistry_integrator,mechanism_identity,schedule_path,evaluate,trim(identity),"amr_tree_les", &
      evaluate_reaction=evaluate_reaction)
  contains
    subroutine evaluate(tree,owners,u,t,h,faces,inlet,inlet_t,depth,fractions,coarse_start,coarse_end,rates,flux,diffusivity,ok, &
        frame,domain_lower,internal_energy)
      type(amr_tree_state_3d), optional, intent(in) :: internal_energy
      type(rotating_frame_3d), optional, intent(in) :: frame
      real(dp), optional, intent(in) :: domain_lower(3)
      type(amr_patch_tree_3d), intent(in) :: tree
      integer, intent(in) :: owners(:),depth
      type(amr_tree_state_3d), intent(in) :: u,t
      real(dp), intent(in) :: h(3),inlet(:,:),inlet_t(:)
      character(len=*), intent(in) :: faces(6)
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
      type(amr_tree_state_3d), intent(inout) :: rates
      type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: diffusivity
      logical, intent(out) :: ok
      call compute_mpi_amr_tree_les_rhs_3d(tree,owners,eos_species,u,t,h,faces,cs,prt,sct,rates,flux,diffusivity,comm,ok, &
        inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
        wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure,model=model,ci=ci, &
        test_model=dynamic_filter,test_ratio=dynamic_ratio,grid_ratio=dynamic_grid,harmonic_mean=dynamic_harmonic, &
        inflow_mach=config%inflow_mach,inflow_pressure=pressure, &
        target_depth=depth,coarse_start=coarse_start,coarse_end=coarse_end,fractions=fractions, &
        frame=frame,domain_lower=domain_lower,internal_energy=internal_energy)
    end subroutine
    subroutine fatal(message)
      character(len=*), intent(in) :: message
      write(*,'(a)') message
      call MPI_Abort(comm,2,status)
      error stop 2
    end subroutine
  end subroutine
end module
