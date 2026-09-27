module les_regular_application_mod
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use reactive_3d_application_mod, only: run_reactive_3d_application,read_regular_les_controls
  use les_smagorinsky_3d_mod, only: compute_les_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_boundary_rhs_3d
  implicit none
  private
  public :: run_regular_les_application_3d
contains
  subroutine run_regular_les_application_3d(input_path,application_label,bundle_sha256,config,species,reactions, &
      transport,base_mole_fractions,chemistry_integrator,evaluate_reaction)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    character(len=*), intent(in) :: input_path,application_label,bundle_sha256
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: base_mole_fractions(:)
    logical :: enabled,harmonic_mean,pinned_heat_flux
    real(dp) :: cs,prt,sct,ci,pressure(6)
    character(len=32) :: model,test_model
    integer :: test_ratio,grid_ratio
    call read_regular_les_controls(input_path,enabled,cs,prt,sct,ci,model,test_model,test_ratio,grid_ratio, &
      harmonic_mean,pinned_heat_flux)
    if(.not.enabled) then
      call run_reactive_3d_application(input_path,application_label,bundle_sha256,config,species,reactions, &
        transport,base_mole_fractions,chemistry_integrator,evaluate_reaction)
      return
    end if
    if(config%time_integrator/='sdc') call fatal('LES requires SDC')
    write(*,'(a,a)') 'LES model: ',trim(model)
    pressure=0
    where(config%inflow_mach>=0) pressure=config%inflow_pressure
    if(model/='dynamic') then
      if(test_model/='box'.or.test_ratio/=2.or.grid_ratio/=1.or.harmonic_mean) call fatal('Filter controls require dynamic LES')
    end if
    call run_reactive_3d_application(input_path,application_label,bundle_sha256,config,species,reactions, &
      transport,base_mole_fractions,chemistry_integrator,evaluate_reaction,evaluate)
  contains
    subroutine evaluate(u,t,h,faces,inlet,inlet_t,rates,boundary_rate,diffusivity,ok)
      real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),h(3),inlet(:,:),inlet_t(:)
      character(len=*), intent(in) :: faces(6)
      real(dp), intent(out) :: rates(:,:,:,:),boundary_rate(:),diffusivity
      logical, intent(out) :: ok
      if(model=='dynamic') then
        if(ci/=0) then
          ok=.false.;return
        end if
        call compute_dynamic_les_boundary_rhs_3d(species,u,t,h,test_model,test_ratio,rates,diffusivity,ok, &
          grid_ratio=grid_ratio,harmonic_mean=harmonic_mean,boundary_faces=faces,boundary_rate=boundary_rate, &
          inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
          wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
          inflow_mach=config%inflow_mach,inflow_pressure=pressure,pinned_heat_flux=pinned_heat_flux)
      else
        call compute_les_rhs_3d(species,u,t,h,cs,prt,sct,rates,diffusivity,ok, &
          boundary_faces=faces,boundary_rate=boundary_rate,model=model,ci=ci, &
          inflow_state=inlet,inflow_temperature=inlet_t,wall_temperature=config%wall_temperature, &
          wall_velocity=config%wall_velocity,outlet_pressure=config%outlet_pressure, &
          inflow_mach=config%inflow_mach,inflow_pressure=pressure)
      end if
    end subroutine
    subroutine fatal(message)
      character(len=*), intent(in) :: message
      write(*,'(a)') message
      error stop 2
    end subroutine
  end subroutine
end module
