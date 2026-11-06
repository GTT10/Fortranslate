module eb_wall_heat_source_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients,simple_bulk_viscosity
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_mass_fraction_component
  use eb_geometry_3d_mod, only: eb_geometry_3d,eb_cut_cell_3d,eb_covered_cell_3d
  use eb_boundary_gradient_3d_mod, only: eb_boundary_gradient_3d, &
    build_eb_quadratic_gradient_3d,apply_eb_boundary_flux_3d,apply_eb_boundary_viscous_flux_3d
  use reactive_transport_3d_mod, only: valid_constant_transport_coefficients
  implicit none
  private
  public :: compute_eb_wall_heat_source_3d,compute_eb_wall_transport_source_3d
contains
  ! Geometry/state include same-resolution donor halos. Only lower:upper cells
  ! contribute to the source/ledger. Missing active donors reject the transaction.
  ! Rates are per fluid volume; the ledger is signed transfer to gas.
  ! Frequency bounds frozen heat/viscous row sums, not nonlinear positivity.
  subroutine compute_eb_wall_transport_source_3d(species,transport,state,temperature,geometry,lower,upper, &
      conduction,noslip,wall_temperature,transport_rate,wall_exchange,frequency,ok, &
      constant_coefficients,internal_energy_density)
    real(dp), optional, intent(in) :: constant_coefficients(:),internal_energy_density(:,:,:)
    real(dp), allocatable :: held
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),wall_temperature
    integer, intent(in) :: lower(3),upper(3)
    logical, intent(in) :: conduction,noslip
    real(dp), intent(inout) :: transport_rate(:,:,:,:),wall_exchange(:),frequency
    logical, intent(out) :: ok
    type(eb_boundary_gradient_3d) :: stencil
    real(dp), allocatable :: rates(:,:,:,:),temps(:,:,:),velocities(:,:,:,:)
    logical, allocatable :: recovered(:,:,:)
    real(dp) :: q(reactive_nprim(size(species))),y(size(species)),diffusion(size(species))
    real(dp) :: spacing(3),origin(3),center(3),centroid(3),normal(3),area,fluid_volume
    real(dp) :: donors(0:2,0:2,0:2),mu,conductivity,sound,t,flux,power,bound,row_sum
    real(dp) :: mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,bulk,traction(3),wall_force(3)
    real(dp) :: velocity_donors(3,0:2,0:2,0:2)
    logical :: covered(0:2,0:2,0:2),valid
    integer :: n(3),i,j,k,a,b,c,s,cell(3),out(3),donor(3),ios
    ok=.false.
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if(.not.geometry%is_valid().or.size(species)<1) return
    n=[geometry%nx,geometry%ny,geometry%nz]
    if(size(transport)/=size(species).or.any(shape(state)/=[reactive_nvar(size(species)),n])) return
    if(any(shape(temperature)/=n).or.any(lower<1).or.any(upper>n).or.any(lower>upper)) return
    if(any(shape(transport_rate)/=[reactive_nvar(size(species)),upper-lower+1])) return
    if(size(wall_exchange)/=reactive_nvar(size(species))) return
    if(conduction) then
      if(.not.ieee_is_finite(wall_temperature).or.wall_temperature<=0) return
    end if
    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=n)) return
      allocate(held)
    endif
    spacing=[geometry%dx,geometry%dy,geometry%dz]
    origin=[geometry%x_lower,geometry%y_lower,geometry%z_lower]
    allocate(rates(reactive_nvar(size(species)),upper(1)-lower(1)+1,upper(2)-lower(2)+1,upper(3)-lower(3)+1), &
      temps(n(1),n(2),n(3)),recovered(n(1),n(2),n(3)),stat=ios)
    if(ios/=0) return
    if(noslip) then
      allocate(velocities(3,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) return
    end if
    rates=0;recovered=.false.;power=0;wall_force=0;bound=0
    do k=lower(3),upper(3)
      do j=lower(2),upper(2)
        do i=lower(1),upper(1)
          if(.not.conduction.and..not.noslip) cycle
          if(geometry%cell_type(i,j,k)/=eb_cut_cell_3d) cycle
          cell=[i,j,k];out=cell-lower+1
          center=origin+(real(cell,dp)-.5_dp)*spacing
          centroid=([geometry%boundary_centroid_x(i,j,k),geometry%boundary_centroid_y(i,j,k), &
            geometry%boundary_centroid_z(i,j,k)]-center)/spacing
          ! Reference area and normal come from the aperture-vector magnitude,
          ! not the sum of triangle areas on a nonplanar cut surface.
          normal=[geometry%boundary_normal_integral_x(i,j,k),geometry%boundary_normal_integral_y(i,j,k), &
            geometry%boundary_normal_integral_z(i,j,k)]
          area=norm2(normal)
          if(area<=0) return
          normal=normal/area
          call build_eb_quadratic_gradient_3d(normal,centroid,area/spacing(1)**2,spacing,stencil,valid)
          if(.not.valid) return
          donors=0;velocity_donors=0;covered=.true.
          do a=0,2
            do b=0,2
              do c=0,2
                if(abs(stencil%weights(a,b,c))<=1e-14_dp) cycle
                donor=cell+stencil%offset+[a,b,c]
                if(any(donor<1).or.any(donor>n)) return
                if(geometry%cell_type(donor(1),donor(2),donor(3))==eb_covered_cell_3d) return
                if(.not.recovered(donor(1),donor(2),donor(3))) then
                  if(allocated(held)) held=internal_energy_density(donor(1),donor(2),donor(3))
                  call reactive_conserved_to_primitive(species,state(:,donor(1),donor(2),donor(3)), &
                    temperature(donor(1),donor(2),donor(3)),q,t,sound,valid,internal_energy_density=held)
                  if(.not.valid) return
                  if(noslip) velocities(:,donor(1),donor(2),donor(3))=q(2:4)
                  temps(donor(1),donor(2),donor(3))=t;recovered(donor(1),donor(2),donor(3))=.true.
                end if
                if(noslip) velocity_donors(:,a,b,c)=velocities(:,donor(1),donor(2),donor(3))
                donors(a,b,c)=temps(donor(1),donor(2),donor(3));covered(a,b,c)=.false.
              end do
            end do
          end do
          if(allocated(held)) held=internal_energy_density(i,j,k)
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),q,t,sound,valid, &
            internal_energy_density=held)
          if(.not.valid) return
          do s=1,size(species)
            y(s)=q(reactive_mass_fraction_component(s))
          end do
          if(present(constant_coefficients)) then
            mu=constant_coefficients(1);conductivity=constant_coefficients(2)
          else
            call mixture_transport_coefficients(species,transport,y,t,q(5),mu,conductivity,diffusion,valid)
            if(.not.valid) return
          endif
          fluid_volume=geometry%volume_fraction(i,j,k)*product(spacing)
          row_sum=sum(abs(stencil%weights),mask=abs(stencil%weights)>1e-14_dp)+abs(stencil%wall_weight)
          if(conduction) then
            call apply_eb_boundary_flux_3d(stencil,donors,covered,wall_temperature,conductivity,flux,valid)
            if(.not.valid) return
            rates(5,out(1),out(2),out(3))=-flux/fluid_volume
            power=power-flux
            call mixture_mass_properties(species,y,t,mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,valid)
            if(.not.valid.or.cv<=0) return
            bound=max(bound,conductivity*row_sum/(fluid_volume*q(1)*cv))
          end if
          if(noslip) then
            bulk=0
            if(present(constant_coefficients)) then
              if(size(constant_coefficients)==4) bulk=constant_coefficients(3)
            else if(allocated(transport(1)%simple_binary)) then
              call simple_bulk_viscosity(species,transport,y,t,bulk,valid)
              if(.not.valid) return
            end if
            call apply_eb_boundary_viscous_flux_3d(stencil,velocity_donors,covered,normal,[0._dp,0._dp,0._dp], &
              mu,bulk,traction,valid)
            if(.not.valid) return
            rates(2:4,out(1),out(2),out(3))=-traction/fluid_volume
            wall_force=wall_force-traction
            ! Stationary no-slip wall has zero mechanical energy transfer.
            bound=max(bound,(4._dp/3._dp*mu+bulk)*row_sum/(fluid_volume*q(1)))
          end if
        end do
      end do
    end do
    if(.not.all(ieee_is_finite(rates)).or..not.all(ieee_is_finite([power,wall_force,bound]))) return
    transport_rate=rates;wall_exchange=0
    wall_exchange(2:4)=wall_force;wall_exchange(5)=power;frequency=bound;ok=.true.
  end subroutine
  ! Compatibility entry for heat-only callers; commit all outputs together.
  subroutine compute_eb_wall_heat_source_3d(species,transport,state,temperature,geometry,lower,upper, &
      wall_temperature,heat_rate,wall_power,frequency,ok, &
      constant_coefficients,internal_energy_density)
    real(dp), optional, intent(in) :: constant_coefficients(:),internal_energy_density(:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),wall_temperature
    integer, intent(in) :: lower(3),upper(3)
    real(dp), intent(inout) :: heat_rate(:,:,:),wall_power,frequency
    logical, intent(out) :: ok
    real(dp), allocatable :: rates(:,:,:,:)
    real(dp) :: exchange(reactive_nvar(size(species))),bound
    integer :: dims(3),ios
    logical :: valid
    ok=.false.;dims=shape(heat_rate)
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if(any(dims/=upper-lower+1)) return
    allocate(rates(reactive_nvar(size(species)),dims(1),dims(2),dims(3)),stat=ios)
    if(ios/=0) return
    call compute_eb_wall_transport_source_3d(species,transport,state,temperature,geometry,lower,upper, &
      .true.,.false.,wall_temperature,rates,exchange,bound,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=internal_energy_density)
    if(.not.valid) return
    heat_rate=rates(5,:,:,:);wall_power=exchange(5);frequency=bound;ok=.true.
  end subroutine
end module
