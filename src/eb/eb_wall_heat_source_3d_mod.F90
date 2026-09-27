module eb_wall_heat_source_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use mixture_transport_mod, only: mixture_transport_coefficients
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_mass_fraction_component
  use eb_geometry_3d_mod, only: eb_geometry_3d,eb_cut_cell_3d,eb_covered_cell_3d
  use eb_boundary_gradient_3d_mod, only: eb_boundary_gradient_3d, &
    build_eb_quadratic_gradient_3d,apply_eb_boundary_flux_3d
  implicit none
  private
  public :: compute_eb_wall_heat_source_3d
contains
  ! Geometry/state include same-resolution donor halos. Only lower:upper cells
  ! contribute to the source/ledger. Missing active donors reject the transaction.
  ! heat_rate is energy per fluid volume per second; wall_power is signed heat
  ! added to the gas. frequency bounds the absolute frozen-conductivity row sum
  ! after division by rho*cv; it is not a nonlinear positivity guarantee.
  subroutine compute_eb_wall_heat_source_3d(species,transport,state,temperature,geometry,lower,upper, &
      wall_temperature,heat_rate,wall_power,frequency,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),wall_temperature
    integer, intent(in) :: lower(3),upper(3)
    real(dp), intent(inout) :: heat_rate(:,:,:),wall_power,frequency
    logical, intent(out) :: ok
    type(eb_boundary_gradient_3d) :: stencil
    real(dp), allocatable :: rates(:,:,:),temps(:,:,:)
    logical, allocatable :: recovered(:,:,:)
    real(dp) :: q(reactive_nprim(size(species))),y(size(species)),diffusion(size(species))
    real(dp) :: spacing(3),origin(3),center(3),centroid(3),normal(3),area,fluid_volume
    real(dp) :: donors(0:2,0:2,0:2),mu,conductivity,sound,t,flux,power,bound,row_sum
    real(dp) :: mw,rgas,cp,cv,gamma,enthalpy,energy,entropy
    logical :: covered(0:2,0:2,0:2),valid
    integer :: n(3),i,j,k,a,b,c,s,cell(3),out(3),donor(3),ios
    ok=.false.
    if(.not.geometry%is_valid().or.size(species)<1) return
    n=[geometry%nx,geometry%ny,geometry%nz]
    if(size(transport)/=size(species).or.any(shape(state)/=[reactive_nvar(size(species)),n])) return
    if(any(shape(temperature)/=n).or.any(lower<1).or.any(upper>n).or.any(lower>upper)) return
    if(any(shape(heat_rate)/=upper-lower+1)) return
    if(.not.ieee_is_finite(wall_temperature).or.wall_temperature<=0) return
    spacing=[geometry%dx,geometry%dy,geometry%dz]
    origin=[geometry%x_lower,geometry%y_lower,geometry%z_lower]
    allocate(rates(upper(1)-lower(1)+1,upper(2)-lower(2)+1,upper(3)-lower(3)+1), &
      temps(n(1),n(2),n(3)),recovered(n(1),n(2),n(3)),stat=ios)
    if(ios/=0) return
    rates=0;recovered=.false.;power=0;bound=0
    do k=lower(3),upper(3)
      do j=lower(2),upper(2)
        do i=lower(1),upper(1)
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
          donors=0;covered=.true.
          do a=0,2
            do b=0,2
              do c=0,2
                if(abs(stencil%weights(a,b,c))<=1e-14_dp) cycle
                donor=cell+stencil%offset+[a,b,c]
                if(any(donor<1).or.any(donor>n)) return
                if(geometry%cell_type(donor(1),donor(2),donor(3))==eb_covered_cell_3d) return
                if(.not.recovered(donor(1),donor(2),donor(3))) then
                  call reactive_conserved_to_primitive(species,state(:,donor(1),donor(2),donor(3)), &
                    temperature(donor(1),donor(2),donor(3)),q,t,sound,valid)
                  if(.not.valid) return
                  temps(donor(1),donor(2),donor(3))=t;recovered(donor(1),donor(2),donor(3))=.true.
                end if
                donors(a,b,c)=temps(donor(1),donor(2),donor(3));covered(a,b,c)=.false.
              end do
            end do
          end do
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),q,t,sound,valid)
          if(.not.valid) return
          do s=1,size(species)
            y(s)=q(reactive_mass_fraction_component(s))
          end do
          call mixture_transport_coefficients(species,transport,y,t,q(5),mu,conductivity,diffusion,valid)
          if(.not.valid) return
          call apply_eb_boundary_flux_3d(stencil,donors,covered,wall_temperature,conductivity,flux,valid)
          if(.not.valid) return
          fluid_volume=geometry%volume_fraction(i,j,k)*product(spacing)
          rates(out(1),out(2),out(3))=-flux/fluid_volume
          power=power-flux
          call mixture_mass_properties(species,y,t,mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,valid)
          if(.not.valid.or.cv<=0) return
          row_sum=sum(abs(stencil%weights),mask=abs(stencil%weights)>1e-14_dp)+abs(stencil%wall_weight)
          bound=max(bound,conductivity*row_sum/(fluid_volume*q(1)*cv))
        end do
      end do
    end do
    if(.not.all(ieee_is_finite(rates)).or..not.all(ieee_is_finite([power,bound]))) return
    heat_rate=rates;wall_power=power;frequency=bound;ok=.true.
  end subroutine
end module
