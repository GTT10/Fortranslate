program test_eb_reactive_wall_flux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: eb_geometry_2d,build_eb_geometry_2d
  use state_indices_mod, only: irho, imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_species_component, reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, build_axis_plane_eb_geometry_3d,extrude_eb_geometry_3d, &
    build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_wall_flux_3d_mod, only: &
    reactive_eb_slip_wall_flux_3d, reactive_eb_slip_wall_source_3d, &
    reactive_eb_flux_divergence_3d,interpolate_reactive_eb_centroid_fluxes_3d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 8
  integer, parameter :: nz = 6
  real(dp), parameter :: pressure = 135000.0_dp
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geometry
  real(dp), allocatable :: primitive(:), changed_primitive(:)
  real(dp), allocatable :: state_cell(:), changed_state(:)
  real(dp), allocatable :: flux(:), changed_flux(:), mass_fractions(:)
  real(dp), allocatable :: state(:, :, :, :), temperature_field(:, :, :)
  real(dp), allocatable :: source(:, :, :, :), rhs(:, :, :, :)
  real(dp), allocatable :: x_flux(:, :, :, :), y_flux(:, :, :, :)
  real(dp), allocatable :: z_flux(:, :, :, :), short_x_flux(:, :, :, :)
  real(dp) :: mole_fractions(7), normal(3)
  real(dp) :: temperature, changed_temperature, sound_speed
  real(dp) :: changed_sound_speed, wall_pressure, changed_pressure
  real(dp) :: tolerance
  logical :: ok
  integer :: i, j, k, nvar, species_first

  call load_h2o2_elementary_thermo(species, ok)
  call require(ok, "thermodynamic database load")
  nvar = reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(changed_primitive(size(primitive)))
  allocate(state_cell(nvar), changed_state(nvar))
  allocate(flux(nvar), changed_flux(nvar), mass_fractions(size(species)))
  allocate(state(nvar, nx, ny, nz), temperature_field(nx, ny, nz))
  allocate(source(nvar, nx, ny, nz), rhs(nvar, nx, ny, nz))
  allocate(x_flux(nvar, 0:nx, ny, nz))
  allocate(y_flux(nvar, nx, 0:ny, nz))
  allocate(z_flux(nvar, nx, ny, 0:nz))

  mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
    1.0e-5_dp, 0.0_dp, 0.55643_dp]
  call mass_fractions_from_mole_fractions( &
    species, mole_fractions, mass_fractions, ok)
  call require(ok, "composition conversion")
  primitive(1:5) = [0.31_dp, 42.0_dp, -27.0_dp, 6.0_dp, pressure]
  do k = 1, size(species)
    primitive(reactive_mass_fraction_component(k)) = mass_fractions(k)
  end do
  call reactive_primitive_to_conserved( &
    species, primitive, state_cell, temperature, sound_speed, ok)
  call require(ok, "reference state construction")
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        state(:, i, j, k) = state_cell
        temperature_field(i, j, k) = temperature
      end do
    end do
  end do

  normal = [2.0_dp, 3.0_dp, 6.0_dp] / 7.0_dp
  call reactive_eb_slip_wall_flux_3d( &
    species, state_cell, temperature, normal, flux, wall_pressure, ok)
  call require(ok, "oblique wall flux")
  tolerance = 2.0e-9_dp * pressure
  call assert_close(wall_pressure, pressure, tolerance, &
    "wall pressure recovery")
  call assert_close(flux(imx), -2.0_dp * pressure / 7.0_dp, &
    tolerance, "wall x-momentum flux")
  call assert_close(flux(imy), -3.0_dp * pressure / 7.0_dp, &
    tolerance, "wall y-momentum flux")
  call assert_close(flux(imz), -6.0_dp * pressure / 7.0_dp, &
    tolerance, "wall z-momentum flux")
  species_first = reactive_species_component(1)
  call require(flux(irho) == 0.0_dp .and. flux(iet) == 0.0_dp .and. &
    maxval(abs(flux(species_first:nvar))) == 0.0_dp, &
    "wall mass energy and species impermeability")

  changed_primitive = primitive
  changed_primitive(2:4) = [-113.0_dp, 87.0_dp, -19.0_dp]
  call reactive_primitive_to_conserved( &
    species, changed_primitive, changed_state, changed_temperature, &
    changed_sound_speed, ok)
  call require(ok, "changed-velocity state construction")
  call reactive_eb_slip_wall_flux_3d( &
    species, changed_state, changed_temperature, normal, changed_flux, &
    changed_pressure, ok)
  call require(ok, "changed-velocity wall flux")
  call require(maxval(abs(changed_flux - flux)) <= tolerance, &
    "slip wall flux velocity independence")

  call verify_plane_balance("x", imx)
  call verify_plane_balance("y", imy)
  call verify_plane_balance("z", imz)
  block
    type(eb_geometry_2d) :: section
    type(eb_geometry_3d) :: fine
    real(dp) :: fine_sphere(0:2*nx,0:2*ny,0:2*nz)
    real(dp) :: phi(0:nx,0:ny),sphere(0:nx,0:ny,0:nz),weighted_residual
    integer :: mode
    do j=0,ny
      do i=0,nx
        phi(i,j)=(real(i,dp)/nx-.5_dp)**2+(real(j,dp)/ny-.5_dp)**2-.29_dp**2
      end do
    end do
    call build_eb_geometry_2d(phi,0._dp,1._dp,0._dp,1._dp,section,ok)
    call require(ok,'curved wall section')
    call extrude_eb_geometry_3d(section,3,nz,0._dp,1._dp,geometry,ok)
    call require(ok,'curved wall extrusion')
    do mode=1,3
      if (mode==2) then
        do k=0,nz
          do j=0,ny
            do i=0,nx
              sphere(i,j,k)=(real(i,dp)/nx-.5_dp)**2+(real(j,dp)/ny-.5_dp)**2 &
                +(real(k,dp)/nz-.5_dp)**2-.29_dp**2
            end do
          end do
        end do
        call build_level_set_eb_geometry_3d(sphere,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],geometry,ok)
        call require(ok,'fully 3D spherical wall')
      end if
      if (mode==3) then
        do k=0,2*nz
          do j=0,2*ny
            do i=0,2*nx
              fine_sphere(i,j,k)=(real(i,dp)/(2*nx)-.5_dp)**2+(real(j,dp)/(2*ny)-.5_dp)**2 &
                +(real(k,dp)/(2*nz)-.5_dp)**2-.29_dp**2
            end do
          end do
        end do
        call build_level_set_eb_geometry_3d(fine_sphere,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],fine,ok)
        call require(ok,'fine spherical wall')
        call coarsen_eb_geometry_3d(fine,2,geometry,ok)
        call require(ok,'restricted spherical wall')
      end if
    ! Reuse the uniform-pressure Cartesian face fluxes from the plane check.
    call reactive_eb_flux_divergence_3d(species,state,temperature_field,geometry,x_flux,y_flux,z_flux,rhs,ok)
    call require(ok,'curved wall pressure divergence')
    weighted_residual=0._dp
    do k=1,nz
      do j=1,ny
        do i=1,nx
          weighted_residual=max(weighted_residual,maxval(abs(rhs(:,i,j,k)))*geometry%volume_fraction(i,j,k) &
            *geometry%dx*geometry%dy*geometry%dz)
        end do
      end do
    end do
    call require(weighted_residual<=2.e-10_dp*pressure*max(geometry%dx*geometry%dy, &
      geometry%dx*geometry%dz,geometry%dy*geometry%dz),'curved wall uniform-pressure balance')
    end do
  end block


  call reactive_eb_slip_wall_flux_3d( &
    species, state_cell, temperature, [2.0_dp, 0.0_dp, 0.0_dp], &
    flux, wall_pressure, ok)
  call require(.not. ok .and. maxval(abs(flux)) == 0.0_dp, &
    "non-unit normal rejection")

  x_flux(irho, 0, 1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call reactive_eb_flux_divergence_3d( &
    species, state, temperature_field, geometry, x_flux, y_flux, z_flux, &
    rhs, ok)
  call require(.not. ok .and. maxval(abs(rhs)) == 0.0_dp, &
    "nonfinite face flux transaction")

  allocate(short_x_flux(nvar, 0:nx - 1, ny, nz))
  short_x_flux = 0.0_dp
  call reactive_eb_flux_divergence_3d( &
    species, state, temperature_field, geometry, short_x_flux, y_flux, &
    z_flux, rhs, ok)
  call require(.not. ok .and. maxval(abs(rhs)) == 0.0_dp, &
    "face-array extent transaction")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "transaction geometry")
  state(irho, 4, 2, 3) = -1.0_dp
  call reactive_eb_slip_wall_source_3d( &
    species, state, temperature_field, geometry, source, ok)
  call require(.not. ok .and. maxval(abs(source)) == 0.0_dp, &
    "failed wall source transaction")

  call check_centroid_interpolation()
  call check_mol_centroid_stencil()
  write(*, '(a)') "test_eb_reactive_wall_flux_3d: PASS"

contains
  subroutine check_centroid_interpolation()
    type(eb_geometry_3d) :: g
    real(dp) :: phi(0:6,0:6,0:6),x(2,0:6,6,6),y(2,6,0:6,6),z(2,6,6,0:6)
    real(dp) :: fx(2,0:4,4,4),fy(2,4,0:4,4),fz(2,4,4,0:4)
    real(dp) :: sx(2,0:4,4,4),sy(2,4,0:4,4),sz(2,4,4,0:4),offset(2),ap,val(2),expected
    integer :: mode,ii,jj,kk,axis,lo(3),hi(3),idx(3),b,c,bilinear_faces
    logical :: valid
    do mode=0,2
      do kk=0,6
      do jj=0,6
      do ii=0,6
        phi(ii,jj,kk)=1
        if(mode>0) phi(ii,jj,kk)=(ii+jj+kk-7.3_dp)*merge(1._dp,-1._dp,mode==1)
      end do
      end do
      end do
      call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[6._dp,6._dp,6._dp],g,valid)
      call require(valid,'centroid interpolation geometry')
      x=ieee_value(0._dp,ieee_quiet_nan);y=ieee_value(0._dp,ieee_quiet_nan);z=ieee_value(0._dp,ieee_quiet_nan)
      do axis=1,3
        lo=1;lo(axis)=0;hi=6;b=merge(2,1,axis==1);c=merge(2,3,axis==3)
        do kk=lo(3),hi(3)
        do jj=lo(2),hi(2)
        do ii=lo(1),hi(1)
          idx=[ii,jj,kk]
          val=[2._dp,1._dp+2*idx(b)+3*idx(c)+.5_dp*idx(b)*idx(c)]
          select case(axis)
          case(1)
            if(g%x_face_fraction(ii,jj,kk)>0) x(:,ii,jj,kk)=val
          case(2)
            if(g%y_face_fraction(ii,jj,kk)>0) y(:,ii,jj,kk)=val
          case(3)
            if(g%z_face_fraction(ii,jj,kk)>0) z(:,ii,jj,kk)=val
          end select
        end do
        end do
        end do
      end do
      fx=-17;fy=-18;fz=-19
      call interpolate_reactive_eb_centroid_fluxes_3d(g,x,y,z,fx,fy,fz,valid)
      call require(valid,'centroid interpolation ignores closed input faces')
      bilinear_faces=0
      do axis=1,3
        lo=2;lo(axis)=1;hi=5;b=merge(2,1,axis==1);c=merge(2,3,axis==3)
        do kk=lo(3),hi(3)
        do jj=lo(2),hi(2)
        do ii=lo(1),hi(1)
          idx=[ii,jj,kk]
          select case(axis)
          case(1)
            ap=g%x_face_fraction(ii,jj,kk);offset=[g%x_face_centroid_y(ii,jj,kk),g%x_face_centroid_z(ii,jj,kk)]
            val=fx(:,ii-1,jj-1,kk-1)
          case(2)
            ap=g%y_face_fraction(ii,jj,kk);offset=[g%y_face_centroid_x(ii,jj,kk),g%y_face_centroid_z(ii,jj,kk)]
            val=fy(:,ii-1,jj-1,kk-1)
          case(3)
            ap=g%z_face_fraction(ii,jj,kk);offset=[g%z_face_centroid_x(ii,jj,kk),g%z_face_centroid_y(ii,jj,kk)]
            val=fz(:,ii-1,jj-1,kk-1)
          end select
          if(ap==0) then
            call require(all(val==0),'closed centroid flux zero')
          else
            expected=1+2*(idx(b)+offset(1))+3*(idx(c)+offset(2))+.5_dp*(idx(b)+offset(1))*(idx(c)+offset(2))
            call require(abs(val(1)-2)<1e-14_dp.and.abs(val(2)-expected)<1e-12_dp, &
              'constant and bilinear field exact at face centroid')
            if(all(abs(offset)>1e-8_dp)) bilinear_faces=bilinear_faces+1
          end if
        end do
        end do
        end do
      end do
      if(mode>0) call require(bilinear_faces>0,'oblique geometry exercises both tangential offsets')
      if(mode==0) then
        sx=fx;sy=fy;sz=fz;x(:,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
        call interpolate_reactive_eb_centroid_fluxes_3d(g,x,y,z,fx,fy,fz,valid)
        call require(.not.valid.and.all(fx==sx).and.all(fy==sy).and.all(fz==sz),'centroid interpolation rollback')
      end if
    end do
  end subroutine check_centroid_interpolation

  subroutine check_mol_centroid_stencil()
    type(eb_geometry_3d) :: g
    real(dp) :: phi(0:6,0:6,0:6),x(1,0:6,6,6),y(1,6,0:6,6),z(1,6,6,0:6)
    real(dp) :: fx(1,0:4,4,4),fy(1,4,0:4,4),fz(1,4,4,0:4),saved(1,4,0:4,4)
    integer :: i,j,k
    logical :: valid
    do k=0,6
    do j=0,6
    do i=0,6
      phi(i,j,k)=real(i,dp)-1.4_dp
    end do
    end do
    end do
    call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[6._dp,6._dp,6._dp],g,valid)
    call require(valid,'MOL face stencil plane')
    x=0;y=0;z=0
    y(1,2,3,3)=2;y(1,3,3,3)=7
    z(1,2,3,3)=3;z(1,3,3,3)=8
    call interpolate_reactive_eb_centroid_fluxes_3d(g,x,y,z,fx,fy,fz,valid,mol_stencil=.true.)
    call require(valid.and.abs(fy(1,1,2,2)-3._dp)<1e-14_dp.and.abs(fz(1,1,2,2)-4._dp)<1e-14_dp, &
      'MOL pinned bilinear face weights before aperture multiplication')
    ! A large finite donor makes the pinned small-coefficient cutoff observable.
    g%y_face_centroid_x(2,3,3)=1e-15_dp
    y=0;y(1,2,3,3)=2;y(1,3,3,3)=7e14_dp
    call interpolate_reactive_eb_centroid_fluxes_3d(g,x,y,z,fx,fy,fz,valid,mol_stencil=.true.)
    call require(valid.and.fy(1,1,2,2)==2*(1-1e-15_dp),'MOL aperture-weighted cutoff')
    saved=fy;y(1,2,3,3)=ieee_value(0._dp,ieee_quiet_nan)
    call interpolate_reactive_eb_centroid_fluxes_3d(g,x,y,z,fx,fy,fz,valid,mol_stencil=.true.)
    call require(.not.valid.and.all(fy==saved),'MOL face stencil rollback')
  end subroutine check_mol_centroid_stencil

  subroutine verify_plane_balance(axis, momentum_component)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: momentum_component

    real(dp) :: cell_volume, integrated_force(3)
    real(dp) :: integrated_residual, maximum_integrated_residual
    logical :: local_ok
    integer :: local_i, local_j, local_k

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.37_dp, geometry, local_ok)
    call require(local_ok .and. geometry%is_valid(), &
      trim(axis)//"-plane geometry")
    call reactive_eb_slip_wall_source_3d( &
      species, state, temperature_field, geometry, source, local_ok)
    call require(local_ok, trim(axis)//"-plane wall source")

    cell_volume = geometry%dx * geometry%dy * geometry%dz
    integrated_force = [ &
      sum(source(imx, :, :, :) * geometry%volume_fraction) * cell_volume, &
      sum(source(imy, :, :, :) * geometry%volume_fraction) * cell_volume, &
      sum(source(imz, :, :, :) * geometry%volume_fraction) * cell_volume]
    call assert_close(integrated_force(momentum_component - imx + 1), &
      pressure, tolerance, trim(axis)//"-plane integrated force")
    integrated_force(momentum_component - imx + 1) = 0.0_dp
    call require(maxval(abs(integrated_force)) <= tolerance, &
      trim(axis)//"-plane transverse force")
    call require(maxval(abs(source(irho, :, :, :))) == 0.0_dp .and. &
      maxval(abs(source(iet, :, :, :))) == 0.0_dp .and. &
      maxval(abs(source(species_first:nvar, :, :, :))) == 0.0_dp, &
      trim(axis)//"-plane impermeable source")

    x_flux = 0.0_dp
    y_flux = 0.0_dp
    z_flux = 0.0_dp
    x_flux(imx, :, :, :) = pressure
    y_flux(imy, :, :, :) = pressure
    z_flux(imz, :, :, :) = pressure
    call reactive_eb_flux_divergence_3d( &
      species, state, temperature_field, geometry, x_flux, y_flux, z_flux, &
      rhs, local_ok)
    call require(local_ok, trim(axis)//"-plane divergence")

    maximum_integrated_residual = 0.0_dp
    do local_k = 1, nz
      do local_j = 1, ny
        do local_i = 1, nx
          if (geometry%cell_type(local_i, local_j, local_k) == &
              eb_covered_cell_3d) then
            call require(maxval(abs(rhs(:, local_i, local_j, local_k))) == &
              0.0_dp, trim(axis)//"-plane covered residual")
            cycle
          end if
          integrated_residual = &
            maxval(abs(rhs(:, local_i, local_j, local_k))) * &
            geometry%volume_fraction(local_i, local_j, local_k) * &
            cell_volume
          maximum_integrated_residual = max( &
            maximum_integrated_residual, integrated_residual)
        end do
      end do
    end do
    call require(maximum_integrated_residual <= &
      2.0e-10_dp * pressure * max( &
        geometry%dx * geometry%dy, geometry%dx * geometry%dz, &
        geometry%dy * geometry%dz), &
      trim(axis)//"-plane uniform-pressure balance")
  end subroutine verify_plane_balance

  subroutine assert_close(actual, expected, local_tolerance, message)
    real(dp), intent(in) :: actual, expected, local_tolerance
    character(len=*), intent(in) :: message

    call require(abs(actual - expected) <= local_tolerance, message)
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_wall_flux_3d
