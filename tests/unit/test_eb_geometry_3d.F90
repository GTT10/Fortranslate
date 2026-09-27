program test_eb_geometry_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: eb_geometry_2d,build_eb_geometry_2d
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    eb_regular_cell_3d, build_axis_plane_eb_geometry_3d,extrude_eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d
  implicit none

  integer, parameter :: nx = 10
  integer, parameter :: ny = 8
  integer, parameter :: nz = 6
  real(dp), parameter :: tolerance = 4.0e-13_dp
  type(eb_geometry_3d) :: geometry
  real(dp) :: fluid_volume, nan_value
  logical :: ok

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", -0.1_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "regular geometry validity")
  call require(all(geometry%cell_type == eb_regular_cell_3d), &
    "regular classification")
  call require(maxval(abs(geometry%volume_fraction - 1.0_dp)) == 0.0_dp, &
    "regular volume fractions")
  call require(maxval(abs(geometry%x_face_fraction - 1.0_dp)) == 0.0_dp &
      .and. maxval(abs(geometry%y_face_fraction - 1.0_dp)) == 0.0_dp &
      .and. maxval(abs(geometry%z_face_fraction - 1.0_dp)) == 0.0_dp, &
    "regular face fractions")
  call require(maxval(abs(geometry%boundary_area)) == 0.0_dp, &
    "regular boundary area")
  call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
    tolerance, "regular geometric identity")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.0_dp, geometry, ok)
  call require(ok .and. all(geometry%cell_type == eb_regular_cell_3d), &
    "domain-lower plane regular collapse")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "z", 1.1_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "covered geometry validity")
  call require(all(geometry%cell_type == eb_covered_cell_3d), &
    "covered classification")
  call require(maxval(abs(geometry%volume_fraction)) == 0.0_dp, &
    "covered volume fractions")
  call require(maxval(abs(geometry%x_face_fraction)) == 0.0_dp .and. &
      maxval(abs(geometry%y_face_fraction)) == 0.0_dp .and. &
      maxval(abs(geometry%z_face_fraction)) == 0.0_dp, &
    "covered face fractions")
  call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
    tolerance, "covered geometric identity")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "z", 1.0_dp, geometry, ok)
  call require(ok .and. all(geometry%cell_type == eb_covered_cell_3d), &
    "domain-upper plane covered collapse")

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok .and. geometry%is_valid(), "x-plane geometry validity")
  fluid_volume = sum(geometry%volume_fraction) * &
    geometry%dx * geometry%dy * geometry%dz
  call assert_close(fluid_volume, 0.63_dp, tolerance, &
    "x-plane fluid volume")
  call require(count(geometry%cell_type == eb_cut_cell_3d) == ny * nz, &
    "x-plane cut-cell count")
  call assert_close(maxval(abs(geometry%volume_fraction(4, :, :) - &
    0.30_dp)), 0.0_dp, tolerance, "x-plane cut fraction")
  call assert_close(maxval(abs(geometry%cell_centroid_x(4, :, :) - &
    0.35_dp)), 0.0_dp, tolerance, "x-plane cell centroid")
  call require(maxval(abs(geometry%cell_centroid_y)) == 0.0_dp .and. &
    maxval(abs(geometry%cell_centroid_z)) == 0.0_dp, &
    "x-plane transverse cell centroids")
  call require(all(geometry%x_face_fraction(3, :, :) == 0.0_dp) .and. &
    all(geometry%x_face_fraction(4, :, :) == 1.0_dp), &
    "x-plane normal face fractions")
  call assert_close(maxval(abs(geometry%y_face_fraction(4, :, :) - &
    0.30_dp)), 0.0_dp, tolerance, "x-plane y-face fraction")
  call assert_close(maxval(abs(geometry%z_face_fraction(4, :, :) - &
    0.30_dp)), 0.0_dp, tolerance, "x-plane z-face fraction")
  call assert_close(maxval(abs(geometry%y_face_centroid_x(4, :, :) - &
    0.35_dp)), 0.0_dp, tolerance, "x-plane y-face centroid")
  call assert_close(maxval(abs(geometry%z_face_centroid_x(4, :, :) - &
    0.35_dp)), 0.0_dp, tolerance, "x-plane z-face centroid")
  call assert_close(sum(geometry%boundary_area), 1.0_dp, tolerance, &
    "x-plane boundary area")
  call assert_close(maxval(abs(geometry%boundary_centroid_x(4, :, :) - &
    0.37_dp)), 0.0_dp, tolerance, "x-plane boundary centroid")
  call assert_close(maxval(abs(geometry%boundary_normal_x(4, :, :) - &
    1.0_dp)), 0.0_dp, tolerance, "x-plane boundary normal")
  call require(maxval(abs(geometry%boundary_normal_y)) == 0.0_dp .and. &
    maxval(abs(geometry%boundary_normal_z)) == 0.0_dp, &
    "x-plane transverse boundary normals")
  call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
    tolerance, "x-plane geometric identity")

  geometry%boundary_normal_integral_x(4, 2, 3) = &
    0.9_dp * geometry%boundary_normal_integral_x(4, 2, 3)
  call require(.not. geometry%is_valid(), &
    "normal-integral corruption rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "x-plane rebuild")
  geometry%y_face_fraction(4, 2, 3) = 0.4_dp
  call require(.not. geometry%is_valid(), "face-divergence corruption rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(ok, "x-plane second rebuild")
  geometry%y_face_centroid_x(4, 2, 3) = 0.6_dp
  call require(.not. geometry%is_valid(), "face-centroid corruption rejection")

  call verify_rotated_plane("y", nx * nz)
  call verify_rotated_plane("z", nx * ny)

  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "q", 0.37_dp, geometry, ok)
  call require(.not. ok, "invalid axis rejection")
  call build_axis_plane_eb_geometry_3d( &
    0, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(.not. ok, "invalid extent rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.37_dp, geometry, ok)
  call require(.not. ok, "invalid bound rejection")
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", 0.50_dp, geometry, ok)
  call require(.not. ok, "face-aligned interior plane rejection")
  nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
  call build_axis_plane_eb_geometry_3d( &
    nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
    "x", nan_value, geometry, ok)
  call require(.not. ok, "nonfinite plane rejection")

  block
    type(eb_geometry_2d) :: section
    type(eb_geometry_3d) :: prism,prior
    real(dp) :: phi(0:10,0:8),expected_volume,expected_area,moment(3),center(3),expected_center(3)
    real(dp) :: cell_volume
    integer :: ii,jj,kk,direction
    do jj=0,8
      do ii=0,10
        phi(ii,jj)=real(ii,dp)/10+real(jj,dp)/8-1.13_dp
      end do
    end do
    call build_eb_geometry_2d(phi,0._dp,1._dp,0._dp,1._dp,section,ok)
    call require(ok,'oblique section')
    expected_volume=.5_dp*.87_dp**2*2.5_dp
    expected_area=sqrt(2._dp)*.87_dp*2.5_dp
    do direction=1,3
      call extrude_eb_geometry_3d(section,direction,3,2._dp,4.5_dp,prism,ok)
      call require(ok.and.prism%is_valid(),'three-axis oblique extrusion')
      call assert_close(sum(prism%volume_fraction)*prism%dx*prism%dy*prism%dz, &
        expected_volume,tolerance,'analytic triangular-prism volume')
      call assert_close(sum(prism%boundary_area),expected_area,tolerance,'analytic oblique wall area')
      moment=0._dp;cell_volume=prism%dx*prism%dy*prism%dz
      do kk=1,prism%nz
        do jj=1,prism%ny
          do ii=1,prism%nx
            center=[prism%x_lower+(real(ii,dp)-.5_dp+prism%cell_centroid_x(ii,jj,kk))*prism%dx, &
              prism%y_lower+(real(jj,dp)-.5_dp+prism%cell_centroid_y(ii,jj,kk))*prism%dy, &
              prism%z_lower+(real(kk,dp)-.5_dp+prism%cell_centroid_z(ii,jj,kk))*prism%dz]
            moment=moment+center*prism%volume_fraction(ii,jj,kk)*cell_volume
          end do
        end do
      end do
      expected_center=.71_dp;expected_center(direction)=3.25_dp
      call assert_close(maxval(abs(moment/expected_volume-expected_center)),0._dp,tolerance, &
        'analytic triangular-prism centroid')
      call assert_close(maximum_geometric_residual(prism),0._dp,tolerance,'extruded closure')
    end do
    do jj=0,8
      do ii=0,10
        phi(ii,jj)=(real(ii,dp)/10-.5_dp)**2+(real(jj,dp)/8-.5_dp)**2-.29_dp**2
      end do
    end do
    call build_eb_geometry_2d(phi,0._dp,1._dp,0._dp,1._dp,section,ok)
    call require(ok,'faceted cylindrical section')
    do direction=1,3
      call extrude_eb_geometry_3d(section,direction,3,2._dp,4.5_dp,prism,ok)
      call require(ok.and.prism%is_valid(),'faceted extrusion')
      call assert_close(sum(prism%volume_fraction)*prism%dx*prism%dy*prism%dz, &
        sum(section%volume_fraction)*section%dx*section%dy*2.5_dp,tolerance,'section volume extrusion')
      call assert_close(sum(prism%boundary_area),sum(section%boundary_length)*2.5_dp,tolerance,'faceted area')
      call assert_close(maximum_geometric_residual(prism),0._dp,tolerance,'faceted closure')
    end do
    prior=prism
    call extrude_eb_geometry_3d(section,0,3,2._dp,4.5_dp,prism,ok)
    call require(.not.ok.and.all(prism%volume_fraction==prior%volume_fraction),'invalid axis preserves geometry')
    call extrude_eb_geometry_3d(section,3,0,2._dp,4.5_dp,prism,ok)
    call require(.not.ok.and.prism%is_valid(),'invalid layers preserves geometry')
    call extrude_eb_geometry_3d(section,3,3,4.5_dp,2._dp,prism,ok)
    call require(.not.ok.and.prism%is_valid(),'invalid bounds preserves geometry')
    prism%boundary_area=.5_dp*prism%boundary_area
    call require(.not.prism%is_valid(),'area below integrated normal rejected')
    section%volume_fraction(1,1)=-1._dp
    call extrude_eb_geometry_3d(section,3,3,2._dp,4.5_dp,prior,ok)
    call require(.not.ok.and.prior%is_valid(),'invalid section preserves geometry')
  end block

  block
    type(eb_geometry_3d) :: symmetric
    real(dp) :: phi(0:8,0:8,0:8),position(3)
    integer :: axis,ii,jj,kk
    do axis=1,3
      do kk=0,8
      do jj=0,8
      do ii=0,8
        position=real([ii,jj,kk],dp)*.0005_dp
        phi(ii,jj,kk)=position(axis)-.0007_dp
      end do
      end do
      end do
      call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[.004_dp,.004_dp,.004_dp],symmetric,ok)
      call require(ok,'axis-invariant nodal geometry')
      select case(axis)
      case(1)
        call require(all(symmetric%x_face_centroid_y==0).and.all(symmetric%x_face_centroid_z==0).and. &
          all(symmetric%y_face_centroid_z==0).and.all(symmetric%z_face_centroid_y==0),'exact x-plane face symmetry')
      case(2)
        call require(all(symmetric%y_face_centroid_x==0).and.all(symmetric%y_face_centroid_z==0).and. &
          all(symmetric%x_face_centroid_z==0).and.all(symmetric%z_face_centroid_x==0),'exact y-plane face symmetry')
      case(3)
        call require(all(symmetric%z_face_centroid_x==0).and.all(symmetric%z_face_centroid_y==0).and. &
          all(symmetric%x_face_centroid_y==0).and.all(symmetric%y_face_centroid_x==0),'exact z-plane face symmetry')
      end select
    end do
  end block

  block
    type(eb_geometry_3d) :: cut,opposite,saved
    real(dp) :: phi(0:10,0:8,0:6),volume,expected,moment(3),position(3),cell_volume
    integer :: ii,jj,kk
    do kk=0,6
      do jj=0,8
        do ii=0,10
          phi(ii,jj,kk)=real(ii,dp)/10+real(jj,dp)/8+real(kk,dp)/6-2.13_dp
        end do
      end do
    end do
    call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],cut,ok)
    call require(ok.and.cut%is_valid(),'fully 3D oblique plane')
    cell_volume=cut%dx*cut%dy*cut%dz
    expected=.87_dp**3/6._dp
    volume=sum(cut%volume_fraction)*cell_volume
    call assert_close(volume,expected,tolerance,'analytic corner tetrahedron volume')
    call assert_close(sum(cut%boundary_area),sqrt(3._dp)*.87_dp**2/2._dp,tolerance,'analytic oblique facet area')
    moment=0._dp
    do kk=1,6
      do jj=1,8
        do ii=1,10
          position=[(real(ii,dp)-.5_dp+cut%cell_centroid_x(ii,jj,kk))*cut%dx, &
            (real(jj,dp)-.5_dp+cut%cell_centroid_y(ii,jj,kk))*cut%dy, &
            (real(kk,dp)-.5_dp+cut%cell_centroid_z(ii,jj,kk))*cut%dz]
          moment=moment+position*cut%volume_fraction(ii,jj,kk)*cell_volume
        end do
      end do
    end do
    call assert_close(maxval(abs(moment/volume-(1._dp-.87_dp/4._dp))),0._dp,tolerance,'analytic 3D centroid')
    do kk=0,6
      do jj=0,8
        do ii=0,10
          phi(ii,jj,kk)=(real(ii,dp)/10-.5_dp)**2+(real(jj,dp)/8-.5_dp)**2 &
            +(real(kk,dp)/6-.5_dp)**2-.29_dp**2
        end do
      end do
    end do
    call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],cut,ok)
    call require(ok.and.cut%is_valid(),'faceted fully 3D sphere')
    call build_level_set_eb_geometry_3d(-phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],opposite,ok)
    call require(ok.and.opposite%is_valid(),'complementary sphere')
    call assert_close(maxval(abs(cut%volume_fraction+opposite%volume_fraction-1._dp)),0._dp, &
      tolerance,'cellwise complement volume')
    call assert_close(maxval(abs(cut%boundary_area-opposite%boundary_area)),0._dp,tolerance,'complement wall area')
    call assert_close(maxval(abs(cut%boundary_normal_integral_x+opposite%boundary_normal_integral_x)), &
      0._dp,tolerance,'complement x normal')
    call assert_close(maxval(abs(cut%boundary_normal_integral_y+opposite%boundary_normal_integral_y)), &
      0._dp,tolerance,'complement y normal')
    call assert_close(maxval(abs(cut%boundary_normal_integral_z+opposite%boundary_normal_integral_z)), &
      0._dp,tolerance,'complement z normal')
    call assert_close(maximum_geometric_residual(cut),0._dp,tolerance,'sphere discrete closure')
    call require(count(cut%cell_type==eb_covered_cell_3d)>0.and.count(cut%cell_type==eb_cut_cell_3d)>0, &
      'sphere has covered and cut cells')
    saved=cut;phi(0,0,0)=ieee_value(0._dp,ieee_quiet_nan)
    call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],cut,ok)
    call require(.not.ok.and.all(cut%volume_fraction==saved%volume_fraction),'bad 3D level set rollback')
    phi=-1._dp;phi(1,0,0)=1._dp;phi(0,1,0)=1._dp
    call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],cut,ok)
    call require(.not.ok.and.all(cut%volume_fraction==saved%volume_fraction),'disconnected cell rejected atomically')
  end block

  call check_coarsening()

  write(*, '(a)') "test_eb_geometry_3d: PASS"

contains

  subroutine check_coarsening()
    type(eb_geometry_3d) :: fine,coarse,direct,saved
    real(dp) :: phi(0:12,0:12,0:12),lower(3),upper(3),vf,vc
    real(dp), allocatable :: coarse_phi(:,:,:)
    integer :: i,j,k,ratio,n
    logical :: valid
    lower=[-1._dp,2._dp,-3._dp]; upper=[2._dp,4._dp,1._dp]
    do k=0,12
      do j=0,12
        do i=0,12
          phi(i,j,k)=real(i+j+k,dp)/12._dp-1.731_dp
        end do
      end do
    end do
    call build_level_set_eb_geometry_3d(phi,lower,upper,fine,valid)
    call require(valid,'fine oblique coarsening geometry')
    do ratio=2,3
      n=12/ratio
      allocate(coarse_phi(0:n,0:n,0:n))
      do k=0,n
        do j=0,n
          do i=0,n
            coarse_phi(i,j,k)=real(i+j+k,dp)/n-1.731_dp
          end do
        end do
      end do
      call build_level_set_eb_geometry_3d(coarse_phi,lower,upper,direct,valid)
      call require(valid,'independent coarse oblique geometry')
      call coarsen_eb_geometry_3d(fine,ratio,coarse,valid)
      call require(valid.and.coarse%is_valid(),'coarsened oblique geometry')
      call require(all(coarse%cell_type==direct%cell_type),'coarse classification')
      call assert_close(maxval(abs(coarse%volume_fraction-direct%volume_fraction)), &
        0._dp,tolerance,'restricted volume_fraction')
      call assert_close(maxval(abs(coarse%cell_centroid_x-direct%cell_centroid_x)), &
        0._dp,tolerance,'restricted cell_centroid_x')
      call assert_close(maxval(abs(coarse%cell_centroid_y-direct%cell_centroid_y)), &
        0._dp,tolerance,'restricted cell_centroid_y')
      call assert_close(maxval(abs(coarse%cell_centroid_z-direct%cell_centroid_z)), &
        0._dp,tolerance,'restricted cell_centroid_z')
      call assert_close(maxval(abs(coarse%x_face_fraction-direct%x_face_fraction)), &
        0._dp,tolerance,'restricted x_face_fraction')
      call assert_close(maxval(abs(coarse%y_face_fraction-direct%y_face_fraction)), &
        0._dp,tolerance,'restricted y_face_fraction')
      call assert_close(maxval(abs(coarse%z_face_fraction-direct%z_face_fraction)), &
        0._dp,tolerance,'restricted z_face_fraction')
      call assert_close(maxval(abs(coarse%x_face_centroid_y-direct%x_face_centroid_y)), &
        0._dp,tolerance,'restricted x_face_centroid_y')
      call assert_close(maxval(abs(coarse%x_face_centroid_z-direct%x_face_centroid_z)), &
        0._dp,tolerance,'restricted x_face_centroid_z')
      call assert_close(maxval(abs(coarse%y_face_centroid_x-direct%y_face_centroid_x)), &
        0._dp,tolerance,'restricted y_face_centroid_x')
      call assert_close(maxval(abs(coarse%y_face_centroid_z-direct%y_face_centroid_z)), &
        0._dp,tolerance,'restricted y_face_centroid_z')
      call assert_close(maxval(abs(coarse%z_face_centroid_x-direct%z_face_centroid_x)), &
        0._dp,tolerance,'restricted z_face_centroid_x')
      call assert_close(maxval(abs(coarse%z_face_centroid_y-direct%z_face_centroid_y)), &
        0._dp,tolerance,'restricted z_face_centroid_y')
      call assert_close(maxval(abs(coarse%boundary_area-direct%boundary_area)), &
        0._dp,tolerance,'restricted boundary_area')
      call assert_close(maxval(abs(coarse%boundary_centroid_x-direct%boundary_centroid_x)), &
        0._dp,tolerance,'restricted boundary_centroid_x')
      call assert_close(maxval(abs(coarse%boundary_centroid_y-direct%boundary_centroid_y)), &
        0._dp,tolerance,'restricted boundary_centroid_y')
      call assert_close(maxval(abs(coarse%boundary_centroid_z-direct%boundary_centroid_z)), &
        0._dp,tolerance,'restricted boundary_centroid_z')
      call assert_close(maxval(abs(coarse%boundary_normal_x-direct%boundary_normal_x)), &
        0._dp,tolerance,'restricted boundary_normal_x')
      call assert_close(maxval(abs(coarse%boundary_normal_y-direct%boundary_normal_y)), &
        0._dp,tolerance,'restricted boundary_normal_y')
      call assert_close(maxval(abs(coarse%boundary_normal_z-direct%boundary_normal_z)), &
        0._dp,tolerance,'restricted boundary_normal_z')
      call assert_close(maxval(abs(coarse%boundary_normal_integral_x-direct%boundary_normal_integral_x)), &
        0._dp,tolerance,'restricted boundary_normal_integral_x')
      call assert_close(maxval(abs(coarse%boundary_normal_integral_y-direct%boundary_normal_integral_y)), &
        0._dp,tolerance,'restricted boundary_normal_integral_y')
      call assert_close(maxval(abs(coarse%boundary_normal_integral_z-direct%boundary_normal_integral_z)), &
        0._dp,tolerance,'restricted boundary_normal_integral_z')
      deallocate(coarse_phi)
    end do
    saved=coarse
    call coarsen_eb_geometry_3d(fine,0,coarse,valid)
    call require(.not.valid.and.all(coarse%volume_fraction==saved%volume_fraction),'zero ratio rollback')
    call coarsen_eb_geometry_3d(fine,5,coarse,valid)
    call require(.not.valid.and.all(coarse%volume_fraction==saved%volume_fraction),'nondivisible ratio rollback')
    call coarsen_eb_geometry_3d(fine,1,coarse,valid)
    call require(valid.and.all(coarse%volume_fraction==fine%volume_fraction),'ratio one identity')
    do k=0,12
      do j=0,12
        do i=0,12
          phi(i,j,k)=(real(i,dp)/12-.5_dp)**2+(real(j,dp)/12-.5_dp)**2+ &
            (real(k,dp)/12-.5_dp)**2-.29_dp**2
        end do
      end do
    end do
    call build_level_set_eb_geometry_3d(phi,lower,upper,fine,valid)
    call require(valid,'curved fine coarsening geometry')
    call coarsen_eb_geometry_3d(fine,2,coarse,valid)
    call require(valid.and.coarse%is_valid(),'curved coarse validity and face closure')
    vf=fine%dx*fine%dy*fine%dz; vc=coarse%dx*coarse%dy*coarse%dz
    call assert_close(sum(fine%volume_fraction)*vf,sum(coarse%volume_fraction)*vc,tolerance, &
      'curved volume conserved')
    call assert_close(sum(fine%boundary_area),sum(coarse%boundary_area),tolerance,'curved wall area conserved')
    call assert_close(sum(fine%boundary_area*fine%boundary_centroid_x), &
      sum(coarse%boundary_area*coarse%boundary_centroid_x),tolerance,'curved wall x moment conserved')
    call assert_close(sum(fine%boundary_area*fine%boundary_centroid_y), &
      sum(coarse%boundary_area*coarse%boundary_centroid_y),tolerance,'curved wall y moment conserved')
    call assert_close(sum(fine%boundary_area*fine%boundary_centroid_z), &
      sum(coarse%boundary_area*coarse%boundary_centroid_z),tolerance,'curved wall z moment conserved')
    saved=coarse
    call coarsen_eb_geometry_3d(fine,12,coarse,valid)
    call require(.not.valid.and.all(coarse%volume_fraction==saved%volume_fraction), &
      'closed embedded surface without representative normal rejected')
    do i=0,12
      phi(i,:,:)=abs(real(i,dp)/12-.5_dp)-.4_dp
    end do
    call build_level_set_eb_geometry_3d(phi,lower,upper,fine,valid)
    call require(valid,'disconnected coarse block fine geometry')
    call coarsen_eb_geometry_3d(fine,12,coarse,valid)
    call require(.not.valid.and.all(coarse%volume_fraction==saved%volume_fraction), &
      'disconnected coarse fluid regions rejected')
    fine%volume_fraction(1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call coarsen_eb_geometry_3d(fine,2,coarse,valid)
    call require(.not.valid.and.all(coarse%volume_fraction==saved%volume_fraction),'invalid fine geometry rollback')
  end subroutine check_coarsening

  subroutine verify_rotated_plane(axis, expected_cut_cells)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: expected_cut_cells

    real(dp) :: local_volume
    logical :: local_ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.37_dp, geometry, local_ok)
    call require(local_ok .and. geometry%is_valid(), &
      trim(axis)//"-plane geometry validity")
    local_volume = sum(geometry%volume_fraction) * &
      geometry%dx * geometry%dy * geometry%dz
    call assert_close(local_volume, 0.63_dp, tolerance, &
      trim(axis)//"-plane fluid volume")
    call require(count(geometry%cell_type == eb_cut_cell_3d) == &
      expected_cut_cells, trim(axis)//"-plane cut-cell count")
    call assert_close(sum(geometry%boundary_area), 1.0_dp, tolerance, &
      trim(axis)//"-plane boundary area")
    call assert_close(maximum_geometric_residual(geometry), 0.0_dp, &
      tolerance, trim(axis)//"-plane geometric identity")
    select case (axis)
    case ("y")
      call assert_close(maxval(abs(geometry%volume_fraction(:, 3, :) - &
        0.04_dp)), 0.0_dp, tolerance, "y-plane cut fraction")
      call assert_close(maxval(abs(geometry%cell_centroid_y(:, 3, :) - &
        0.48_dp)), 0.0_dp, tolerance, "y-plane cell centroid")
      call assert_close(maxval(abs(geometry%x_face_centroid_y(:, 3, :) - &
        0.48_dp)), 0.0_dp, tolerance, "y-plane x-face centroid")
      call assert_close(maxval(abs(geometry%z_face_centroid_y(:, 3, :) - &
        0.48_dp)), 0.0_dp, tolerance, "y-plane z-face centroid")
      call assert_close(maxval(abs(geometry%boundary_normal_y(:, 3, :) - &
        1.0_dp)), 0.0_dp, tolerance, "y-plane boundary normal")
    case ("z")
      call assert_close(maxval(abs(geometry%volume_fraction(:, :, 3) - &
        0.78_dp)), 0.0_dp, tolerance, "z-plane cut fraction")
      call assert_close(maxval(abs(geometry%cell_centroid_z(:, :, 3) - &
        0.11_dp)), 0.0_dp, tolerance, "z-plane cell centroid")
      call assert_close(maxval(abs(geometry%x_face_centroid_z(:, :, 3) - &
        0.11_dp)), 0.0_dp, tolerance, "z-plane x-face centroid")
      call assert_close(maxval(abs(geometry%y_face_centroid_z(:, :, 3) - &
        0.11_dp)), 0.0_dp, tolerance, "z-plane y-face centroid")
      call assert_close(maxval(abs(geometry%boundary_normal_z(:, :, 3) - &
        1.0_dp)), 0.0_dp, tolerance, "z-plane boundary normal")
    end select
  end subroutine verify_rotated_plane

  real(dp) function maximum_geometric_residual(local_geometry) result(error)
    type(eb_geometry_3d), intent(in) :: local_geometry

    real(dp) :: residual(3)
    integer :: i, j, k

    error = 0.0_dp
    do k = 1, local_geometry%nz
      do j = 1, local_geometry%ny
        do i = 1, local_geometry%nx
          residual = [ &
            local_geometry%dy * local_geometry%dz * ( &
              local_geometry%x_face_fraction(i, j, k) - &
              local_geometry%x_face_fraction(i - 1, j, k)) - &
                local_geometry%boundary_normal_integral_x(i, j, k), &
            local_geometry%dx * local_geometry%dz * ( &
              local_geometry%y_face_fraction(i, j, k) - &
              local_geometry%y_face_fraction(i, j - 1, k)) - &
                local_geometry%boundary_normal_integral_y(i, j, k), &
            local_geometry%dx * local_geometry%dy * ( &
              local_geometry%z_face_fraction(i, j, k) - &
              local_geometry%z_face_fraction(i, j, k - 1)) - &
                local_geometry%boundary_normal_integral_z(i, j, k)]
          error = max(error, maxval(abs(residual)))
        end do
      end do
    end do
  end function maximum_geometric_residual

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

end program test_eb_geometry_3d
