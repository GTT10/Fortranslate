program eb_guarded_quartic_probe
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use eb_guarded_quartic_mod, only: guarded_children,admit_family,curvature_screen
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_quartic_stencil_geometry_mod, only: prepare_connected_quartic_stencil
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: neighborhood,fine_geometry,parent_geometry
  real(dp), allocatable :: prepared(:,:,:,:),saved_prepared(:,:,:,:)
  real(dp) :: volumes(8),weights(8),physical_total,rotated_volumes(8)
  logical :: available(-2:2,-2:2,-2:2),axes(3),saved_axes(3)
  real(dp), allocatable :: stencil(:,:,:,:),children(:,:),base(:,:),high(:,:),gas(:),q(:),parent(:),saved(:,:)
  real(dp) :: theta,t,cs,h,width,angle,mean_u,mean_u2,exact,error,delta,saved_theta
  integer :: nv,n2,c,i,j,k,n,other,checks,grid,phase,iw
  real(dp) :: transition_width,transition_offset
  logical :: ok,smooth
  checks=0
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species));n2=0
  do c=1,size(species)
    if(species(c)%name=='N2') n2=c
  enddo
  call require(n2>0,'nitrogen')
  allocate(stencil(nv+1,-2:2,-2:2,-2:2),children(nv+1,8),base(nv+1,8),high(nv+1,8), &
    gas(nv),q(reactive_nprim(size(species))),parent(nv+1),saved(nv+1,8), &
    prepared(nv+1,-2:2,-2:2,-2:2),saved_prepared(nv+1,-2:2,-2:2,-2:2))
  q=0;q(1)=1;q(5)=universal_gas_constant/species(n2)%molecular_weight*900;q(5+n2)=1
  call reactive_primitive_to_conserved(species,q,gas,t,cs,ok)
  call require(ok,'initial EOS')
  h=gas(5);parent(:nv)=gas;parent(nv+1)=h
  width=2*acos(-1._dp)/16
  do k=-2,2
  do j=-2,2
  do i=-2,2
    angle=(4+j-.5_dp)*width
    mean_u=10*sin(angle)*sin(width/2)/(width/2)
    mean_u2=50*(1-cos(2*angle)*sin(width)/width)
    stencil(:,i,j,k)=parent
    stencil(2,i,j,k)=mean_u
    stencil(5,i,j,k)=h+.5_dp*mean_u2
  enddo
  enddo
  enddo
  theta=-1;children=-777
  call guarded_children(species,stencil,900._dp,children,theta,smooth,ok)
  call require(ok.and.smooth.and.theta==1,'smooth shear retains quartic')
  exact=10*sin(7.5_dp*width/2)*sin(width/4)/(width/4)
  error=abs(children(2,3)/children(1,3)-exact)
  call require(error<.000174_dp,'shear error reduced')
  do c=1,nv+1
    call require(abs(sum(children(c,:))/8-stencil(c,0,0,0))<1e-12_dp*max(1._dp,abs(stencil(c,0,0,0))), &
      'smooth conserved family mean')
  enddo
  print *, 'shear_velocity_error_m_s=',error
  ! A near-zero transverse component must not set its own cancellation scale.
  saved_prepared=stencil
  stencil(4,0:2,:,:)=5.82e-16_dp
  call require(.not.curvature_screen(stencil),'scalar scale reproduces spurious transverse rejection')
  call require(curvature_screen(stencil,.true.),'momentum vector scale admits roundoff transverse perturbation')
  theta=-1
  call guarded_children(species,stencil,900._dp,children,theta,smooth,ok,momentum_roundoff_scale=.true.)
  call require(ok.and.smooth.and.theta==1,'roundoff-scale candidate retains EOS-admissible family')
  stencil=stencil*1.e-6_dp
  call require(curvature_screen(stencil,.true.),'vector screen invariant under conserved-unit scaling')
  stencil=saved_prepared
  stencil(4,0:2,:,:)=1.e-6_dp
  call require(.not.curvature_screen(stencil,.true.),'resolved transverse jump remains rejected')
  stencil=saved_prepared
  stencil(1,0:2,:,:)=1.01_dp
  call require(.not.curvature_screen(stencil,.true.),'density jump remains rejected')
  stencil=saved_prepared
  ! A rigid z-axis rotation preserves physical volume weights and admissibility.
  stencil(4,0:2,:,:)=5.82e-16_dp
  volumes(:4)=.44_dp;volumes(5:)=1._dp
  call guarded_children(species,stencil,900._dp,base,theta,smooth,ok,volumes,.true.)
  call require(ok.and.smooth.and.theta==1,'weighted roundoff family before rotation')
  do k=-2,2
  do j=-2,2
  do i=-2,2
    prepared(:,i,j,k)=stencil(:,j,-i,k)
    prepared(2,i,j,k)=-stencil(3,j,-i,k)
    prepared(3,i,j,k)=stencil(2,j,-i,k)
  enddo
  enddo
  enddo
  call guarded_children(species,prepared,900._dp,high,theta,smooth,ok,volumes,.true.)
  call require(ok.and.smooth.and.theta==1,'rotated weighted roundoff family admitted')
  do k=0,1
  do j=0,1
  do i=0,1
    n=1+i+2*j+4*k;other=1+j+2*(1-i)+4*k
    children(:,n)=base(:,other)
    children(2,n)=-base(3,other);children(3,n)=base(2,other)
    call require(maxval(abs(high(:,n)-children(:,n))/max(1._dp,abs(children(:,n))))< &
      128*epsilon(1._dp),'weighted quartic family commutes with axis rotation')
  enddo
  enddo
  enddo
  stencil=saved_prepared


  call build_axis_plane_eb_geometry_3d(5,5,5,-2._dp,3._dp,-2._dp,3._dp,-2._dp,3._dp, &
    'z',.28_dp,neighborhood,ok)
  call require(ok,'coarse stencil geometry')
  call build_axis_plane_eb_geometry_3d(2,2,2,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    'z',.28_dp,fine_geometry,ok)
  call require(ok,'fine child geometry')
  call coarsen_eb_geometry_3d(fine_geometry,2,parent_geometry,ok)
  call require(ok,'physical parent geometry')
  volumes=reshape(fine_geometry%volume_fraction,[8])*fine_geometry%dx*fine_geometry%dy*fine_geometry%dz
  physical_total=parent_geometry%volume_fraction(1,1,1)*parent_geometry%dx*parent_geometry%dy*parent_geometry%dz
  call require(abs(sum(volumes)-physical_total)<1e-14_dp,'physical parent child volume sum')
  available=.true.
  stencil(:,:,:,-2:-1)=ieee_value(0._dp,ieee_quiet_nan)
  call prepare_connected_quartic_stencil(stencil,neighborhood,available,prepared,axes,ok)
  call require(ok.and.all(axes.eqv.[.true.,.true.,.false.]),'wall retains connected tangent plane')
  call require(all(prepared(:,:,:,-2)==prepared(:,:,:,2)),'blocked normal direction collapsed without solid reads')
  call guarded_children(species,prepared,900._dp,children,theta,smooth,ok,volumes)
  call require(ok.and.smooth.and.theta==1,'weighted cut family admits smooth shear')
  call require(abs(children(2,3)/children(1,3)-exact)<.000174_dp,'weighted cut shear retains accuracy')
  do c=1,nv+1
    call require(abs(sum(volumes*children(c,:))-physical_total*prepared(c,0,0,0))< &
      1e-12_dp*max(1._dp,abs(physical_total*prepared(c,0,0,0))), 'physical volume weighted inventory')
  enddo
  print *, 'physical_parent_volume=',physical_total
  ! Nonuniform algebraic weights exercise recentering and covered children;
  ! this is not a general cut-cell polynomial-accuracy claim.
  weights=[0._dp,.1_dp,.2_dp,.3_dp,.4_dp,.5_dp,.6_dp,1._dp]
  call guarded_children(species,prepared,900._dp,children,theta,smooth,ok,weights)
  call require(ok,'unequal weights with a covered child')
  do c=1,nv+1
    call require(abs(sum(weights*children(c,:))/sum(weights)-prepared(c,0,0,0))< &
      1e-12_dp*max(1._dp,abs(prepared(c,0,0,0))), 'unequal weighted mean')
  enddo
  saved=children;saved_theta=theta;weights=-1
  call guarded_children(species,prepared,900._dp,children,theta,smooth,ok,weights)
  call require(.not.ok.and.all(children==saved).and.theta==saved_theta,'negative weights rollback')
  weights=0
  call guarded_children(species,prepared,900._dp,children,theta,smooth,ok,weights)
  call require(.not.ok.and.all(children==saved),'fully covered family rollback')
  weights=1;weights(1)=ieee_value(0._dp,ieee_quiet_nan)
  call guarded_children(species,prepared,900._dp,children,theta,smooth,ok,weights)
  call require(.not.ok.and.all(children==saved),'nonfinite weights rollback')
  available(1,1,0)=.false.
  call prepare_connected_quartic_stencil(stencil,neighborhood,available,prepared,axes,ok)
  call require(ok.and..not.any(axes),'missing retained diagonal selects PCM')
  call require(all(prepared(2,:,:,:)==stencil(2,0,0,0)),'PCM fallback preserves center')
  available=.true.;saved_prepared=prepared;saved_axes=axes
  stencil(1,1,1,0)=ieee_value(0._dp,ieee_quiet_nan)
  call prepare_connected_quartic_stencil(stencil,neighborhood,available,prepared,axes,ok)
  call require(.not.ok.and.all(prepared==saved_prepared).and.all(axes.eqv.saved_axes),'invalid active stencil rollback')

  ! A velocity jump is thermodynamically valid but should fail smooth screening.
  do k=-2,2
  do j=-2,2
  do i=-2,2
    stencil(:,i,j,k)=parent
    if(j>0) then
      stencil(2,i,j,k)=10;stencil(5,i,j,k)=h+50
    endif
  enddo
  enddo
  enddo
  call guarded_children(species,stencil,900._dp,children,theta,smooth,ok)
  call require(ok.and..not.smooth.and.all(children(2,:)==0),'jump returns limited baseline')
  ! Shared contraction must preserve all parent means, including energy.
  base=spread(parent,2,8);high=base
  do n=1,8
    delta=merge(-1._dp,1._dp,mod(n,2)==0)
    high(2,n)=delta*1.e6_dp
  enddo
  call admit_family(species,base,high,900._dp,children,theta,ok)
  call require(ok.and.theta>0.and.theta<1,'kinetic energy contracts entire family')
  do c=1,nv+1
    call require(abs(sum(children(c,:))/8-parent(c))<1e-12_dp*max(1._dp,abs(parent(c))), 'contracted mean')
  enddo
  print *, 'kinetic_contraction=',theta
  high=base
  do n=1,8
    delta=merge(-1._dp,1._dp,mod(n,2)==0)*10*abs(h)
    high(5,n)=high(5,n)+delta;high(nv+1,n)=high(nv+1,n)+delta
  enddo
  call admit_family(species,base,high,900._dp,children,theta,ok)
  call require(ok.and.theta>0.and.theta<1,'thermal EOS contraction')
  print *, 'thermal_contraction=',theta
  high=base
  do n=1,8
    delta=merge(-.01_dp,.01_dp,mod(n,2)==0)
    high(6,n)=delta;high(5+n2,n)=high(5+n2,n)-delta
  enddo
  call admit_family(species,base,high,900._dp,children,theta,ok)
  call require(ok.and.theta==0.and.all(children==base),'absent species forces exact baseline')
  saved=children;saved_theta=theta;high(1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call admit_family(species,base,high,900._dp,children,theta,ok)
  call require(.not.ok.and.all(children==saved).and.theta==saved_theta,'nonfinite rollback')
  high=base;high(2,1)=high(2,1)+1
  call admit_family(species,base,high,900._dp,children,theta,ok)
  call require(.not.ok.and.all(children==saved),'nonconservative input rejected')
  ! Smooth inflections, unlike extrema, have resolved curvature sign changes.
  do grid=16,32,16
    width=2*acos(-1._dp)/grid
    do k=-2,2
    do j=-2,2
    do i=-2,2
      angle=(1+j-.5_dp)*width
      stencil(:,i,j,k)=parent
      stencil(2,i,j,k)=10*sin(angle)*sin(width/2)/(width/2)
      stencil(5,i,j,k)=h+25*(1-cos(2*angle)*sin(width)/width)
    enddo
    enddo
    enddo
    call require(.not.curvature_screen(stencil,.true.),'default inflection fallback preserved')
    call require(curvature_screen(stencil,.true.,.true.),'resolved shear inflection admitted')
    theta=-1
    call guarded_children(species,stencil,900._dp,children,theta,smooth,ok, &
      momentum_roundoff_scale=.true.,inflection_hierarchy=.true.)
    call require(ok.and.smooth.and.theta==1,'inflection conserved EOS family admitted')
    do c=1,nv+1
      call require(abs(sum(children(c,:))/8-stencil(c,0,0,0))<1e-12_dp*max(1._dp, &
        abs(stencil(c,0,0,0))),'inflection family conserves each component')
    enddo
    call require(curvature_screen(stencil*1.e-6_dp,.true.,.true.),'inflection unit-scale invariant')
  enddo
  ! Finite-volume tanh transitions: four widths, nine subcell positions.
  do iw=1,4
    transition_width=2._dp**(iw-4)
    do phase=-4,4
      transition_offset=real(phase,dp)/8
      do j=-2,2
        stencil(:, :,j,:)=spread(spread(parent,2,5),3,5)
        stencil(2,:,j,:)=transition_width*( &
          log(cosh((j+.5_dp-transition_offset)/transition_width))- &
          log(cosh((j-.5_dp-transition_offset)/transition_width)))
      enddo
      call require(.not.curvature_screen(stencil,.true.,.true.),'underresolved tanh rejected')
    enddo
  enddo
  ! Bounded tensor family at the residual y=3/4 band.
  width=2*acos(-1._dp)/16
  do k=-2,2
  do j=-2,2
  do i=-2,2
    angle=(2+j-.5_dp)*width
    stencil(:,i,j,k)=parent
    stencil(2,i,j,k)=10*sin(angle)*sin(width/2)/(width/2)
    stencil(5,i,j,k)=h+25*(1-cos(2*angle)*sin(width)/width)
  enddo
  enddo
  enddo
  theta=-1
  call guarded_children(species,stencil,900._dp,children,theta,smooth,ok,bounded_correction=.true.)
  call require(ok.and.theta==1,'bounded tensor gas EOS admission')
  exact=10*sin(2.5_dp*width/2)*sin(width/4)/(width/4)
  call require(abs(children(2,1)-exact)<.000741_dp,'bounded tensor residual-band improvement')
  volumes=[.2_dp,.3_dp,.4_dp,.5_dp,.6_dp,.7_dp,.8_dp,1._dp]
  call guarded_children(species,stencil,900._dp,children,theta,smooth,ok,volumes,bounded_correction=.true.)
  call require(ok,'bounded tensor weighted EOS admission')
  weights=volumes/sum(volumes)
  do c=1,nv+1
    call require(abs(sum(weights*children(c,:))-stencil(c,0,0,0))< &
      1e-12_dp*max(1._dp,abs(stencil(c,0,0,0))),'bounded tensor weighted mean')
  enddo
  base=children
  do k=-2,2
  do j=-2,2
  do i=-2,2
    prepared(:,i,j,k)=stencil(:,j,-i,k)
    prepared(2,i,j,k)=-stencil(3,j,-i,k)
    prepared(3,i,j,k)=stencil(2,j,-i,k)
  enddo
  enddo
  enddo
  do k=0,1
  do j=0,1
  do i=0,1
    n=1+i+2*j+4*k;other=1+j+2*(1-i)+4*k
    rotated_volumes(n)=volumes(other)
  enddo
  enddo
  enddo
  call guarded_children(species,prepared,900._dp,high,theta,smooth,ok, &
    rotated_volumes,bounded_correction=.true.)
  call require(ok.and.theta==1,'bounded weighted rotated gas EOS')
  do k=0,1
  do j=0,1
  do i=0,1
    n=1+i+2*j+4*k;other=1+j+2*(1-i)+4*k
    children(:,n)=base(:,other)
    children(2,n)=-base(3,other);children(3,n)=base(2,other)
    call require(maxval(abs(high(:,n)-children(:,n))/max(1._dp,abs(children(:,n))))< &
      128*epsilon(1._dp),'bounded weighted gas rotation')
  enddo
  enddo
  enddo
  saved=children;saved_theta=theta
  stencil(1,0,0,0)=-1
  call guarded_children(species,stencil,900._dp,children,theta,smooth,ok,volumes,bounded_correction=.true.)
  call require(.not.ok.and.all(children==saved).and.theta==saved_theta,'bounded tensor failure rollback')
  print *, 'PASS guarded candidate checks=',checks
contains
  subroutine require(test,label)
    logical, intent(in) :: test
    character(len=*), intent(in) :: label
    if(.not.test) then
      print *, 'FAIL ',label
      error stop 1
    endif
    checks=checks+1
  end subroutine
end program
