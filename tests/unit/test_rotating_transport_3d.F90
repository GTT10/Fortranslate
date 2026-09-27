program test_rotating_transport_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    compute_sdc_eb_tree_transport_3d,compute_eb_tree_wall_transport_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: boxes(0)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: g,extended(1)
  type(eb_reactive_tree_3d) :: ordinary,rotated
  type(amr_tree_state_3d) :: a,b
  type(eb_integrated_flux_3d), allocatable :: fa(:),fb(:)
  type(rotating_frame_3d) :: frame
  real(dp), allocatable :: q(:),u(:,:,:,:),t(:,:,:),exchange(:)
  real(dp) :: sound,potential,position(3),ta,tb,error,worst,frequency
  integer :: nv,i,j,k,s,mode
  logical :: ok,flags(4)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  nv=reactive_nvar(size(species));allocate(q(reactive_nprim(size(species))),u(nv,4,4,4),t(4,4,4))
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],boxes,topology,ok);call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-10._dp,g,ok)
  call require(ok,'geometry')
  do k=1,4
  do j=1,4
  do i=1,4
    q=0;q(1)=1+.01_dp*i;q(2:4)=[.1_dp*i,.2_dp*j,-.1_dp*k];q(5)=3e5_dp+1000*j
    do s=1,size(species)
      select case(species(s)%name)
      case('H2');q(5+s)=.02_dp+.001_dp*i
      case('O2');q(5+s)=.2_dp
      case('N2');q(5+s)=.78_dp-.001_dp*i
      end select
    end do
    call reactive_primitive_to_conserved(species,q,u(:,i,j,k),t(i,j,k),sound,ok);call require(ok,'initial state')
  end do
  end do
  end do
  call initialize_reactive_eb_tree_3d(species,topology,g,u,t,ordinary,ok);call require(ok,'ordinary tree')
  rotated=ordinary;frame%angular_speed=100;frame%axis_point=[10._dp,10._dp,0._dp]
  do k=1,4
  do j=1,4
  do i=1,4
    position=(real([i,j,k],dp)-.5_dp)/4
    call rotating_potential_3d(frame,position,potential,ok);call require(ok,'potential')
    rotated%nodes(1)%state(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
  end do
  end do
  end do
  call require(maxval(rotated%nodes(1)%state(5,:,:,:))<0,'negative rotating energy')
  worst=0
  do mode=1,4
    flags=.false.;flags(min(mode,3))=.true.;flags(4)=mode==4
    call compute_sdc_eb_tree_transport_3d(species,transport,1e-9_dp,flags(1),flags(2),flags(3),flags(4), &
      ordinary,a,fa,ta,ok);call require(ok,'ordinary transport')
    call compute_sdc_eb_tree_transport_3d(species,transport,1e-9_dp,flags(1),flags(2),flags(3),flags(4), &
      rotated,b,fb,tb,ok,frame=frame);call require(ok,'rotating transport')
    error=maxval(abs(a%fields(1)%values-b%fields(1)%values)/max(1._dp,abs(a%fields(1)%values)))
    error=max(error,maxval(abs(fa(1)%x-fb(1)%x)),maxval(abs(fa(1)%y-fb(1)%y)),maxval(abs(fa(1)%z-fb(1)%z)))
    call require(error<3e-11_dp.and.ta==tb,'transport energy convention invariance')
    call require(maxval(abs(a%fields(1)%values))>1e-8_dp,'nonzero transport fixture')
    call require(all(b%fields(1)%values(1,:,:,:)==0),'zero diffusive total mass rate')
    worst=max(worst,error)
  end do
  call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.31_dp,g,ok)
  call require(ok,'cut geometry')
  call build_axis_plane_eb_geometry_3d(8,8,8,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp, &
    'x',.31_dp,extended(1),ok);call require(ok,'wall halo geometry')
  call initialize_reactive_eb_tree_3d(species,topology,g,u,t,ordinary,ok);call require(ok,'cut tree')
  rotated=ordinary
  do k=1,4
  do j=1,4
  do i=1,4
    if(g%volume_fraction(i,j,k)<=0) cycle
    position=(real([i,j,k],dp)-.5_dp)/4
    call rotating_potential_3d(frame,position,potential,ok);call require(ok,'wall potential')
    rotated%nodes(1)%state(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
  end do
  end do
  end do
  allocate(exchange(nv))
  call compute_eb_tree_wall_transport_3d(species,transport,rotated,extended,.true.,.true.,400._dp, &
    b,exchange,frequency,ok,frame=frame)
  call require(ok,'rotating wall heat and traction')
  call require(exchange(5)<0.and.exchange(3)<0.and.frequency>0,'wall heat and drag direction')
  error=0
  do s=1,nv
    error=max(error,abs(exchange(s)-g%dx*g%dy*g%dz*sum(g%volume_fraction*b%fields(1)%values(s,:,:,:))))
  end do
  call require(error<1e-12_dp,'wall source integral')
  print *, 'Wall energy rate, y momentum rate, inventory residual:',exchange(5),exchange(3),error
  print *, 'PASS: rotating viscosity, conduction, diffusion and barodiffusion; maximum error:',worst
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
