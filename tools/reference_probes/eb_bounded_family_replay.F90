! Saved prepared-state replay; geometry rebuilt from the original plane input.
program eb_bounded_family_replay
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
  use lidryer_mechanism_mod, only: load_lidryer_thermo,lidryer_nspecies
  use reactive_1d_mod, only: reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_guarded_quartic_mod, only: guarded_children
  implicit none
  integer, parameter :: nv=5+lidryer_nspecies
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geom
  real(dp) :: q(nv+1,-2:2,-2:2,-2:2),children(nv+1,8),vf(8),weights(8)
  real(dp) :: dx(3),plane,theta,t,cs,prim(nv),residual
  integer :: root,src(3),ios,i,j,k,n,c,eos
  logical :: ok,smooth
  read(*,*) root,dx,plane,eos
  call load_lidryer_thermo(species,ok)
  if(.not.ok) error stop 'thermo load'
  call configure_nasa7_eos_policy(species,.true.,ok,converged_inverse=eos==1)
  if(.not.ok) error stop 'EOS policy'
  call build_axis_plane_eb_geometry_3d(2*root,2*root,2*root,0._dp,root*dx(1), &
    0._dp,root*dx(2),0._dp,root*dx(3),'z',plane,geom,ok)
  if(.not.ok) error stop 'geometry'
  do
    read(*,*,iostat=ios) src
    if(ios<0) exit
    if(ios/=0) error stop 'source read'
    do k=-2,2
    do j=-2,2
    do i=-2,2
      read(*,*) q(:,i,j,k)
    enddo
    enddo
    enddo
    do k=0,1
    do j=0,1
    do i=0,1
      n=1+i+2*j+4*k
      vf(n)=geom%volume_fraction(2*src(1)-1+i,2*src(2)-1+j,2*src(3)-1+k)
    enddo
    enddo
    enddo
    call reactive_conserved_to_primitive(species,q(:nv,0,0,0),900._dp,prim,t,cs,ok, &
      internal_energy_density=q(nv+1,0,0,0))
    if(.not.ok) error stop 'parent held EOS'
    children=-777;theta=-1
    call guarded_children(species,q,t,children,theta,smooth,ok,vf,bounded_correction=.true.)
    if(.not.ok) error stop 'family admission'
    weights=vf/sum(vf);residual=0
    do c=1,nv+1
      residual=max(residual,abs(sum(weights*children(c,:))-q(c,0,0,0))/max(1._dp,abs(q(c,0,0,0))))
    enddo
    if(residual>1e-12_dp) error stop 'weighted mean'
    write(*,'(3(i0,1x),11(es25.16e3,1x))') src,theta,residual,t,vf
    do n=1,8
      write(*,'(*(es25.16e3,1x))') children(:,n)
    enddo
  enddo
end program
