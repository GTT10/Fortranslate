program test_amr_regrid_3d
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d,prolong_pcm_3d,composite_integrals_amr_3d
  use amr_reactive_3d_mod, only: regrid_amr_reactive_3d,temperature_tag_patch_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config
  use simulation_config_amr_reactive_3d_mod, only: amr_reactive_3d_config,validate_amr_reactive_3d_configuration
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patch,next_patch,bad
  type(reactive_3d_config) :: base
  type(amr_reactive_3d_config) :: config
  real(dp) :: c(12,8,8,8),t(8,8,8),before_c(12,8,8,8),q(12),sound,expected,tag(8,8,8)
  real(dp), allocatable :: f(:,:,:,:),ft(:,:,:),before(:,:,:,:)
  integer :: i,j,k
  logical :: ok,changed
  character(len=256) :: message
  character(len=32) :: faces(6)
  real(dp) :: initial(12),final(12),inlet(12,6),inlet_t(6)
  call load_h2o2_elementary_thermo(species,ok); call require(ok,'thermo')
  call initialize_amr_patch_3d(8,8,8,2,4,2,4,2,4,2,patch,ok); call require(ok,'old patch')
  do k=1,8
    do j=1,8
      do i=1,8
        q=0; q(1)=.3_dp+.01_dp*i+.005_dp*j+.002_dp*k; q(5)=101325*q(1)/.3_dp
        q(2)=20; q(6)=.03_dp; q(9)=.22_dp; q(12)=.75_dp
        call reactive_primitive_to_conserved(species,q,c(:,i,j,k),t(i,j,k),sound,ok)
        call require(ok,'initial cells')
      end do
    end do
  end do
  allocate(f(12,6,6,6),ft(6,6,6))
  call prolong_pcm_3d(c,patch,f,ok); call require(ok,'initial prolongation')
  ft=t(1,1,1)
  do i=1,6
    f(:,i,:,:)=f(:,i,:,:)*(1.0_dp+.01_dp*(-1.0_dp)**i)
  end do
  before=f
  call initialize_amr_patch_3d(8,8,8,3,6,3,6,3,6,2,next_patch,ok); call require(ok,'next patch')
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'characteristic_plm','mc')
  call require(ok,'conservative expand and move')
  call require(all(f(:,1:4,1:4,1:4)==before(:,3:6,3:6,3:6)),'exact overlap')
  expected=.3_dp+.01_dp*(5-.25_dp)+.005_dp*(5-.25_dp)+.002_dp*(5-.25_dp)
  call require(abs(f(1,5,5,5)-expected)<1e-13_dp,'linear prolongation')
  before=f; before_c=c; bad=patch; bad%coarse_i_lower=0
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,bad,ok,'pcm','mc')
  call require(.not.ok.and.all(f==before).and.all(c==before_c),'rollback')
  call initialize_amr_patch_3d(8,8,8,2,3,2,3,2,3,2,next_patch,ok)
  call require(ok,'shrink patch')
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'pcm','mc')
  call require(ok.and.size(f,2)==4,'conservative shrink')
  tag=1000; tag(6,5,4)=1200
  call temperature_tag_patch_3d(patch,tag,1100.0_dp,1,next_patch,changed,ok)
  call require(ok.and.changed.and.next_patch%coarse_i_lower==5.and.next_patch%coarse_i_upper==7,'tag bounds')
  call temperature_tag_patch_3d(patch,tag,1300.0_dp,1,next_patch,changed,ok)
  call require(ok.and..not.changed,'empty tags retain current patch')
  base%nx=8; base%ny=8; base%nz=8; base%boundary_faces(1:2)='outflow'
  call validate_amr_reactive_3d_configuration(base,config,ok,message)
  call require(.not.ok,'physical face override rejected')
  base%chemistry_enabled=.false.
  base%boundary_faces=''; config%regrid_interval_steps=1; config%checkpoint_interval_steps=1
  call validate_amr_reactive_3d_configuration(base,config,ok,message)
  call require(ok,'dynamic checkpoint supported')
  faces='outflow';faces(1)='slip_wall';faces(5)='no_slip_wall'
  call composite_integrals_amr_3d(c,f,patch,1._dp,1._dp,1._dp,initial,ok)
  call require(ok,'physical initial inventory')
  before=f
  call initialize_amr_patch_3d(8,8,8,1,3,1,3,1,3,2,next_patch,ok)
  call require(ok,'physical lower corner patch')
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'characteristic_plm','mc',faces)
  call require(ok,'physical lower corner growth')
  call require(all(f(:,3:6,3:6,3:6)==before),'physical overlap preserved exactly')
  call composite_integrals_amr_3d(c,f,patch,1._dp,1._dp,1._dp,final,ok)
  call require(ok.and.all(abs(final-initial)<=5e-12_dp*max(1._dp,abs(initial))),'physical grow inventory')
  before=f;before_c=c
  faces(1)='inflow';inlet=spread(c(:,1,1,1),2,6);inlet_t=t(1,1,1)
  inlet(1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'characteristic_plm','mc',faces,inlet,inlet_t)
  call require(.not.ok.and.all(f==before).and.all(c==before_c),'physical invalid reservoir rollback')
  faces(1)='slip_wall'
  call initialize_amr_patch_3d(8,8,8,2,3,2,3,2,3,2,next_patch,ok)
  call require(ok,'physical shrink target')
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'pcm','mc',faces)
  call require(ok.and.all(f==before(:,3:6,3:6,3:6)),'physical shrink overlap')
  call initialize_amr_patch_3d(8,8,8,6,8,6,8,6,8,2,next_patch,ok)
  call require(ok,'physical upper corner target')
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'characteristic_plm','mc',faces)
  call require(ok,'physical upper corner move')
  call composite_integrals_amr_3d(c,f,patch,1._dp,1._dp,1._dp,final,ok)
  call require(ok.and.all(abs(final-initial)<=5e-12_dp*max(1._dp,abs(initial))),'physical move inventory')
  tag=1000;tag(1,8,1)=1200
  call temperature_tag_patch_3d(patch,tag,1100._dp,1,next_patch,changed,ok,faces)
  call require(ok.and.changed.and.next_patch%coarse_i_lower==1.and.next_patch%coarse_j_upper==8.and. &
    next_patch%coarse_k_lower==1,'physical boundary tags included')
  faces='periodic'
  call regrid_amr_reactive_3d(species,patch,c,t,f,ft,next_patch,ok,'pcm','mc',faces)
  call require(.not.ok,'periodic seam remap rejected')
  print *, 'PASS: conservative remap, linear prolongation, exact overlap, shrink, tags and rollback'
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *,label
      error stop 1
    end if
  end subroutine
end program
