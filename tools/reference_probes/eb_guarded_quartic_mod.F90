! Experimental volume-weighted family safeguard, not wired into EB halo filling.
! Curvature-sign screening is a diagnostic heuristic, not a proven shock sensor.
module eb_guarded_quartic_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use eb_reactive_transfer_3d_mod, only: interpolate_eb_linear_halo_3d
  use eb_quartic_kernel_mod, only: child_average
  use eb_bounded_tensor_mod, only: bounded_child_average
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: guarded_children,admit_family,curvature_screen
contains
  logical function curvature_screen(q,momentum_roundoff_scale,inflection_hierarchy) result(smooth)
    real(dp), intent(in) :: q(:,-2:,-2:,-2:)
    logical, optional, intent(in) :: momentum_roundoff_scale
    logical, optional, intent(in) :: inflection_hierarchy
    logical :: vector_scale,allow_inflection
    real(dp) :: line(-2:2),d2(3),noise,line_scale
    integer :: c,a,i,j,k,v,idx(3)
    smooth=.false.
    if(any(shape(q(1,:,:,:))/=[5,5,5]).or..not.all(ieee_is_finite(q))) return
    allow_inflection=.false.
    if(present(inflection_hierarchy)) allow_inflection=inflection_hierarchy
    vector_scale=.false.
    if(present(momentum_roundoff_scale)) vector_scale=momentum_roundoff_scale
    if(vector_scale.and.size(q,1)<4) return
    do c=1,size(q,1)
    do a=1,3
    do j=-2,2
    do i=-2,2
      line_scale=0
      do v=-2,2
        select case(a)
        case(1);idx=[v,i,j]
        case(2);idx=[i,v,j]
        case(3);idx=[i,j,v]
        end select
        line(v)=q(c,idx(1),idx(2),idx(3))
        if(vector_scale.and.c>=2.and.c<=4) &
          line_scale=max(line_scale,norm2(q(2:4,idx(1),idx(2),idx(3))))
      enddo
      do k=1,3
        d2(k)=(line(k-1)-line(k-2))-(line(k-2)-line(k-3))
      enddo
      ! Only a floating-point cancellation allowance, no fitted physical scale.
      noise=16*epsilon(1._dp)*max(maxval(abs(line)),line_scale)
      if(all(abs(d2)<=noise)) cycle
      if(all(d2>noise).or.all(d2< -noise)) cycle
      if(.not.allow_inflection) return
      if(.not.monotone_inflection(line,noise)) return
    enddo
    enddo
    enddo
    enddo
    smooth=.true.
  end function

  ! Experimental difference hierarchy, not a proven discontinuity detector.
  pure logical function monotone_inflection(line,noise) result(admit)
    real(dp), intent(in) :: line(-2:2),noise
    real(dp) :: d1(4),d2(3),d3(2),low,high
    d1=line(-1:2)-line(-2:1);d2=d1(2:4)-d1(1:3);d3=d2(2:3)-d2(1:2)
    admit=.false.
    if(.not.(all(d1>noise).or.all(d1< -noise))) return
    if(maxval(abs(d3))>minval(abs(d1))+noise) return
    low=(-3*line(-2)+22*line(-1)+128*line(0)-22*line(1)+3*line(2))/128
    high=2*line(0)-low
    if(low<min(line(-1),line(0))-noise.or.low>max(line(-1),line(0))+noise) return
    if(high<min(line(0),line(1))-noise.or.high>max(line(0),line(1))+noise) return
    admit=.true.
  end function

  subroutine guarded_children(species,q,temperature,children,theta,smooth,ok,child_volume,momentum_roundoff_scale, &
      inflection_hierarchy,bounded_correction,continuous_correction)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: q(:,-2:,-2:,-2:),temperature
    real(dp), intent(inout) :: children(:,:),theta
    logical, intent(out) :: smooth,ok
    logical, optional, intent(in) :: momentum_roundoff_scale
    real(dp), optional, intent(in) :: child_volume(8)
    real(dp) :: weights(8)
    real(dp) :: base(size(q,1),8),high(size(q,1),8),neighbors(size(q,1),6),value(size(q,1)),offset(3)
    integer :: i,j,k,c,n,nv
    logical :: valid,bounded,continuous
    logical, optional, intent(in) :: inflection_hierarchy,bounded_correction,continuous_correction
    ok=.false.;smooth=.false.;nv=reactive_nvar(size(species))
    bounded=.false.;continuous=.false.
    if(present(bounded_correction)) bounded=bounded_correction
    if(present(continuous_correction)) continuous=continuous_correction
    bounded=bounded.or.continuous
    if(size(q,1)/=nv+1.or.any(shape(q(1,:,:,:))/=[5,5,5])) return
    if(any(shape(children)/=[nv+1,8]).or..not.all(ieee_is_finite(q))) return
    call normalize_weights(child_volume,weights,valid)
    if(.not.valid) return
    neighbors(:,1)=q(:,-1,0,0);neighbors(:,2)=q(:,1,0,0)
    neighbors(:,3)=q(:,0,-1,0);neighbors(:,4)=q(:,0,1,0)
    neighbors(:,5)=q(:,0,0,-1);neighbors(:,6)=q(:,0,0,1)
    n=0
    do k=0,1
    do j=0,1
    do i=0,1
      n=n+1;value=q(:,0,0,0);offset=.5_dp*real([i,j,k],dp)-.25_dp
      call interpolate_eb_linear_halo_3d(value,neighbors,offset,valid,separate_momentum=.true.)
      if(.not.valid) return
      base(:,n)=value
      do c=1,nv+1
        if(bounded) then
          high(c,n)=bounded_child_average(q(c,:,:,:),i,j,k,continuous)
        else
          high(c,n)=child_average(q(c,:,:,:),i,j,k)
        endif
      enddo
    enddo
    enddo
    enddo
    ! Center both families on the physical parent average. This is a
    ! conservation correction, not a cut-cell accuracy/centroid reconstruction.
    do c=1,nv+1
      base(c,:)=base(c,:)+(q(c,0,0,0)-sum(weights*base(c,:)))
      high(c,:)=high(c,:)+(q(c,0,0,0)-sum(weights*high(c,:)))
    enddo
    do n=1,8
      if(weights(n)>0) cycle
      base(:,n)=q(:,0,0,0);high(:,n)=q(:,0,0,0)
    enddo
    if(.not.valid_family(species,base,temperature,weights)) base=spread(q(:,0,0,0),2,8)
    if(bounded) then
      ! Candidate is limited continuously; smooth means no binary fallback.
      smooth=.true.
    else
      smooth=curvature_screen(q,momentum_roundoff_scale,inflection_hierarchy)
    endif
    if(.not.smooth) high=base
    call admit_family(species,base,high,temperature,children,theta,ok,child_volume)
  end subroutine

  subroutine admit_family(species,base,high,temperature,children,theta,ok,child_volume)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: base(:,:),high(:,:),temperature
    real(dp), intent(inout) :: children(:,:),theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: child_volume(8)
    real(dp) :: weights(8)
    logical :: good
    real(dp) :: trial(size(base,1),size(base,2)),weight,scale
    integer :: attempt,c
    ok=.false.
    if(any(shape(base)/=[reactive_nvar(size(species))+1,8])) return
    if(any(shape(base)/=shape(high)).or.any(shape(base)/=shape(children))) return
    call normalize_weights(child_volume,weights,good)
    if(.not.good) return
    if(.not.all(ieee_is_finite(high)).or..not.valid_family(species,base,temperature,weights)) return
    do c=1,size(base,1)
      scale=max(1._dp,maxval(abs(base(c,:))),maxval(abs(high(c,:))))
      if(abs(sum(weights*base(c,:))-sum(weights*high(c,:)))>128*epsilon(1._dp)*scale) return
    enddo
    weight=1
    do attempt=0,54
      if(attempt==54) weight=0
      trial=base+weight*(high-base)
      if(valid_family(species,trial,temperature,weights)) then
        children=trial;theta=weight;ok=.true.;return
      endif
      weight=.5_dp*weight
    enddo
  end subroutine

  subroutine normalize_weights(volume,weights,ok)
    real(dp), optional, intent(in) :: volume(8)
    real(dp), intent(out) :: weights(8)
    logical, intent(out) :: ok
    ok=.false.;weights=1._dp/8
    if(present(volume)) then
      if(.not.all(ieee_is_finite(volume))) return
      if(any(volume<0).or.maxval(volume)<=0) return
      weights=volume/maxval(volume)
      weights=weights/sum(weights)
    endif
    ok=.true.
  end subroutine

  logical function valid_family(species,u,temperature,weights) result(valid)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: u(:,:),temperature,weights(8)
    real(dp) :: primitive(reactive_nprim(size(species))),t,cs
    integer :: nv,n
    logical :: good
    valid=.false.;nv=reactive_nvar(size(species))
    if(.not.ieee_is_finite(temperature).or.temperature<=0) return
    if(size(u,1)/=nv+1.or..not.all(ieee_is_finite(u))) return
    do n=1,size(u,2)
      if(weights(n)<=0) cycle
      if(u(1,n)<=0.or.any(u(6:nv,n)<0)) return
      call reactive_conserved_to_primitive(species,u(:nv,n),temperature,primitive,t,cs,good, &
        internal_energy_density=u(nv+1,n))
      if(.not.good) return
      ! Independent total-energy admission catches excessive kinetic energy;
      ! held-energy EOS alone would not see it. No sign bound on formation energy.
      call reactive_conserved_to_primitive(species,u(:nv,n),temperature,primitive,t,cs,good)
      if(.not.good) return
    enddo
    valid=.true.
  end function
end module eb_guarded_quartic_mod
