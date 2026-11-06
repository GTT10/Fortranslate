module les_eb_sdc_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nprim,reactive_conserved_to_primitive
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,eb_regular_cell_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_ghosted_3d
  use les_filter_3d_mod, only: les_filter_weights
  implicit none
  private
  public :: compute_regular_eb_les_spatial_sdc_3d
contains
  ! Uses the same per-patch operator/crops as the MPI spatial callback, with
  ! serial stage/ancestor donors supplied by the EB SDC coordinator.
  subroutine compute_regular_eb_les_spatial_sdc_3d(tree,geometry,species,old_halos, &
      interval,depth,cs,prt,sct,rates,spatial_rates,flux,maximum_diffusivity,ok, &
      model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,frame)
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: old_halos
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), intent(in) :: interval,cs,prt,sct
    type(amr_tree_state_3d), intent(inout) :: rates,spatial_rates
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: maximum_diffusivity
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    type(amr_tree_state_3d) :: candidate,spatial
    type(eb_integrated_flux_3d), allocatable :: result(:)
    real(dp), allocatable :: weights(:),gas(:,:,:,:),temp(:,:,:),rhs(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: largest,diffusion,sound,isotropic
    real(dp) :: position(3),spacing(3),origin(3),physical_energy
    real(dp) :: primitive(reactive_nprim(size(species)))
    integer :: node,n(3),m(3),nv,ios,halo,needed,lo(3),hi(3),i,j,k
    integer :: ratio,grid
    character(len=32) :: selected,filter
    logical :: valid,dynamic,harmonic
    ok=.false.;nv=5+size(species);selected='smagorinsky';filter='box';ratio=2;grid=1;harmonic=.false.;isotropic=0
    if(present(model)) selected=model
    if(present(test_model)) filter=test_model
    if(present(test_ratio)) ratio=test_ratio
    if(present(grid_ratio)) grid=grid_ratio
    if(present(harmonic_mean)) harmonic=harmonic_mean
    if(present(ci)) isotropic=ci
    if(.not.tree%topology%is_valid().or..not.allocated(tree%nodes)) return
    if(size(tree%nodes)/=size(tree%topology%nodes).or.size(geometry)/=size(tree%nodes)) return
    if(.not.all(ieee_is_finite([interval,cs,prt,sct,isotropic]))) return
    if(interval<=0.or.depth< -1.or.depth>=tree%topology%level_count()) return
    if(.not.any(selected==[character(len=32)::'smagorinsky','wale','vreman','dynamic'])) return
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    endif
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(.not.geometry(node)%is_valid()) return
      if(any([geometry(node)%nx,geometry(node)%ny,geometry(node)%nz]/=n+12)) return
    enddo
    valid=.true.
    dynamic=selected=='dynamic';needed=2
    if(dynamic) then
      call les_filter_weights(filter,ratio,weights,valid)
      if(valid) needed=(size(weights)-1)/2+5
    endif
    valid=valid.and.allocated(old_halos%fields)
    if(valid) valid=size(old_halos%fields)==size(tree%nodes)
    if(.not.valid) return
    allocate(candidate%fields(size(tree%nodes)),spatial%fields(size(tree%nodes)),result(size(tree%nodes)),stat=ios)
    valid=ios==0;if(.not.valid) return
    largest=0
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape;m=n+6
      allocate(candidate%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      if(depth>=0.and.tree%topology%nodes(node)%level/=depth) cycle
      valid=all(geometry(node)%cell_type==eb_regular_cell_3d).and. &
        all(geometry(node)%volume_fraction==1).and.all(geometry(node)%x_face_fraction==1).and. &
        all(geometry(node)%y_face_fraction==1).and.all(geometry(node)%z_face_fraction==1)
      if(.not.valid.or..not.allocated(old_halos%fields(node)%values)) then
        valid=.false.;exit
      endif
      halo=(size(old_halos%fields(node)%values,2)-n(1))/2
      valid=halo>=needed+3.and.all(shape(old_halos%fields(node)%values)==[nv,n+2*halo])
      if(.not.valid) exit
      lo=halo-needed-2;hi=halo+n+needed+3
      lo=lo+[lbound(old_halos%fields(node)%values,2),lbound(old_halos%fields(node)%values,3), &
        lbound(old_halos%fields(node)%values,4)]-1
      hi=hi+[lbound(old_halos%fields(node)%values,2),lbound(old_halos%fields(node)%values,3), &
        lbound(old_halos%fields(node)%values,4)]-1
      gas=old_halos%fields(node)%values(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      allocate(temp(m(1)+2*needed,m(2)+2*needed,m(3)+2*needed),rhs(nv,m(1),m(2),m(3)), &
        fx(nv,0:m(1),m(2),m(3)),fy(nv,m(1),0:m(2),m(3)),fz(nv,m(1),m(2),0:m(3)), &
        spatial%fields(node)%values(nv,n(1)+12,n(2)+12,n(3)+12),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      spacing=[geometry(node)%dx,geometry(node)%dy,geometry(node)%dz]
      origin=[geometry(node)%x_lower,geometry(node)%y_lower,geometry(node)%z_lower]+6*spacing
      do k=1,size(temp,3)
      do j=1,size(temp,2)
      do i=1,size(temp,1)
        if(present(frame)) then
          position=origin+(real([i,j,k]-needed-3,dp)-.5_dp)*spacing
          physical_energy=gas(5,i,j,k)
          call convert_rotating_energy_3d(frame,position,gas(1,i,j,k),gas(5,i,j,k),.false.,physical_energy,valid)
          if(.not.valid) exit
          gas(5,i,j,k)=physical_energy
        endif
        call reactive_conserved_to_primitive(species,gas(:,i,j,k),300._dp,primitive,temp(i,j,k),sound,valid)
        if(.not.valid) exit
      enddo
      if(.not.valid) exit
      enddo
      if(.not.valid) exit
      enddo
      if(.not.valid) exit
      if(dynamic) then
        call compute_dynamic_les_rhs_ghosted_3d(species,gas,temp, &
          [geometry(node)%dx,geometry(node)%dy,geometry(node)%dz],filter,ratio,needed, &
          rhs,fx,fy,fz,valid,grid,harmonic,diffusion)
      else
        call compute_les_ghosted_rhs_3d(species,gas,temp,m, &
          [geometry(node)%dx,geometry(node)%dy,geometry(node)%dz],cs,prt,sct,rhs,diffusion,valid, &
          model=selected,ci=isotropic,fx=fx,fy=fy,fz=fz)
      endif
      if(.not.valid) exit
      if(present(frame)) then
        ! Pinned SGS diffusion transports no total mass. Hence the published
        ! energy rate and flux equal their physical-energy counterparts.
        ! A future mass flux must implement its potential transport explicitly.
        valid=all(rhs(1,:,:,:)==0).and.all(fx(1,:,:,:)==0).and. &
          all(fy(1,:,:,:)==0).and.all(fz(1,:,:,:)==0)
        if(.not.valid) exit
      endif
      allocate(result(node)%x(nv,0:n(1),n(2),n(3)),result(node)%y(nv,n(1),0:n(2),n(3)), &
        result(node)%z(nv,n(1),n(2),0:n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      candidate%fields(node)%values=rhs(:,4:n(1)+3,4:n(2)+3,4:n(3)+3)
      spatial%fields(node)%values=0
      spatial%fields(node)%values(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)=rhs
      result(node)%x=interval*fx(:,3:n(1)+3,4:n(2)+3,4:n(3)+3)
      result(node)%y=interval*fy(:,4:n(1)+3,3:n(2)+3,4:n(3)+3)
      result(node)%z=interval*fz(:,4:n(1)+3,4:n(2)+3,3:n(3)+3)
      valid=all(ieee_is_finite(rhs)).and.all(ieee_is_finite(result(node)%x)).and. &
        all(ieee_is_finite(result(node)%y)).and.all(ieee_is_finite(result(node)%z))
      if(.not.valid) exit
      valid=ieee_is_finite(diffusion).and.diffusion>=0
      if(.not.valid) exit
      largest=max(largest,diffusion)
      deallocate(gas,temp,rhs,fx,fy,fz)
    enddo
    if(.not.valid) return
    call move_alloc(candidate%fields,rates%fields);call move_alloc(spatial%fields,spatial_rates%fields)
    call move_alloc(result,flux);maximum_diffusivity=largest;ok=.true.
  end subroutine compute_regular_eb_les_spatial_sdc_3d
end module les_eb_sdc_mod
