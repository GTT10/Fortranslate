module les_mpi_eb_sdc_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,eb_regular_cell_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use reactive_1d_mod, only: reactive_conserved_to_primitive,reactive_nprim
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_ghosted_3d
  use les_filter_3d_mod, only: les_filter_weights
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d
  use les_mpi_amr_tree_mod, only: compute_mpi_amr_tree_les_rhs_3d
  implicit none
  private
  public :: compute_mpi_regular_eb_les_sdc_3d,compute_mpi_regular_eb_les_spatial_sdc_3d
contains
  ! The pinned LES model applies only to a completely regular geometry.
  ! Convert separate EB temporal donors and zero-based integrated face fluxes.
  subroutine compute_mpi_regular_eb_les_sdc_3d(tree,owners,geometry,species,state,temperature,spacing,faces, &
      interval,depth,cs,prt,sct,rates,flux,maximum_diffusivity,comm,ok, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      wall_temperature,wall_velocity,frame,domain_lower)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    real(dp), intent(in) :: spacing(3),interval,cs,prt,sct
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: maximum_diffusivity
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions(:)
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6),wall_velocity(:,:)
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    type(amr_tree_state_3d), allocatable :: packed_start,packed_end
    type(amr_tree_state_3d) :: candidate
    type(amr_tree_flux_3d), allocatable :: density_flux(:)
    type(eb_integrated_flux_3d), allocatable :: result(:)
    real(dp), allocatable :: stage_fractions(:)
    real(dp) :: largest,low_dt,high_dt
    integer :: rank,status,node,n(3),nv,ios,controls(6),low(6),high(6)
    logical :: valid,ready
    ok=.false.;nv=5+size(species)
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    controls=[depth,merge(1,0,present(ancestor_start)),merge(1,0,present(ancestor_end)), &
      merge(1,0,present(temperature_start)),merge(1,0,present(temperature_end)),merge(1,0,present(fractions))]
    call MPI_Allreduce(controls,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low/=high).or.depth< -1.or.depth>=tree%level_count()) return
    valid=ieee_is_finite(interval).and.interval>0
    do node=1,size(geometry)
      valid=valid.and.all(geometry(node)%cell_type==eb_regular_cell_3d)
      valid=valid.and.all(geometry(node)%volume_fraction==1)
      valid=valid.and.all(geometry(node)%x_face_fraction==1).and.all(geometry(node)%y_face_fraction==1).and. &
        all(geometry(node)%z_face_fraction==1)
    end do
    call agree(valid)
    if(.not.valid) return
    call MPI_Allreduce(interval,low_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(interval,high_dt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(low_dt/=high_dt) return
    if(depth>0) then
      if(any(controls(2:)/=1)) return
      valid=size(fractions)>=depth
      call agree(valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_start,nv,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_end,nv,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature_start,1,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature_end,1,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      allocate(packed_start,packed_end,stat=ios)
      valid=ios==0;call agree(valid)
      if(.not.valid) return
      allocate(packed_start%fields(size(owners)),packed_end%fields(size(owners)),stage_fractions(depth),stat=ios)
      valid=ios==0;call agree(valid)
      if(.not.valid) return
      stage_fractions=fractions(:depth)
      do node=1,size(owners)
        if(owners(node)/=rank.or.tree%nodes(node)%level>=depth) cycle
        n=tree%nodes(node)%shape
        allocate(packed_start%fields(node)%values(nv+1,n(1),n(2),n(3)), &
          packed_end%fields(node)%values(nv+1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) exit
        packed_start%fields(node)%values(:nv,:,:,:)=ancestor_start%fields(node)%values
        packed_start%fields(node)%values(nv+1,:,:,:)=temperature_start%fields(node)%values(1,:,:,:)
        packed_end%fields(node)%values(:nv,:,:,:)=ancestor_end%fields(node)%values
        packed_end%fields(node)%values(nv+1,:,:,:)=temperature_end%fields(node)%values(1,:,:,:)
      end do
      valid=ios==0;call agree(valid)
      if(.not.valid) return
    end if
    call compute_mpi_amr_tree_les_rhs_3d(tree,owners,species,state,temperature,spacing,faces,cs,prt,sct, &
      candidate,density_flux,largest,comm,ready,inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
      wall_temperature=wall_temperature,wall_velocity=wall_velocity,model=model,ci=ci,test_model=test_model, &
      test_ratio=test_ratio,grid_ratio=grid_ratio,harmonic_mean=harmonic_mean,inflow_mach=inflow_mach, &
      inflow_pressure=inflow_pressure,target_depth=depth,coarse_start=packed_start,coarse_end=packed_end, &
      fractions=stage_fractions,frame=frame,domain_lower=domain_lower)
    if(.not.ready) return
    allocate(result(size(owners)),stat=ios)
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      allocate(result(node)%x(nv,0:n(1),n(2),n(3)),result(node)%y(nv,n(1),0:n(2),n(3)), &
        result(node)%z(nv,n(1),n(2),0:n(3)),stat=ios)
      if(ios/=0) exit
      result(node)%x=interval*density_flux(node)%x
      result(node)%y=interval*density_flux(node)%y
      result(node)%z=interval*density_flux(node)%z
    end do
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    call move_alloc(candidate%fields,rates%fields)
    call move_alloc(result,flux)
    maximum_diffusivity=largest;ok=.true.
  contains
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: global
      integer :: ierr
      call MPI_Allreduce(condition,global,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if(ierr==MPI_SUCCESS) condition=global
    end subroutine
  end subroutine
  ! Evaluate grow3 directly from supplied stage halos. With frame present,
  ! published rotating energy is converted privately for thermodynamics.
  ! Dynamic filtering needs its full filter stencil in addition to grow3;
  ! insufficient halos reject transactionally, rather than interpolate RHS.
  subroutine compute_mpi_regular_eb_les_spatial_sdc_3d(tree,owners,geometry,species,old_halos, &
      interval,depth,cs,prt,sct,rates,spatial_rates,flux,maximum_diffusivity,comm,ok, &
      model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,frame)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: old_halos
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), intent(in) :: interval,cs,prt,sct
    type(amr_tree_state_3d), intent(inout) :: rates,spatial_rates
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: maximum_diffusivity
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    type(amr_tree_state_3d) :: candidate,spatial
    type(eb_integrated_flux_3d), allocatable :: result(:)
    real(dp), allocatable :: weights(:),gas(:,:,:,:),temp(:,:,:),rhs(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp) :: controls(9),reference_controls(9),largest,diffusion,global_diffusion,sound,isotropic
    real(dp) :: position(3),spacing(3),origin(3),physical_energy
    real(dp) :: primitive(reactive_nprim(size(species)))
    integer :: flags(6),reference_flags(6),rank,status,node,n(3),m(3),nv,ios,halo,needed,lo(3),hi(3),i,j,k
    integer :: ratio,grid
    character(len=32) :: selected,filter,reference_names(2),names(2)
    logical :: valid,ready,dynamic,harmonic
    ok=.false.;nv=5+size(species);selected='smagorinsky';filter='box';ratio=2;grid=1;harmonic=.false.;isotropic=0
    if(present(model)) selected=model
    if(present(test_model)) filter=test_model
    if(present(test_ratio)) ratio=test_ratio
    if(present(grid_ratio)) grid=grid_ratio
    if(present(harmonic_mean)) harmonic=harmonic_mean
    if(present(ci)) isotropic=ci
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    names=[selected,filter];reference_names=names
    flags=[depth,ratio,grid,merge(1,0,harmonic),merge(1,0,present(frame)),0]
    controls=[interval,cs,prt,sct,isotropic,0._dp,0._dp,0._dp,0._dp]
    reference_flags=flags
    call MPI_Bcast(reference_flags,6,MPI_INTEGER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=flags(5)==reference_flags(5)
    call agree(valid);if(.not.valid) return
    if(present(frame)) then
      valid=frame%is_valid()
      call agree(valid);if(.not.valid) return
      flags(6)=frame%axis;controls(6:9)=[frame%angular_speed,frame%axis_point]
    endif
    reference_flags=flags;reference_controls=controls
    call MPI_Bcast(reference_names,64,MPI_CHARACTER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Bcast(reference_flags,6,MPI_INTEGER,0,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Bcast(reference_controls,9,MPI_DOUBLE_PRECISION,0,comm,status)
    if(status/=MPI_SUCCESS) return
    valid=all(names==reference_names).and.all(flags==reference_flags).and. &
      all(ieee_is_finite(controls)).and.all(ieee_is_finite(reference_controls))
    call agree(valid);if(.not.valid) return
    valid=all(controls==reference_controls).and.interval>0.and.depth>=-1.and.depth<tree%level_count()
    valid=valid.and.any(selected==[character(len=32)::'smagorinsky','wale','vreman','dynamic'])
    call agree(valid);if(.not.valid) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready,padding=6)
    if(.not.ready) return
    dynamic=selected=='dynamic';needed=2
    if(dynamic) then
      call les_filter_weights(filter,ratio,weights,valid)
      if(valid) needed=(size(weights)-1)/2+5
    endif
    valid=valid.and.allocated(old_halos%fields)
    if(valid) valid=size(old_halos%fields)==size(owners)
    call agree(valid);if(.not.valid) return
    allocate(candidate%fields(size(owners)),spatial%fields(size(owners)),result(size(owners)),stat=ios)
    valid=ios==0;call agree(valid);if(.not.valid) return
    largest=0
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;m=n+6
      allocate(candidate%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp,stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
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
    call agree(valid);if(.not.valid) return
    call MPI_Allreduce(largest,global_diffusion,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    call move_alloc(candidate%fields,rates%fields);call move_alloc(spatial%fields,spatial_rates%fields)
    call move_alloc(result,flux);maximum_diffusivity=global_diffusion;ok=.true.
  contains
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: agreed
      integer :: ierr
      agreed=.false.
      call MPI_Allreduce(condition,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=ierr==MPI_SUCCESS.and.agreed
    end subroutine agree
  end subroutine compute_mpi_regular_eb_les_spatial_sdc_3d
end module
