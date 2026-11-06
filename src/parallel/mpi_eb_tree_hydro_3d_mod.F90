module mpi_eb_tree_hydro_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d
  use nasa7_thermo_mod, only: nasa7_species
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d, &
    advance_reactive_eb_plm_ctu_3d,compute_reactive_eb_flattening_3d,eb_plm_boundary_3d,compute_sdc_eb_hydro_3d, &
    limit_reactive_eb_flux_correction_3d,correct_reactive_eb_spatial_hydro_flux_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,validate_mpi_eb_tree_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d,fill_sdc_mpi_eb_tree_source_halo_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: fetch_mpi_amr_tree_cells_3d,validate_mpi_amr_tree_local_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d
  use mpi_eb_tree_redistribution_3d_mod, only: compute_state_redist_mpi_level_stage_3d
  use eb_reactive_redistribution_3d_mod, only: post_eb_spatial_density_rate_3d
  implicit none
  private
  public :: compute_godunov_state_redist_mpi_eb_tree_hydro_3d
  public :: advance_mpi_eb_tree_hydro_3d,compute_sdc_mpi_eb_tree_hydro_3d
  public :: compute_mol_mpi_eb_tree_hydro_3d
  public :: fill_mpi_eb_flux_allowance_halo_3d
  public :: apply_mpi_eb_flux_correction_3d
contains
  ! baseline/flux contain dt-integrated intensive faces. hydro_rate already
  ! includes the original flux divergence and wall sources. source is the
  ! complete additional forward rate, including any chemistry predictor.
  ! Shared input faces must agree. This reconciles same-level allowances;
  ! coarse/fine reflux needs its own composite donor budget before acceptance.
  subroutine apply_mpi_eb_flux_correction_3d(tree,owners,geometry,state,source,baseline,dt,hydro_rate,flux, &
      comm,ok,periodic,target_depth, &
      internal_energy_baseline,internal_energy_rate,internal_energy_flux)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,source
    type(eb_integrated_flux_3d), intent(in) :: baseline(:)
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    integer, optional, intent(in) :: target_depth
    type(eb_integrated_flux_3d), optional, intent(in) :: internal_energy_baseline(:)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    type(amr_tree_state_3d) :: held_candidate
    type(eb_integrated_flux_3d), allocatable :: held_limited(:)
    type(eb_integrated_flux_3d) :: held_low,held_original,held_delta
    real(dp) :: held_increment
    type(amr_tree_state_3d) :: candidate,allowance,halo
    type(eb_integrated_flux_3d), allocatable :: original(:),limited(:)
    type(eb_integrated_flux_3d) :: low,delta
    real(dp), allocatable :: endpoint(:)
    real(dp) :: controls(2),minimum(2),maximum(2),volume,areas(3)
    integer :: nv,nc,rank,status,ios,p,n(3),i,j,k,c,depth,held_options,held_min,held_max
    logical :: valid,ready
    ok=.false.;depth=-1
    if(present(target_depth)) depth=target_depth
    held_options=merge(1,0,present(internal_energy_baseline))+merge(2,0,present(internal_energy_rate))+ &
      merge(4,0,present(internal_energy_flux))
    call MPI_Allreduce(held_options,held_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_options,held_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_min/=held_max) return
    if(held_options/=0.and.held_options/=7) return
    controls=[dt,real(depth,dp)]
    valid=all(ieee_is_finite(controls)).and.dt>0.and.depth>=-1
    call agree(valid,ready);if(.not.ready) return
    call MPI_Allreduce(controls,minimum,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    valid=allocated(hydro_rate%fields).and.allocated(flux)
    call agree(valid,ready);if(.not.ready) return
    valid=size(hydro_rate%fields)==size(owners).and.size(flux)==size(owners).and.size(baseline)==size(owners)
    call agree(valid,ready);if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    nv=0
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(allocated(hydro_rate%fields(p)%values)) nv=size(hydro_rate%fields(p)%values,1)
      exit
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,nv,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    call validate_mpi_amr_tree_local_3d(tree,owners,hydro_rate,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,source,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call agree(nv>=6.and.depth<tree%level_count(),ready);if(.not.ready) return
    if(present(internal_energy_rate)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy_rate,1,comm,ready)
      if(.not.ready) return
      call agree(allocated(internal_energy_flux),ready);if(.not.ready) return
      call agree(size(internal_energy_flux)==size(owners).and.size(internal_energy_baseline)==size(owners),ready)
      if(.not.ready) return
      held_candidate=internal_energy_rate
      allocate(held_limited(size(owners)),stat=ios)
      call agree(ios==0,ready);if(.not.ready) return
      valid=.true.
      do p=1,size(owners)
        if(owners(p)/=rank) cycle
        if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
        valid=valid.and.valid_eb_integrated_flux_3d(internal_energy_baseline(p),geometry(p),1).and. &
          valid_eb_integrated_flux_3d(internal_energy_flux(p),geometry(p),1)
      enddo
      call agree(valid,ready);if(.not.ready) return
    endif
    nc=nv-4;candidate=hydro_rate
    allocate(allowance%fields(size(owners)),original(size(owners)),limited(size(owners)),endpoint(nv),stat=ios)
    call agree(ios==0,ready);if(.not.ready) return
    valid=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      allocate(allowance%fields(p)%values(nc,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      allowance%fields(p)%values=1
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      original(p)=flux(p);low=baseline(p)
      valid=allocated(original(p)%x).and.allocated(original(p)%y).and.allocated(original(p)%z).and. &
        allocated(low%x).and.allocated(low%y).and.allocated(low%z)
      if(.not.valid) exit
      original(p)%x=original(p)%x/dt;original(p)%y=original(p)%y/dt;original(p)%z=original(p)%z/dt
      low%x=low%x/dt;low%y=low%y/dt;low%z=low%z/dt
      limited(p)=original(p)
      call limit_reactive_eb_flux_correction_3d(geometry(p),state%fields(p)%values,low,limited(p),dt,valid, &
        source%fields(p)%values,allowance%fields(p)%values)
      if(.not.valid) exit
    enddo
    call agree(valid,ready);if(.not.ready) return
    call fill_mpi_eb_flux_allowance_halo_3d(tree,owners,allowance,nc,halo,comm,ready,periodic)
    if(.not.ready) return
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      low=baseline(p);low%x=low%x/dt;low%y=low%y/dt;low%z=low%z/dt
      limited(p)=original(p)
      if(present(internal_energy_rate)) then
        held_low=internal_energy_baseline(p)
        held_low%x=held_low%x/dt;held_low%y=held_low%y/dt;held_low%z=held_low%z/dt
        held_original=internal_energy_flux(p)
        held_original%x=held_original%x/dt;held_original%y=held_original%y/dt;held_original%z=held_original%z/dt
        held_limited(p)=held_original
        call limit_reactive_eb_flux_correction_3d(geometry(p),state%fields(p)%values,low,limited(p),dt,valid, &
          source=source%fields(p)%values,peer_alpha=halo%fields(p)%values, &
          internal_energy_baseline=held_low,internal_energy_flux=held_limited(p))
      else
      call limit_reactive_eb_flux_correction_3d(geometry(p),state%fields(p)%values,low,limited(p),dt,valid, &
        source=source%fields(p)%values,peer_alpha=halo%fields(p)%values)
      endif
      if(.not.valid) exit
      where(limited(p)%x==original(p)%x)
        limited(p)%x=flux(p)%x
      elsewhere
        limited(p)%x=dt*limited(p)%x
      endwhere
      where(limited(p)%y==original(p)%y)
        limited(p)%y=flux(p)%y
      elsewhere
        limited(p)%y=dt*limited(p)%y
      endwhere
      where(limited(p)%z==original(p)%z)
        limited(p)%z=flux(p)%z
      elsewhere
        limited(p)%z=dt*limited(p)%z
      endwhere
      if(present(internal_energy_rate)) then
        where(held_limited(p)%x==held_original%x)
          held_limited(p)%x=internal_energy_flux(p)%x
        elsewhere
          held_limited(p)%x=dt*held_limited(p)%x
        endwhere
        where(held_limited(p)%y==held_original%y)
          held_limited(p)%y=internal_energy_flux(p)%y
        elsewhere
          held_limited(p)%y=dt*held_limited(p)%y
        endwhere
        where(held_limited(p)%z==held_original%z)
          held_limited(p)%z=internal_energy_flux(p)%z
        elsewhere
          held_limited(p)%z=dt*held_limited(p)%z
        endwhere
        if(.not.valid_eb_integrated_flux_3d(held_limited(p),geometry(p),1)) then
          valid=.false.;exit
        endif
        held_delta=held_limited(p)
        held_delta%x=held_delta%x-internal_energy_flux(p)%x
        held_delta%y=held_delta%y-internal_energy_flux(p)%y
        held_delta%z=held_delta%z-internal_energy_flux(p)%z
      endif
      delta=limited(p);delta%x=delta%x-flux(p)%x;delta%y=delta%y-flux(p)%y;delta%z=delta%z-flux(p)%z
      n=tree%nodes(p)%shape
      associate(g=>geometry(p))
      areas=[g%dy*g%dz,g%dx*g%dz,g%dx*g%dy]
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(g%volume_fraction(i,j,k)<=0) cycle
        volume=dt*g%dx*g%dy*g%dz*g%volume_fraction(i,j,k)
        endpoint=( &
          areas(1)*(g%x_face_fraction(i-1,j,k)*delta%x(:,i-1,j,k)-g%x_face_fraction(i,j,k)*delta%x(:,i,j,k))+ &
          areas(2)*(g%y_face_fraction(i,j-1,k)*delta%y(:,i,j-1,k)-g%y_face_fraction(i,j,k)*delta%y(:,i,j,k))+ &
          areas(3)*(g%z_face_fraction(i,j,k-1)*delta%z(:,i,j,k-1)-g%z_face_fraction(i,j,k)*delta%z(:,i,j,k)))/volume
        if(any(endpoint/=0)) candidate%fields(p)%values(:,i,j,k)=hydro_rate%fields(p)%values(:,i,j,k)+endpoint
        ! Slip walls carry no density or species flux. Recompute changed
        ! components directly: high divergence plus its opposite correction
        ! can leave a negative rounding residual when a species flux is zero.
        ! Retain the original rate for unchanged components and wall terms.
        do c=1,nv
          if(c>=2.and.c<=5) cycle
          if(endpoint(c)==0) cycle
          candidate%fields(p)%values(c,i,j,k)=( &
            areas(1)*(g%x_face_fraction(i-1,j,k)*limited(p)%x(c,i-1,j,k)- &
              g%x_face_fraction(i,j,k)*limited(p)%x(c,i,j,k))+ &
            areas(2)*(g%y_face_fraction(i,j-1,k)*limited(p)%y(c,i,j-1,k)- &
              g%y_face_fraction(i,j,k)*limited(p)%y(c,i,j,k))+ &
            areas(3)*(g%z_face_fraction(i,j,k-1)*limited(p)%z(c,i,j,k-1)- &
              g%z_face_fraction(i,j,k)*limited(p)%z(c,i,j,k)))/volume
        enddo
        if(present(internal_energy_rate)) then
          held_increment=( &
            areas(1)*(g%x_face_fraction(i-1,j,k)*held_delta%x(1,i-1,j,k)- &
              g%x_face_fraction(i,j,k)*held_delta%x(1,i,j,k))+ &
            areas(2)*(g%y_face_fraction(i,j-1,k)*held_delta%y(1,i,j-1,k)- &
              g%y_face_fraction(i,j,k)*held_delta%y(1,i,j,k))+ &
            areas(3)*(g%z_face_fraction(i,j,k-1)*held_delta%z(1,i,j,k-1)- &
              g%z_face_fraction(i,j,k)*held_delta%z(1,i,j,k)))/volume
          if(held_increment/=0) held_candidate%fields(p)%values(1,i,j,k)= &
            internal_energy_rate%fields(p)%values(1,i,j,k)+held_increment
          if(.not.ieee_is_finite(held_candidate%fields(p)%values(1,i,j,k))) then
            valid=.false.;exit
          endif
        endif
        endpoint=state%fields(p)%values(:,i,j,k)+dt*(candidate%fields(p)%values(:,i,j,k)+source%fields(p)%values(:,i,j,k))
        valid=valid.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
      enddo
      enddo
      enddo
      end associate
      if(.not.valid) exit
    enddo
    call agree(valid,ready);if(.not.ready) return
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      if(depth>=0.and.tree%nodes(p)%level/=depth) cycle
      call move_alloc(limited(p)%x,flux(p)%x);call move_alloc(limited(p)%y,flux(p)%y)
      call move_alloc(limited(p)%z,flux(p)%z)
      if(present(internal_energy_rate)) then
        call move_alloc(held_limited(p)%x,internal_energy_flux(p)%x)
        call move_alloc(held_limited(p)%y,internal_energy_flux(p)%y)
        call move_alloc(held_limited(p)%z,internal_energy_flux(p)%z)
      endif
    enddo
    if(present(internal_energy_rate)) call move_alloc(held_candidate%fields,internal_energy_rate%fields)
    call move_alloc(candidate%fields,hydro_rate%fields);ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
  end subroutine apply_mpi_eb_flux_correction_3d

  ! Exchange same-level face donors, including periodic seams. Coarse/fine
  ! and physical boundary donors are left at one; their budgets must be
  ! supplied by the caller's separate coarse/fine or boundary policy.
  ! Interior allowances are indexed 1:N; returned owner-local halos 0:N+1.
  subroutine fill_mpi_eb_flux_allowance_halo_3d(tree,owners,allowance,ncomponents,halo,comm,ok,periodic)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),ncomponents
    type(amr_tree_state_3d), intent(in) :: allowance
    type(amr_tree_state_3d), intent(inout) :: halo
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    type(amr_tree_state_3d) :: candidate
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:),requests(:,:),records(:,:)
    real(dp), allocatable :: remote(:,:)
    integer(int64) :: face_total,global_cell(3)
    integer :: rank,status,ios,p,q,n(3),level,axis,side,i,j,k,cell(3),neighbor(3),peer,nrequests,r
    integer :: flags(3),low(3),high(3)
    logical :: valid,ready,wrap(3)
    ok=.false.;wrap=.false.
    if(present(periodic)) wrap=periodic
    call validate_mpi_amr_tree_local_3d(tree,owners,allowance,ncomponents,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    flags=merge(1,0,wrap)
    call MPI_Allreduce(flags,low,3,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,3,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    valid=ncomponents>=1;face_total=0
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      associate(a=>allowance%fields(p)%values)
      valid=valid.and.all(lbound(a)==1).and.all(ieee_is_finite(a))
      valid=valid.and.all(a>=0).and.all(a<=1)
      end associate
      face_total=face_total+2_int64*(int(n(1),int64)*n(2)+int(n(1),int64)*n(3)+int(n(2),int64)*n(3))
      if(face_total>huge(1)/5) then
        valid=.false.;exit
      endif
    enddo
    call agree(valid,ready);if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready);if(.not.ready) return
    allocate(candidate%fields(size(owners)),requests(4,int(face_total)),records(4,int(face_total)),stat=ios)
    call agree(ios==0,ready);if(.not.ready) return
    valid=.true.
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape
      allocate(candidate%fields(p)%values(ncomponents,0:n(1)+1,0:n(2)+1,0:n(3)+1),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      endif
      candidate%fields(p)%values=1
      candidate%fields(p)%values(:,1:n(1),1:n(2),1:n(3))=allowance%fields(p)%values
    enddo
    call agree(valid,ready);if(.not.ready) return
    nrequests=0
    do p=1,size(owners)
      if(owners(p)/=rank) cycle
      n=tree%nodes(p)%shape;level=tree%nodes(p)%level+1
      do axis=1,3
      do side=0,1
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        cell=[i,j,k]
        if(side==0.and.cell(axis)/=1) cycle
        if(side==1.and.cell(axis)/=n(axis)) cycle
        cell(axis)=cell(axis)+2*side-1
        global_cell=origin(:,p)+int(cell-1,int64)
        if(wrap(axis)) global_cell(axis)=modulo(global_cell(axis),domain(axis,level))
        peer=0
        do q=1,size(owners)
          if(tree%nodes(q)%level/=tree%nodes(p)%level) cycle
          if(all(global_cell>=origin(:,q)).and.all(global_cell<origin(:,q)+int(tree%nodes(q)%shape,int64))) then
            peer=q;exit
          endif
        enddo
        if(peer==0) cycle
        neighbor=int(global_cell-origin(:,peer))+1
        nrequests=nrequests+1
        requests(:,nrequests)=[peer,neighbor];records(:,nrequests)=[p,cell]
      enddo
      enddo
      enddo
      enddo
      enddo
    enddo
    call fetch_mpi_amr_tree_cells_3d(tree,owners,allowance,ncomponents,requests(:,:nrequests),remote,comm,ready)
    if(.not.ready) return
    valid=all(ieee_is_finite(remote)).and.all(remote>=0).and.all(remote<=1)
    call agree(valid,ready);if(.not.ready) return
    do r=1,nrequests
      p=records(1,r);cell=records(2:4,r)
      candidate%fields(p)%values(:,cell(1),cell(2),cell(3))=remote(:,r)
    enddo
    call move_alloc(candidate%fields,halo%fields);ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine agree
  end subroutine fill_mpi_eb_flux_allowance_halo_3d

  ! Common interval with PCM or an explicit six-layer PLM CTU geometry cache.
  subroutine advance_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,solver,dt,state,temperature, &
      boundary_change,wall_change,comm,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(eb_geometry_3d) :: center,working
    type(eb_plm_boundary_3d) :: physical_boundary
    type(amr_tree_state_3d) :: metrics
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: flat(:,:,:)
    logical :: ctu,reference_ctu
    real(dp) :: coefficient
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,ghost,gt,vf
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: source(:,:,:,:)
    real(dp) :: work_stamp
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),local_wall(size(wall_change)),low,high
    logical, allocatable :: exposed(:,:,:)
    integer :: nv,rank,status,node,child,n(3),lo(3),hi(3),c,ios
    character(len=32) :: reference
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(dt).and.ieee_is_finite(coefficient).and.coefficient>=0
    if (valid) valid=dt>0._dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.len_trim(solver)<=len(reference)
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high) return
    call MPI_Allreduce(coefficient,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    reference=solver
    call MPI_Bcast(reference,len(reference),MPI_CHARACTER,0,comm,status)
    call agree(status==MPI_SUCCESS.and.trim(solver)==trim(reference),ready)
    if (.not.ready) return
    ctu=present(extended_geometry);reference_ctu=ctu
    call MPI_Bcast(reference_ctu,1,MPI_LOGICAL,0,comm,status)
    call agree(status==MPI_SUCCESS.and.(ctu.eqv.reference_ctu),ready)
    if(.not.ready) return
    if(ctu) then
      call agree(solver=='pelec',ready)
      if(.not.ready) return
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
      if(.not.ready) return
    end if
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    if(ctu) then
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      if(valid) valid=all(domain<=int(huge(1)-4,int64))
      call agree(valid,ready)
      if(.not.ready) return
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      allocate(metrics%fields(size(owners)))
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
        if(valid) valid=matching_eb_metrics_3d(center,geometry(node))
        if(.not.valid) exit
        if(owners(node)/=rank) cycle
        allocate(metrics%fields(node)%values(1,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
        metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction
      end do
      call agree(valid,ready)
      if(.not.ready) return
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state, &
        ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=metrics)
    else if(coefficient>0) then
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf)
    else
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state)
    end if
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    candidate=snapshot;ct=st
    allocate(flux(size(owners)),stat=ios)
    valid=ios==0
    call agree(valid,ready)
    if (.not.ready) return
    local_wall=0._dp
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      work_stamp=tree_work_start()
      if(ctu) then
        call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,working,valid)
        if(.not.valid) exit
        allocate(flat(0:n(1)+7,0:n(2)+7,0:n(3)+7))
        call compute_reactive_eb_flattening_3d(species,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
          vf%fields(node)%values(1,:,:,:),n+6,flat,valid)
        if(.not.valid) exit
        physical_boundary%lower=4-int(origin(:,node))
        physical_boundary%upper=3+int(domain(:,tree%nodes(node)%level+1)-origin(:,node))
        physical_boundary%faces='outflow'
        if(present(boundary_faces)) physical_boundary%faces=boundary_faces
        if(present(inflow_mach)) then
          where(inflow_mach>=0) physical_boundary%faces='inflow'
        end if
        call advance_reactive_eb_plm_ctu_3d(species,working,ghost%fields(node)%values, &
          gt%fields(node)%values(1,:,:,:),vf%fields(node)%values(1,:,:,:),dt,flat,candidate%fields(node)%values, &
          ct%fields(node)%values(1,:,:,:),valid,difmag=coefficient,integrated_flux=flux(node),boundary=physical_boundary)
        if(.not.valid) exit
        do c=1,nv
          where(geometry(node)%volume_fraction<=0) candidate%fields(node)%values(c,:,:,:)=snapshot%fields(node)%values(c,:,:,:)
        end do
        where(geometry(node)%volume_fraction<=0) ct%fields(node)%values(1,:,:,:)=st%fields(node)%values(1,:,:,:)
        deallocate(flat)
      else if(coefficient>0) then
      call advance_reactive_eb_euler_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),solver,dt,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),flux(node), &
        difmag=coefficient,ghost_volume_fraction=vf%fields(node)%values(1,:,:,:))
      else
      call advance_reactive_eb_euler_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),solver,dt,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),flux(node),difmag=coefficient)
      end if
      call tree_work_stop(node,work_stamp)
      if (.not.valid) exit
      allocate(source(nv,n(1),n(2),n(3)),exposed(n(1),n(2),n(3)),stat=ios)
      valid=ios==0
      if (.not.valid) exit
      exposed=.true.
      do child=node+1,size(owners)
        if (tree%nodes(child)%parent/=node) cycle
        associate(p=>tree%nodes(child)%patch)
        lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
        hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
        end associate
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call reactive_eb_slip_wall_source_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),source,valid)
      if (.not.valid) exit
      associate(g=>geometry(node))
      do c=1,nv
        local_wall(c)=local_wall(c)+dt*g%dx*g%dy*g%dz*sum(source(c,:,:,:)*g%volume_fraction,mask=exposed)
      end do
      end associate
      deallocate(source,exposed)
    end do
    call agree(valid.and.all(ieee_is_finite(local_wall)),ready)
    if (.not.ready) return
    call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,flux,comm,ready, &
      boundary_faces=boundary_faces)
    if (.not.ready) return
    boundary=0._dp
    if (rank==owners(1)) then
      n=tree%nodes(1)%shape
      associate(g=>geometry(1),f=>flux(1))
      do c=1,nv
        boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
          g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
          g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
          g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
          g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
          g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
      end do
      end associate
    end if
    call MPI_Bcast(boundary,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(local_wall,wall,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.all(ieee_is_finite(wall)),ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
  ! SDC CTU stage: owner-local rates and integrated fluxes, without reflux or restriction.
  subroutine compute_godunov_state_redist_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,state,temperature,predictor, &
      extended_geometry,dt,hydro_rate,flux,comm,ok,difmag,boundary_faces,inflow_state, &
      inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      source_start,source_end,source_fraction,predictor_halo,baseline_flux,spatial_rate_halos,old_state_halos, &
      registers,max_order,target_volume_fraction,correct_forward_flux, &
      internal_energy_state,internal_energy_rate,internal_energy_flux,internal_energy_registers, &
      internal_energy_baseline_flux,internal_energy_spatial_rate_halos,internal_energy_old_state_halos, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,internal_energy_pressure_work,spatial_flux_weight)
    real(dp), optional, intent(in) :: spatial_flux_weight
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_state,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_spatial_rate_halos, &
      internal_energy_old_state_halos,internal_energy_pressure_work
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:),internal_energy_baseline_flux(:)
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: internal_energy_registers(:)
    type(amr_tree_state_3d) :: pending_held_rate,pending_held_source,pending_held_old,pending_held_work
    type(eb_integrated_flux_3d), allocatable :: pending_held_flux(:),pending_held_baseline(:)
    type(state_redist_register_3d), allocatable :: pending_held_registers(:)
    integer :: held_options
    logical :: held,held_valid
    type(amr_tree_state_3d), optional, intent(inout) :: spatial_rate_halos,old_state_halos
    logical, optional, intent(in) :: correct_forward_flux
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,predictor
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: baseline_flux(:)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    integer :: selected_depth,min_depth,max_depth
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end,predictor_halo
    real(dp), optional, intent(in) :: source_fraction(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(state_redist_register_3d), allocatable, intent(inout) :: registers(:)
    integer, optional, intent(in) :: max_order
    real(dp), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: pending_rate,pending_source,pending_old
    type(eb_integrated_flux_3d), allocatable :: pending_flux(:),pending_baseline(:)
    type(state_redist_register_3d), allocatable :: pending_registers(:)
    type(eb_geometry_3d), allocatable :: redist_geometry(:)
    integer :: p,n(3),nv,first_depth,last_depth,depth,padding
    logical :: good
    integer :: ios,status,presence,minimum_presence,maximum_presence
    logical :: ready
    ok=.false.;selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    held_options=count([present(internal_energy_state),present(internal_energy_rate), &
      present(internal_energy_flux),present(internal_energy_registers)])
    held=held_options==4
    held_valid=held_options==0.or.held_options==4
    held_valid=held_valid.and.(held.or..not.(present(internal_energy_baseline_flux).or. &
      present(internal_energy_spatial_rate_halos).or.present(internal_energy_old_state_halos).or. &
      present(internal_energy_pressure_work).or.present(ancestor_internal_energy_start).or. &
      present(ancestor_internal_energy_end)))
    held_valid=held_valid.and.(present(ancestor_internal_energy_start).eqv.present(ancestor_internal_energy_end))
    if(held) held_valid=held_valid.and.(present(ancestor_fraction).eqv.present(ancestor_internal_energy_start))
    presence=merge(1,0,present(baseline_flux))+2*merge(1,0,present(spatial_rate_halos))+ &
      4*merge(1,0,present(old_state_halos))+8*merge(1,0,held)+ &
      16*merge(1,0,present(internal_energy_baseline_flux))+32*merge(1,0,present(internal_energy_spatial_rate_halos))+ &
      64*merge(1,0,present(internal_energy_old_state_halos))+128*merge(1,0,present(internal_energy_pressure_work))+ &
      256*merge(1,0,present(ancestor_internal_energy_start))+512*merge(1,0,present(ancestor_internal_energy_end))
    call MPI_Allreduce(held_valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    call MPI_Allreduce(presence,minimum_presence,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(presence,maximum_presence,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.minimum_presence/=maximum_presence) return
    if(present(baseline_flux)) then
      if(held) then
      call compute_sdc_mpi_eb_tree_hydro_3d(species=species,tree=tree,owners=owners,geometry=geometry, &
        state=state,temperature=temperature,predictor=predictor,extended_geometry=extended_geometry,dt=dt, &
        hydro_rate=pending_rate,flux=pending_flux,comm=comm,ok=good,difmag=difmag, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
        temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,source_start=source_start, &
        source_end=source_end,source_fraction=source_fraction,predictor_halo=predictor_halo, &
        baseline_flux=pending_baseline,spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        internal_energy_state=internal_energy_state,internal_energy_rate=pending_held_rate, &
        internal_energy_flux=pending_held_flux,internal_energy_baseline_flux=pending_held_baseline, &
        internal_energy_spatial_rate_halos=pending_held_source,internal_energy_old_state_halos=pending_held_old, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end,internal_energy_pressure_work=pending_held_work, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_mpi_eb_tree_hydro_3d(species=species,tree=tree,owners=owners,geometry=geometry, &
        state=state,temperature=temperature,predictor=predictor,extended_geometry=extended_geometry,dt=dt, &
        hydro_rate=pending_rate,flux=pending_flux,comm=comm,ok=good,difmag=difmag, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
        temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,source_start=source_start, &
        source_end=source_end,source_fraction=source_fraction,predictor_halo=predictor_halo, &
        baseline_flux=pending_baseline,spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        spatial_flux_weight=spatial_flux_weight)
      endif
    else
      if(held) then
      call compute_sdc_mpi_eb_tree_hydro_3d(species=species,tree=tree,owners=owners,geometry=geometry, &
        state=state,temperature=temperature,predictor=predictor,extended_geometry=extended_geometry,dt=dt, &
        hydro_rate=pending_rate,flux=pending_flux,comm=comm,ok=good,difmag=difmag, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
        temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,source_start=source_start, &
        source_end=source_end,source_fraction=source_fraction,predictor_halo=predictor_halo, &
        spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        internal_energy_state=internal_energy_state,internal_energy_rate=pending_held_rate, &
        internal_energy_flux=pending_held_flux,internal_energy_baseline_flux=pending_held_baseline, &
        internal_energy_spatial_rate_halos=pending_held_source,internal_energy_old_state_halos=pending_held_old, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end,internal_energy_pressure_work=pending_held_work, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_mpi_eb_tree_hydro_3d(species=species,tree=tree,owners=owners,geometry=geometry, &
        state=state,temperature=temperature,predictor=predictor,extended_geometry=extended_geometry,dt=dt, &
        hydro_rate=pending_rate,flux=pending_flux,comm=comm,ok=good,difmag=difmag, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
        temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,source_start=source_start, &
        source_end=source_end,source_fraction=source_fraction,predictor_halo=predictor_halo, &
        spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        spatial_flux_weight=spatial_flux_weight)
      endif
    endif
    if(.not.good) return
    nv=reactive_nvar(size(species));first_depth=0
    if(present(target_depth)) first_depth=target_depth
    last_depth=first_depth
    if(.not.present(target_depth)) last_depth=tree%level_count()-1
    allocate(redist_geometry(size(owners)),stat=ios)
    call MPI_Allreduce(ios==0,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    do p=1,size(owners)
      n=tree%nodes(p)%shape
      padding=(extended_geometry(p)%nx-n(1))/2
      call extract_eb_geometry_3d(extended_geometry(p),[padding-5,padding-5,padding-5], &
        n+padding+6,redist_geometry(p),good)
      if(.not.good) exit
    enddo
    call MPI_Allreduce(good,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    ! Pinned Hydro.cpp uses factor1 for ML StateRedist. The caller owns
    ! provisional coarse accumulation and final-iteration fine/flux acceptance.
    do depth=first_depth,last_depth
      call compute_state_redist_mpi_level_stage_3d(tree,owners,redist_geometry,pending_old,pending_source, &
        nv,depth,dt,pending_rate,pending_registers,comm,good,boundary_faces,target_volume_fraction, &
        max_order,register_scale=1._dp)
      if(.not.good) return
      if(held) then
        call compute_state_redist_mpi_level_stage_3d(tree,owners,redist_geometry, &
          pending_held_old,pending_held_source,1,depth,dt,pending_held_rate,pending_held_registers,comm,good, &
          boundary_faces,target_volume_fraction,max_order,register_scale=1._dp)
        if(.not.good) return
      endif
    enddo
    good=.true.
    do p=1,size(owners)
      if(.not.allocated(pending_rate%fields(p)%values)) cycle
      if(selected_depth>=0.and.tree%nodes(p)%level/=selected_depth) cycle
      call post_eb_spatial_density_rate_3d(geometry(p),pending_rate%fields(p)%values,good)
      if(.not.good) exit
    enddo
    call MPI_Allreduce(good,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    if(held) then
      call move_alloc(pending_held_rate%fields,internal_energy_rate%fields)
      call move_alloc(pending_held_flux,internal_energy_flux)
      call move_alloc(pending_held_registers,internal_energy_registers)
      if(present(internal_energy_baseline_flux)) call move_alloc(pending_held_baseline,internal_energy_baseline_flux)
      if(present(internal_energy_spatial_rate_halos)) &
        call move_alloc(pending_held_source%fields,internal_energy_spatial_rate_halos%fields)
      if(present(internal_energy_old_state_halos)) &
        call move_alloc(pending_held_old%fields,internal_energy_old_state_halos%fields)
      if(present(internal_energy_pressure_work)) call move_alloc(pending_held_work%fields,internal_energy_pressure_work%fields)
    endif
    call move_alloc(pending_rate%fields,hydro_rate%fields)
    call move_alloc(pending_flux,flux)
    call move_alloc(pending_registers,registers)
    if(present(baseline_flux)) call move_alloc(pending_baseline,baseline_flux)
    ! These remain the raw spatial RHS and actual old state used by StateRedist.
    if(present(spatial_rate_halos)) call move_alloc(pending_source%fields,spatial_rate_halos%fields)
    if(present(old_state_halos)) call move_alloc(pending_old%fields,old_state_halos%fields)
    ok=.true.
  end subroutine compute_godunov_state_redist_mpi_eb_tree_hydro_3d

  subroutine compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,state,temperature,predictor, &
      extended_geometry,dt,hydro_rate,flux,comm,ok,difmag,boundary_faces,inflow_state, &
      inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      source_start,source_end,source_fraction,predictor_halo,baseline_flux,spatial_rate_halos,old_state_halos, &
      correct_forward_flux, &
      internal_energy_state,internal_energy_rate,internal_energy_flux,internal_energy_baseline_flux, &
      internal_energy_spatial_rate_halos,internal_energy_old_state_halos, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,internal_energy_pressure_work,spatial_flux_weight)
    real(dp), optional, intent(in) :: spatial_flux_weight
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_state,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_spatial_rate_halos, &
      internal_energy_old_state_halos,internal_energy_pressure_work
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:),internal_energy_baseline_flux(:)
    type(amr_tree_state_3d) :: held_rates,held_spatial,held_old,held_work
    type(eb_integrated_flux_3d), allocatable :: held_faces(:),held_baselines(:),held_local_baseline
    type(eb_integrated_flux_3d) :: held_full_flux
    real(dp), allocatable :: held_full_rate(:,:,:),held_full_work(:,:,:)
    integer :: held_options
    logical :: held,held_valid
    type(amr_tree_state_3d), optional, intent(inout) :: spatial_rate_halos,old_state_halos
    logical, optional, intent(in) :: correct_forward_flux
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,predictor
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: baseline_flux(:)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    integer :: selected_depth,min_depth,max_depth
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end,predictor_halo
    real(dp), optional, intent(in) :: source_fraction(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: rates,metrics,source_metrics,ghost,gt,vf,source,pending_halos,pending_old
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:)
    type(eb_integrated_flux_3d) :: full_flux
    type(amr_tree_state_3d) :: held_ghost
    type(eb_integrated_flux_3d), allocatable :: candidate_baseline(:),local_baseline
    type(eb_geometry_3d), allocatable :: halo_geometry(:)
    type(eb_geometry_3d) :: center,working,correction_geometry
    type(eb_plm_boundary_3d) :: physical_boundary
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: flat(:,:,:),full_rate(:,:,:,:)
    real(dp) :: coefficient,controls(3),low(3),high(3),work_stamp
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: nv,rank,status,ios,node,n(3),flag,loflag,hiflag,i,j,k,offset(3),width,grow,available_padding,offset_geometry
    integer :: predictor_padding,crop_offset(3)
    character(len=32) :: faces(6)
    logical :: valid,ready,corrected
    ok=.false.;nv=reactive_nvar(size(species))
    held_options=count([present(internal_energy_state),present(internal_energy_rate),present(internal_energy_flux)])
    held=held_options==3
    held_valid=held_options==0.or.held_options==3
    held_valid=held_valid.and.(held.or..not.(present(internal_energy_baseline_flux).or. &
      present(internal_energy_spatial_rate_halos).or.present(internal_energy_old_state_halos).or. &
      present(internal_energy_pressure_work).or.present(ancestor_internal_energy_start).or. &
      present(ancestor_internal_energy_end)))
    held_valid=held_valid.and.(present(ancestor_internal_energy_start).eqv.present(ancestor_internal_energy_end))
    if(held) held_valid=held_valid.and.(present(ancestor_fraction).eqv.present(ancestor_internal_energy_start))
    corrected=.false.
    if(present(correct_forward_flux)) corrected=correct_forward_flux
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    controls=[dt,coefficient,1._dp]
    if(present(spatial_flux_weight)) controls(3)=spatial_flux_weight
    call agree(held_valid.and.all(ieee_is_finite(controls)),ready)
    if(.not.ready) return
    call agree(dt>0.and.coefficient>=0.and.controls(3)>=0.and.controls(3)<=1,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    flag=merge(1,0,present(predictor_halo))+2*merge(1,0,present(baseline_flux))+ &
      4*merge(1,0,present(spatial_rate_halos))+8*merge(1,0,present(old_state_halos))+16*merge(1,0,corrected)+ &
      32*merge(1,0,held)+64*merge(1,0,present(internal_energy_baseline_flux))+ &
      128*merge(1,0,present(internal_energy_spatial_rate_halos))+ &
      256*merge(1,0,present(internal_energy_old_state_halos))+512*merge(1,0,present(internal_energy_pressure_work))+ &
      1024*merge(1,0,present(ancestor_internal_energy_start))+2048*merge(1,0,present(ancestor_internal_energy_end))+ &
      4096*merge(1,0,present(spatial_flux_weight))
    call MPI_Allreduce(flag,loflag,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,hiflag,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.loflag/=hiflag) return
    width=6
    if(present(spatial_rate_halos).or.present(old_state_halos).or.present(internal_energy_spatial_rate_halos).or. &
      present(internal_energy_old_state_halos)) width=9
    if(corrected) width=10
    grow=width-6
    selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    valid=selected_depth>=-1
    if(present(target_depth)) valid=selected_depth>=0
    if(present(ancestor_start).or.present(ancestor_end).or.present(temperature_start).or. &
      present(temperature_end).or.present(ancestor_fraction)) valid=valid.and.selected_depth>=0
    if(present(source_start).or.present(source_end).or.present(source_fraction)) valid=valid.and.selected_depth>=0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(selected_depth,min_depth,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selected_depth,max_depth,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.min_depth/=max_depth) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call agree(size(extended_geometry)==size(owners),ready)
    if(.not.ready) return
    allocate(halo_geometry(size(owners)),stat=ios)
    call agree(ios==0,ready);if(.not.ready) return
    valid=.true.
    do node=1,size(owners)
      valid=extended_geometry(node)%is_valid()
      if(.not.valid) exit
      n=tree%nodes(node)%shape;available_padding=(extended_geometry(node)%nx-n(1))/2
      valid=available_padding>=width.and.all([extended_geometry(node)%nx, &
        extended_geometry(node)%ny,extended_geometry(node)%nz]==n+2*available_padding)
      if(.not.valid) exit
      offset_geometry=available_padding-width
      call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]+offset_geometry, &
        n+2*width+offset_geometry,halo_geometry(node),valid)
      if(.not.valid) exit
    enddo
    call agree(valid,ready);if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,halo_geometry,nv,comm,ready,padding=width)
    if(.not.ready) return
    call agree(selected_depth<tree%level_count(),ready)
    if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(valid) valid=all(domain<=int(huge(1)-(width-2),int64))
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    if(held) then
      valid=stored_field_valid(internal_energy_state)
      if(present(ancestor_internal_energy_start)) valid=valid.and. &
        stored_field_valid(ancestor_internal_energy_start).and.stored_field_valid(ancestor_internal_energy_end)
      call agree(valid,ready);if(.not.ready) return
    endif
    allocate(metrics%fields(size(owners)),source_metrics%fields(size(owners)), &
      rates%fields(size(owners)),candidate_flux(size(owners)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    if(width>=9) then
      allocate(pending_halos%fields(size(owners)),pending_old%fields(size(owners)),stat=ios)
      call agree(ios==0,ready);if(.not.ready) return
    endif
    if(present(baseline_flux)) then
      allocate(candidate_baseline(size(owners)),stat=ios)
      call agree(ios==0,ready);if(.not.ready) return
    endif
    if(present(baseline_flux).or.corrected.or.present(internal_energy_baseline_flux)) then
      allocate(local_baseline,stat=ios)
      call agree(ios==0,ready);if(.not.ready) return
    endif
    if(held) then
      allocate(held_rates%fields(size(owners)),held_faces(size(owners)),stat=ios)
      call agree(ios==0,ready);if(.not.ready) return
      if(present(internal_energy_baseline_flux)) then
        allocate(held_baselines(size(owners)),stat=ios)
        call agree(ios==0,ready);if(.not.ready) return
      endif
      if(present(internal_energy_baseline_flux).or.corrected) then
        allocate(held_local_baseline ,stat=ios)
        call agree(ios==0,ready);if(.not.ready) return
      endif
      if(width>=9) then
        allocate(held_spatial%fields(size(owners)),held_old%fields(size(owners)),stat=ios)
        call agree(ios==0,ready);if(.not.ready) return
      endif
      if(present(internal_energy_pressure_work)) then
        allocate(held_work%fields(size(owners)),stat=ios)
        call agree(ios==0,ready);if(.not.ready) return
      endif
    endif
    do node=1,size(owners)
      n=tree%nodes(node)%shape
      call extract_eb_geometry_3d(halo_geometry(node),[1,1,1]+width,n+width,center,valid)
      if(valid) valid=matching_eb_metrics_3d(center,geometry(node))
      if(.not.valid) exit
      if(owners(node)/=rank) cycle
      allocate(metrics%fields(node)%values(1,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
        source_metrics%fields(node)%values(1,n(1)+2*width-6,n(2)+2*width-6,n(3)+2*width-6), &
        rates%fields(node)%values(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      metrics%fields(node)%values(1,:,:,:)=halo_geometry(node)%volume_fraction
      source_metrics%fields(node)%values(1,:,:,:)=halo_geometry(node)%volume_fraction( &
        4:n(1)+2*width-3,4:n(2)+2*width-3,4:n(3)+2*width-3)
      rates%fields(node)%values=0
      if(held) then
        allocate(held_rates%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        held_rates%fields(node)%values=0
        if(present(internal_energy_pressure_work)) then
          allocate(held_work%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
          if(ios/=0) then
            valid=.false.;exit
          endif
          held_work%fields(node)%values=0
        endif
      endif
    end do
    call agree(valid,ready)
    if(.not.ready) return
    if(held) then
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf,halo_width=width, &
      target_volume_fraction=metrics,restrict_snapshot=.false., &
      ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
      temperature_end=temperature_end,ancestor_fraction=ancestor_fraction, &
      internal_energy=internal_energy_state,ghost_internal_energy=held_ghost, &
      ancestor_internal_energy_start=ancestor_internal_energy_start, &
      ancestor_internal_energy_end=ancestor_internal_energy_end)
    else
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf,halo_width=width, &
      target_volume_fraction=metrics,restrict_snapshot=.false., &
      ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
      temperature_end=temperature_end,ancestor_fraction=ancestor_fraction)
    endif
    if(.not.ready) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(present(inflow_mach)) then
      where(inflow_mach>=0) faces='inflow'
    end if
    call fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,geometry,predictor,nv,width-3,source,comm,ready, &
      boundary_faces=faces,target_volume_fraction=source_metrics, &
      ancestor_start=source_start,ancestor_end=source_end,ancestor_fraction=source_fraction)
    if(.not.ready) return
    if(present(predictor_halo)) then
      valid=allocated(predictor_halo%fields)
      if(valid) valid=size(predictor_halo%fields)==size(owners)
      call agree(valid,ready);if(.not.ready) return
      do node=1,size(owners)
        valid=allocated(predictor_halo%fields(node)%values).eqv.(owners(node)==rank)
        if(.not.valid) exit
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        predictor_padding=(size(predictor_halo%fields(node)%values,2)-n(1))/2
        valid=predictor_padding>=width-3.and.predictor_padding<=7
        valid=valid.and.all(shape(predictor_halo%fields(node)%values)==[nv,n+2*predictor_padding])
        if(.not.valid) exit
        valid=all(ieee_is_finite(predictor_halo%fields(node)%values))
        if(.not.valid) exit
        offset=[lbound(predictor_halo%fields(node)%values,2),lbound(predictor_halo%fields(node)%values,3), &
          lbound(predictor_halo%fields(node)%values,4)]-1
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          valid=valid.and.all(predictor_halo%fields(node)%values(:,i+predictor_padding+offset(1), &
            j+predictor_padding+offset(2),k+predictor_padding+offset(3))== &
            source%fields(node)%values(:,i+width-3,j+width-3,k+width-3))
        end do
        end do
        end do
        if(.not.valid) exit
        crop_offset=offset+predictor_padding-(width-3)
        ! Actual same-level recipients, including periodic images, use the
        ! exchanged valid source. Cached virtual/ghost population corrections
        ! apply where no same-level valid cell exists. This keeps both copies
        ! of a peer face on the same predictor without changing the cache.
        do k=1,n(3)+2*(width-3)
        do j=1,n(2)+2*(width-3)
        do i=1,n(1)+2*(width-3)
          if(same_level_source([i,j,k])) cycle
          source%fields(node)%values(:,i,j,k)=predictor_halo%fields(node)%values(:, &
            i+crop_offset(1),j+crop_offset(2),k+crop_offset(3))
        end do
        end do
        end do
      end do
      call agree(valid,ready);if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(selected_depth>=0.and.tree%nodes(node)%level/=selected_depth) cycle
      n=tree%nodes(node)%shape
      call extract_eb_geometry_3d(halo_geometry(node),[4,4,4],n+2*width-3,working,valid)
      if(.not.valid) exit
      allocate(flat(0:n(1)+2*width-5,0:n(2)+2*width-5,0:n(3)+2*width-5), &
        full_rate(nv,n(1)+2*grow,n(2)+2*grow,n(3)+2*grow),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      call compute_reactive_eb_flattening_3d(species,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
        vf%fields(node)%values(1,:,:,:),n+2*width-6,flat,valid)
      if(.not.valid) exit
      physical_boundary%lower=width-2-int(origin(:,node))
      physical_boundary%upper=width-3+int(domain(:,tree%nodes(node)%level+1)-origin(:,node))
      physical_boundary%faces=faces
      work_stamp=tree_work_start()
      if(held) then
        allocate(held_full_rate(n(1)+2*grow,n(2)+2*grow,n(3)+2*grow), &
          held_full_work(n(1)+2*grow,n(2)+2*grow,n(3)+2*grow),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
      call compute_sdc_eb_hydro_3d(species,working,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
        vf%fields(node)%values(1,:,:,:),dt,flat,source%fields(node)%values,full_rate, &
        full_flux,boundary,wall,valid,difmag=coefficient,boundary=physical_boundary, &
        baseline_flux=local_baseline, &
        internal_energy_state=held_ghost%fields(node)%values(1,:,:,:),internal_energy_rate=held_full_rate, &
        internal_energy_integrated_flux=held_full_flux,internal_energy_baseline_flux=held_local_baseline, &
        internal_energy_pressure_work=held_full_work, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_eb_hydro_3d(species,working,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
        vf%fields(node)%values(1,:,:,:),dt,flat,source%fields(node)%values,full_rate, &
        full_flux,boundary,wall,valid,difmag=coefficient,boundary=physical_boundary, &
        baseline_flux=local_baseline, &
        spatial_flux_weight=spatial_flux_weight)
      endif
      call tree_work_stop(node,work_stamp)
      if(.not.valid) exit
      if(corrected) then
        call extract_eb_geometry_3d(working,[4,4,4],n+2*grow+3,correction_geometry,valid)
        if(.not.valid) exit
        if(held) then
        call correct_reactive_eb_spatial_hydro_flux_3d(correction_geometry, &
          ghost%fields(node)%values(:,1-grow:n(1)+grow,1-grow:n(2)+grow,1-grow:n(3)+grow), &
          source%fields(node)%values(:,4:n(1)+2*grow+3,4:n(2)+2*grow+3,4:n(3)+2*grow+3), &
          dt,full_rate,full_flux,local_baseline,valid, &
          internal_energy_rate=held_full_rate,internal_energy_flux=held_full_flux, &
          internal_energy_baseline=held_local_baseline)
        else
        call correct_reactive_eb_spatial_hydro_flux_3d(correction_geometry, &
          ghost%fields(node)%values(:,1-grow:n(1)+grow,1-grow:n(2)+grow,1-grow:n(3)+grow), &
          source%fields(node)%values(:,4:n(1)+2*grow+3,4:n(2)+2*grow+3,4:n(3)+2*grow+3), &
          dt,full_rate,full_flux,local_baseline,valid)
        endif
        if(.not.valid) exit
      endif
      rates%fields(node)%values=full_rate(:,grow+1:n(1)+grow,grow+1:n(2)+grow,grow+1:n(3)+grow)
      call crop_flux(full_flux,candidate_flux(node),n,grow,valid)
      if(.not.valid) exit
      if(present(baseline_flux)) then
        call crop_flux(local_baseline,candidate_baseline(node),n,grow,valid)
        if(.not.valid) exit
      endif
      if(held) then
        held_rates%fields(node)%values(1,:,:,:)= &
          held_full_rate(grow+1:n(1)+grow,grow+1:n(2)+grow,grow+1:n(3)+grow)
        call crop_flux(held_full_flux,held_faces(node),n,grow,valid)
        if(.not.valid) exit
        if(present(internal_energy_baseline_flux)) then
          call crop_flux(held_local_baseline,held_baselines(node),n,grow,valid)
          if(.not.valid) exit
        endif
        if(present(internal_energy_pressure_work)) held_work%fields(node)%values(1,:,:,:)= &
          held_full_work(grow+1:n(1)+grow,grow+1:n(2)+grow,grow+1:n(3)+grow)
        if(width>=9) then
          allocate(held_spatial%fields(node)%values(1,n(1)+12,n(2)+12,n(3)+12), &
            held_old%fields(node)%values(1,n(1)+12,n(2)+12,n(3)+12),stat=ios)
          valid=ios==0;if(.not.valid) exit
          held_spatial%fields(node)%values=0
          held_spatial%fields(node)%values(1,4:n(1)+9,4:n(2)+9,4:n(3)+9)= &
            held_full_rate(grow-2:n(1)+grow+3,grow-2:n(2)+grow+3,grow-2:n(3)+grow+3)
          held_old%fields(node)%values(1,:,:,:)=held_ghost%fields(node)%values(1,1-width+grow:n(1)+width-grow, &
            1-width+grow:n(2)+width-grow,1-width+grow:n(3)+width-grow)
        endif
        deallocate(held_full_rate,held_full_work)
      endif
      if(width>=9) then
        allocate(pending_halos%fields(node)%values(nv,n(1)+12,n(2)+12,n(3)+12), &
          pending_old%fields(node)%values(nv,n(1)+12,n(2)+12,n(3)+12),stat=ios)
        valid=ios==0;if(.not.valid) exit
        pending_halos%fields(node)%values=0
        pending_halos%fields(node)%values(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)= &
          full_rate(:,grow-2:n(1)+grow+3,grow-2:n(2)+grow+3,grow-2:n(3)+grow+3)
        pending_old%fields(node)%values=ghost%fields(node)%values(:,1-width+grow:n(1)+width-grow, &
          1-width+grow:n(2)+width-grow,1-width+grow:n(3)+width-grow)
      endif
      deallocate(flat,full_rate)
    end do
    call agree(valid,ready)
    if(.not.ready) return
    if(width>=9) then
      call refresh_peer_source_halos(ready)
      if(.not.ready) return
    endif
    if(present(spatial_rate_halos)) call move_alloc(pending_halos%fields,spatial_rate_halos%fields)
    if(present(old_state_halos)) call move_alloc(pending_old%fields,old_state_halos%fields)
    call move_alloc(rates%fields,hydro_rate%fields)
    if(present(baseline_flux)) call move_alloc(candidate_baseline,baseline_flux)
    if(held) then
      call move_alloc(held_rates%fields,internal_energy_rate%fields)
      call move_alloc(held_faces,internal_energy_flux)
      if(present(internal_energy_baseline_flux)) call move_alloc(held_baselines,internal_energy_baseline_flux)
      if(present(internal_energy_spatial_rate_halos)) &
        call move_alloc(held_spatial%fields,internal_energy_spatial_rate_halos%fields)
      if(present(internal_energy_old_state_halos)) call move_alloc(held_old%fields,internal_energy_old_state_halos%fields)
      if(present(internal_energy_pressure_work)) call move_alloc(held_work%fields,internal_energy_pressure_work%fields)
    endif
    call move_alloc(candidate_flux,flux);ok=.true.
  contains
    logical function stored_field_valid(field) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer :: p,d(3)
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(owners).or.lbound(field%fields,1)/=1) return
      do p=1,size(owners)
        if(owners(p)/=rank) then
          if(allocated(field%fields(p)%values)) return
          cycle
        endif
        if(.not.allocated(field%fields(p)%values)) return
        d=tree%nodes(p)%shape
        if(any(shape(field%fields(p)%values)/=[1,d]).or.any(lbound(field%fields(p)%values)/=1)) return
        if(any(.not.ieee_is_finite(field%fields(p)%values(1,:,:,:)).and.geometry(p)%volume_fraction>0)) return
      enddo
      good=.true.
    end function stored_field_valid
    logical function same_level_source(index) result(found)
      integer, intent(in) :: index(3)
      integer(int64) :: cell(3)
      integer :: a,p,level
      found=.false.;level=tree%nodes(node)%level+1
      cell=origin(:,node)+int(index-(width-3)-1,int64)
      do a=1,3
        if(faces(2*a-1)=='periodic') cell(a)=modulo(cell(a),domain(a,level))
      end do
      if(any(cell<0).or.any(cell>=domain(:,level))) return
      do p=1,size(owners)
        if(tree%nodes(p)%level/=level-1) cycle
        if(any(cell<origin(:,p)).or.any(cell>=origin(:,p)+int(tree%nodes(p)%shape,int64))) cycle
        found=.true.;return
      end do
    end function same_level_source

    subroutine crop_flux(full,cropped,n,g,success)
      type(eb_integrated_flux_3d), intent(in) :: full
      type(eb_integrated_flux_3d), intent(out) :: cropped
      integer, intent(in) :: n(3),g
      logical, intent(out) :: success
      integer :: allocation_status,nc
      nc=size(full%x,1)
      allocate(cropped%x(nc,0:n(1),n(2),n(3)),cropped%y(nc,n(1),0:n(2),n(3)), &
        cropped%z(nc,n(1),n(2),0:n(3)),stat=allocation_status)
      success=allocation_status==0;if(.not.success) return
      cropped%x=full%x(:,g:n(1)+g,g+1:n(2)+g,g+1:n(3)+g)
      cropped%y=full%y(:,g+1:n(1)+g,g:n(2)+g,g+1:n(3)+g)
      cropped%z=full%z(:,g+1:n(1)+g,g+1:n(2)+g,g:n(3)+g)
    end subroutine crop_flux
    subroutine refresh_peer_source_halos(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:),requests(:,:),destinations(:,:)
      real(dp), allocatable :: values(:,:)
      integer(int64) :: cell(3),total
      integer :: p,q,level,i,j,k,index(3),a,m(3),stage,r,allocation_status
      character(len=32) :: faces(6)
      logical :: inside,good
      faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
      call agree(good,success)
      if(.not.success) return
      total=0
      do stage=1,2
        r=0;good=.true.
        do p=1,size(owners)
          if(.not.allocated(pending_halos%fields(p)%values)) cycle
          level=tree%nodes(p)%level;m=tree%nodes(p)%shape
          do k=4,m(3)+9
          do j=4,m(2)+9
          do i=4,m(1)+9
            cell=origin(:,p)+int([i,j,k]-7,int64);inside=.true.
            do a=1,3
              if(faces(2*a-1)=='periodic') then
                cell(a)=modulo(cell(a),domain(a,level+1))
              else if(cell(a)<0.or.cell(a)>=domain(a,level+1)) then
                inside=.false.
              endif
            enddo
            if(.not.inside) cycle
            do q=1,size(owners)
              if(tree%nodes(q)%level/=level) cycle
              if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%nodes(q)%shape,int64))) cycle
              index=int(cell-origin(:,q))+1
              if(owners(q)==rank) then
                if(stage==2) then
                  pending_halos%fields(p)%values(:,i,j,k)=rates%fields(q)%values(:,index(1),index(2),index(3))
                  if(held) held_spatial%fields(p)%values(:,i,j,k)= &
                    held_rates%fields(q)%values(:,index(1),index(2),index(3))
                endif
              else if(stage==1) then
                total=total+1_int64
              else
                r=r+1;requests(:,r)=[q,index];destinations(:,r)=[p,i,j,k]
              endif
              exit
            enddo
          enddo
          enddo
          enddo
        enddo
        if(stage==1) then
          call agree(total<=int(huge(1)/max(4,nv),int64),success)
          if(.not.success) return
          allocate(requests(4,int(total)),destinations(4,int(total)),stat=allocation_status)
          call agree(allocation_status==0,success)
          if(.not.success) return
        endif
      enddo
      call fetch_mpi_amr_tree_cells_3d(tree,owners,rates,nv,requests,values,comm,success)
      if(.not.success) return
      do r=1,size(requests,2)
        p=destinations(1,r);index=destinations(2:4,r)
        pending_halos%fields(p)%values(:,index(1),index(2),index(3))=values(:,r)
      enddo
      if(held) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,held_rates,1,requests,values,comm,success)
        if(.not.success) return
        do r=1,size(requests,2)
          p=destinations(1,r);index=destinations(2:4,r)
          held_spatial%fields(p)%values(:,index(1),index(2),index(3))=values(:,r)
        enddo
      endif
    end subroutine refresh_peer_source_halos
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine compute_sdc_mpi_eb_tree_hydro_3d
  ! Spatial MOL evaluation only: no update, redistribution, restriction or reflux.
  ! Extended metrics contain three cells on each side; publication is collective.
  subroutine compute_mol_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,extended_geometry,state,temperature, &
      dt,order,hydro_rate,flux,comm,ok,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction,lax_friedrichs,frame, &
      spatial_rate_halos,old_state_halos,internal_energy_rate,internal_energy_flux, &
      internal_energy_spatial_rate_halos,internal_energy_state,internal_energy_old_state_halos, &
      ancestor_internal_energy_start,ancestor_internal_energy_end, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions)
    logical :: preserve_seed
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_spatial_rate_halos
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_old_state_halos
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_state,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: spatial_rate_halos,old_state_halos
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),order
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    logical, optional, intent(in) :: lax_friedrichs
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: rates,metrics,ghost,gt,vf,alpha,pending_halos,pending_old,auxiliary_rates, &
      pending_aux_halos,pending_aux_old,aux_ghost
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:),raw_flux(:),auxiliary_flux(:),auxiliary_raw(:)
    type(eb_geometry_3d) :: center
    real(dp), allocatable :: patch_halo(:,:,:,:),patch_aux_halo(:,:,:),held_input(:,:,:)
    real(dp) :: low,high,stamp,boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: nv,rank,status,node,n(3),ios,depth,halo_width,controls(10),minimum(10),maximum(10)
    logical :: valid,ready,use_lax,use_total_energy
    preserve_seed=.not.present(frame)
    if(present(separate_momentum)) preserve_seed=preserve_seed.and..not.separate_momentum
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(target_depth)) depth=target_depth
    valid=ieee_is_finite(dt).and.dt>0.and.(order==1.or.order==2).and.depth>=-1
    if(present(target_depth)) valid=valid.and.depth>=0
    if(present(ancestor_start).or.present(ancestor_end).or.present(temperature_start).or. &
      present(temperature_end).or.present(ancestor_fraction)) valid=valid.and.depth>=0
    valid=valid.and.(present(internal_energy_rate).eqv.present(internal_energy_flux))
    valid=valid.and.(.not.present(internal_energy_spatial_rate_halos).or.present(internal_energy_rate))
    valid=valid.and.(present(internal_energy_state).eqv.present(internal_energy_old_state_halos))
    valid=valid.and.(.not.present(internal_energy_state).or.present(internal_energy_rate))
    valid=valid.and.(present(ancestor_internal_energy_start).eqv.present(ancestor_internal_energy_end))
    valid=valid.and.(.not.present(ancestor_internal_energy_start).or.present(internal_energy_state))
    if(present(internal_energy_state)) &
      valid=valid.and.(present(ancestor_fraction).eqv.present(ancestor_internal_energy_start))
    if(present(separate_momentum)) &
      valid=valid.and.(.not.separate_momentum.or.present(internal_energy_old_state_halos))
    call agree(valid,ready)
    if(.not.ready) return
    controls=[order,depth,0,merge(1,0,present(spatial_rate_halos)),merge(1,0,present(old_state_halos)), &
      merge(1,0,present(internal_energy_rate)),merge(1,0,present(internal_energy_spatial_rate_halos)), &
      merge(1,0,present(internal_energy_old_state_halos)),0,merge(1,0,present(frame))]
    if(present(separate_momentum)) controls(9)=merge(1,0,separate_momentum)
    if(present(consistent_velocity)) controls(9)=controls(9)+merge(2,0,consistent_velocity)
    if(present(paired_pressure)) controls(9)=controls(9)+merge(4,0,paired_pressure)
    if(present(connected_directions)) controls(9)=controls(9)+merge(8,0,connected_directions)
    if(present(lax_friedrichs)) controls(3)=merge(1,0,lax_friedrichs)
    call MPI_Allreduce(controls,minimum,10,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,maximum,10,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(minimum/=maximum)) return
    if(iand(controls(9),3)==2) return
    if(btest(controls(9),3).and..not.btest(controls(9),0)) return
    if(btest(controls(9),2).and.iand(controls(9),7)/=7) return
    ! Match serial primitive and wall energy selection on every owner.
    use_total_energy=controls(10)==0.and..not.btest(controls(9),0)
    use_lax=controls(3)==1
    halo_width=3
    if(any(controls(4:5)==1).or.controls(7)==1.or.controls(8)==1) halo_width=6
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=halo_width)
    if(.not.ready) return
    call agree(depth<tree%level_count(),ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(rates%fields(size(owners)),metrics%fields(size(owners)),candidate_flux(size(owners)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    if(present(internal_energy_rate)) then
      allocate(auxiliary_rates%fields(size(owners)),auxiliary_flux(size(owners)),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
      if(use_lax) then
        allocate(auxiliary_raw(size(owners)),stat=ios)
        call agree(ios==0,ready)
        if(.not.ready) return
      endif
    endif
    if(halo_width==6) then
      allocate(pending_halos%fields(size(owners)),pending_old%fields(size(owners)),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
    endif
    if(present(internal_energy_spatial_rate_halos)) then
      allocate(pending_aux_halos%fields(size(owners)),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
    endif
    if(present(internal_energy_old_state_halos)) then
      allocate(pending_aux_old%fields(size(owners)),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
    endif
    if(use_lax) then
      allocate(alpha%fields(size(owners)),raw_flux(size(owners)),stat=ios)
      call agree(ios==0,ready)
      if(.not.ready) return
    end if
    do node=1,size(owners)
      n=tree%nodes(node)%shape
      call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]+halo_width,n+halo_width,center,valid)
      if(valid) valid=matching_eb_metrics_3d(center,geometry(node))
      if(.not.valid) exit
      if(owners(node)/=rank) cycle
      allocate(metrics%fields(node)%values(1,1-halo_width:n(1)+halo_width, &
        1-halo_width:n(2)+halo_width,1-halo_width:n(3)+halo_width), &
        rates%fields(node)%values(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction
      rates%fields(node)%values=0
      if(present(internal_energy_rate)) then
        allocate(auxiliary_rates%fields(node)%values(1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        endif
        auxiliary_rates%fields(node)%values=0
      endif
      if(use_lax) then
        allocate(alpha%fields(node)%values(nv-5,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
        alpha%fields(node)%values=1
      end if
    end do
    call agree(valid,ready)
    if(.not.ready) return
    if(present(internal_energy_old_state_halos)) then
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=halo_width,target_volume_fraction=metrics,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
        temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,frame=frame, &
        internal_energy=internal_energy_state,ghost_internal_energy=aux_ghost, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_geometry=extended_geometry, &
          separate_momentum=separate_momentum, &
        consistent_velocity=consistent_velocity,paired_pressure=paired_pressure,pressure_level=target_depth, &
      connected_directions=connected_directions,preserve_temperature_seed=preserve_seed)
    else
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=halo_width,target_volume_fraction=metrics,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
        temperature_end=temperature_end,ancestor_fraction=ancestor_fraction,frame=frame)
    endif
    if(.not.ready) return
    if(present(internal_energy_old_state_halos)) then
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        call move_alloc(aux_ghost%fields(node)%values,pending_aux_old%fields(node)%values)
      enddo
    endif
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      stamp=tree_work_start()
      if(present(internal_energy_old_state_halos)) held_input=pending_aux_old%fields(node)%values(1,:,:,:)
      if(halo_width==6) then
        if(use_lax) then
          if(present(internal_energy_rate)) then
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
              gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
              lax_friedrichs,frame,raw_flux(node),alpha%fields(node)%values,.true.,patch_halo, &
              internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
              internal_energy_integrated_flux=auxiliary_flux(node),internal_energy_raw_face_flux=auxiliary_raw(node), &
              internal_energy_spatial_rate_halo=patch_aux_halo,internal_energy_density=held_input, &
              primitive_from_total_energy=use_total_energy)
          else
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
              gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
              lax_friedrichs,frame,raw_flux(node),alpha%fields(node)%values,.true.,patch_halo, &
              primitive_from_total_energy=use_total_energy)
          endif
        else
          if(present(internal_energy_rate)) then
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
              gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
              lax_friedrichs,frame,spatial_rate_halo=patch_halo, &
              internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
              internal_energy_integrated_flux=auxiliary_flux(node), &
              internal_energy_spatial_rate_halo=patch_aux_halo,internal_energy_density=held_input, &
              primitive_from_total_energy=use_total_energy)
          else
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
              gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
              lax_friedrichs,frame,spatial_rate_halo=patch_halo, &
              primitive_from_total_energy=use_total_energy)
          endif
        endif
        if(valid) then
          call move_alloc(patch_halo,pending_halos%fields(node)%values)
          if(present(internal_energy_spatial_rate_halos)) then
            allocate(pending_aux_halos%fields(node)%values(1,size(patch_aux_halo,1), &
              size(patch_aux_halo,2),size(patch_aux_halo,3)),stat=ios)
            valid=ios==0
            if(valid) pending_aux_halos%fields(node)%values(1,:,:,:)=patch_aux_halo
          endif
          if(allocated(patch_aux_halo)) deallocate(patch_aux_halo)
          if(valid) then
            n=tree%nodes(node)%shape
            allocate(pending_old%fields(node)%values(nv,n(1)+12,n(2)+12,n(3)+12),stat=ios)
            valid=ios==0
            if(valid) pending_old%fields(node)%values=ghost%fields(node)%values
          endif
        endif
      else
      if(use_lax) then
        if(present(internal_energy_rate)) then
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
            gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
            lax_friedrichs,frame,raw_flux(node),alpha%fields(node)%values,.true., &
            internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
            internal_energy_integrated_flux=auxiliary_flux(node),internal_energy_raw_face_flux=auxiliary_raw(node), &
            internal_energy_density=held_input, &
              primitive_from_total_energy=use_total_energy)
        else
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
            gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
            lax_friedrichs,frame,raw_flux(node),alpha%fields(node)%values,.true., &
              primitive_from_total_energy=use_total_energy)
        endif
      else
        if(present(internal_energy_rate)) then
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
            gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
            lax_friedrichs,frame,internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
            internal_energy_integrated_flux=auxiliary_flux(node),internal_energy_density=held_input, &
              primitive_from_total_energy=use_total_energy)
        else
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost%fields(node)%values, &
            gt%fields(node)%values(1,:,:,:),dt,order,rates%fields(node)%values,candidate_flux(node),boundary,wall,valid, &
            lax_friedrichs,frame, &
              primitive_from_total_energy=use_total_energy)
        endif
      end if
      endif
      if(allocated(held_input)) deallocate(held_input)
      call tree_work_stop(node,stamp)
      if(.not.valid) exit
    end do
    call agree(valid,ready)
    if(.not.ready) return
    if(use_lax) then
      call limit_boundary_faces(ready)
      if(.not.ready) return
    end if
    if(present(internal_energy_spatial_rate_halos)) then
      do node=1,size(owners)
        if(.not.allocated(pending_aux_halos%fields(node)%values)) cycle
        n=tree%nodes(node)%shape
        pending_aux_halos%fields(node)%values(1,7:n(1)+6,7:n(2)+6,7:n(3)+6)= &
          auxiliary_rates%fields(node)%values(1,:,:,:)
      enddo
    endif
    if(halo_width==6) then
      call refresh_peer_source_halos(ready)
      if(.not.ready) return
    endif
    call move_alloc(rates%fields,hydro_rate%fields);call move_alloc(candidate_flux,flux)
    if(present(internal_energy_rate)) then
      call move_alloc(auxiliary_rates%fields,internal_energy_rate%fields)
      call move_alloc(auxiliary_flux,internal_energy_flux)
    endif
    if(present(spatial_rate_halos)) call move_alloc(pending_halos%fields,spatial_rate_halos%fields)
    if(present(old_state_halos)) call move_alloc(pending_old%fields,old_state_halos%fields)
    if(present(internal_energy_spatial_rate_halos)) &
      call move_alloc(pending_aux_halos%fields,internal_energy_spatial_rate_halos%fields)
    if(present(internal_energy_old_state_halos)) &
      call move_alloc(pending_aux_old%fields,internal_energy_old_state_halos%fields)
    ok=.true.
  contains
    ! Refresh only accepted same-level valid donors after face limiting.
    ! Locally evaluated coarse/fine ghost RHS remains at its own resolution.
    ! Gas, full halos, and registers are never gathered onto a rank.
    subroutine refresh_peer_source_halos(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:),requests(:,:),destinations(:,:)
      real(dp), allocatable :: values(:,:),aux_values(:,:)
      integer(int64) :: cell(3),total
      integer :: p,q,level,i,j,k,index(3),a,m(3),stage,r,allocation_status
      character(len=32) :: faces(6)
      logical :: inside,good
      faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
      call agree(good,success)
      if(.not.success) return
      total=0
      do stage=1,2
        r=0;good=.true.
        do p=1,size(owners)
          if(.not.allocated(pending_halos%fields(p)%values)) cycle
          level=tree%nodes(p)%level;m=tree%nodes(p)%shape
          do k=4,m(3)+9
          do j=4,m(2)+9
          do i=4,m(1)+9
            cell=origin(:,p)+int([i,j,k]-7,int64);inside=.true.
            do a=1,3
              if(faces(2*a-1)=='periodic') then
                cell(a)=modulo(cell(a),domain(a,level+1))
              else if(cell(a)<0.or.cell(a)>=domain(a,level+1)) then
                inside=.false.
              endif
            enddo
            if(.not.inside) cycle
            do q=1,size(owners)
              if(tree%nodes(q)%level/=level) cycle
              if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%nodes(q)%shape,int64))) cycle
              index=int(cell-origin(:,q))+1
              if(owners(q)==rank) then
                if(stage==2) then
                  pending_halos%fields(p)%values(:,i,j,k)=rates%fields(q)%values(:,index(1),index(2),index(3))
                  if(present(internal_energy_spatial_rate_halos)) &
                    pending_aux_halos%fields(p)%values(1,i,j,k)= &
                    auxiliary_rates%fields(q)%values(1,index(1),index(2),index(3))
                endif
              else if(stage==1) then
                total=total+1_int64
              else
                r=r+1;requests(:,r)=[q,index];destinations(:,r)=[p,i,j,k]
              endif
              exit
            enddo
          enddo
          enddo
          enddo
        enddo
        if(stage==1) then
          call agree(total<=int(huge(1)/max(4,nv),int64),success)
          if(.not.success) return
          allocate(requests(4,int(total)),destinations(4,int(total)),stat=allocation_status)
          call agree(allocation_status==0,success)
          if(.not.success) return
        endif
      enddo
      call fetch_mpi_amr_tree_cells_3d(tree,owners,rates,nv,requests,values,comm,success)
      if(.not.success) return
      if(present(internal_energy_spatial_rate_halos)) then
        call fetch_mpi_amr_tree_cells_3d(tree,owners,auxiliary_rates,1,requests,aux_values,comm,success)
        if(.not.success) return
      endif
      do r=1,size(requests,2)
        p=destinations(1,r);index=destinations(2:4,r)
        pending_halos%fields(p)%values(:,index(1),index(2),index(3))=values(:,r)
        if(present(internal_energy_spatial_rate_halos)) &
          pending_aux_halos%fields(p)%values(1,index(1),index(2),index(3))=aux_values(1,r)
      enddo
    end subroutine refresh_peer_source_halos

    subroutine limit_boundary_faces(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:),records(:,:),requests(:,:),query_for_record(:)
      real(dp), allocatable :: remote(:,:)
      real(dp) :: neighbor_alpha(nv-5),delta(nv),new_flux(nv),delta_aux,theta,face_area,fluid_volume,sign_face
      integer(int64) :: global_cell(3),coarse_cell(3),face_total
      integer :: face(3),inside(3),neighbor(3),level,p,axis,side,i,j,k,r,q,peer,max_faces,nrecords,nrequests
      character(len=16) :: faces(6)
      logical :: found
      success=.false.;faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,found)
      call agree(found,success)
      if(.not.success) return
      face_total=0;found=.true.
      do p=1,size(owners)
        if(owners(p)/=rank.or.(depth>=0.and.tree%nodes(p)%level/=depth)) cycle
        n=tree%nodes(p)%shape
        face_total=face_total+2_int64*(int(n(2),int64)*n(3)+int(n(1),int64)*n(3)+int(n(1),int64)*n(2))
        if(face_total>huge(1)) then
          found=.false.;exit
        end if
      end do
      call agree(found,success)
      if(.not.success) return
      max_faces=int(face_total)
      allocate(records(5,max_faces),requests(4,max_faces),query_for_record(max_faces),stat=ios)
      call agree(ios==0,success)
      if(.not.success) return
      nrecords=0;nrequests=0
      do p=1,size(owners)
        if(owners(p)/=rank.or.(depth>=0.and.tree%nodes(p)%level/=depth)) cycle
        n=tree%nodes(p)%shape;level=tree%nodes(p)%level+1
        do axis=1,3
        do side=0,1
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            face=[i,j,k]
            if(side==0.and.face(axis)/=1) cycle
            if(side==1.and.face(axis)/=n(axis)) cycle
            if(side==0) face(axis)=0
            select case(axis)
            case(1)
              if(all(candidate_flux(p)%x(6:,face(1),face(2),face(3))== &
                raw_flux(p)%x(6:,face(1),face(2),face(3)))) cycle
            case(2)
              if(all(candidate_flux(p)%y(6:,face(1),face(2),face(3))== &
                raw_flux(p)%y(6:,face(1),face(2),face(3)))) cycle
            case(3)
              if(all(candidate_flux(p)%z(6:,face(1),face(2),face(3))== &
                raw_flux(p)%z(6:,face(1),face(2),face(3)))) cycle
            end select
            inside=face
            if(side==0) inside(axis)=1
            nrecords=nrecords+1;records(:,nrecords)=[p,axis,face];query_for_record(nrecords)=0
            global_cell=origin(:,p)+int(inside-1,int64)
            global_cell(axis)=global_cell(axis)+int(2*side-1,int64)
            if(faces(2*axis-1)=='periodic'.and.faces(2*axis)=='periodic') &
              global_cell(axis)=modulo(global_cell(axis),domain(axis,level))
            peer=0
            do q=1,size(owners)
              if(tree%nodes(q)%level/=tree%nodes(p)%level) cycle
              if(all(global_cell>=origin(:,q)).and. &
                all(global_cell<origin(:,q)+int(tree%nodes(q)%shape,int64))) then
                peer=q;exit
              end if
            end do
            if(peer==0.and.level>1) then
              if(all(global_cell>=0_int64).and.all(global_cell<domain(:,level))) then
                coarse_cell=global_cell/int(ratios(level),int64)
                do q=1,size(owners)
                  if(tree%nodes(q)%level/=tree%nodes(p)%level-1) cycle
                  if(all(coarse_cell>=origin(:,q)).and. &
                    all(coarse_cell<origin(:,q)+int(tree%nodes(q)%shape,int64))) then
                    peer=q;exit
                  end if
                end do
              end if
            end if
            if(peer==0) cycle
            if(tree%nodes(peer)%level==tree%nodes(p)%level) then
              neighbor=int(global_cell-origin(:,peer))+1
            else
              neighbor=int(coarse_cell-origin(:,peer))+1
            end if
            nrequests=nrequests+1;requests(:,nrequests)=[peer,neighbor]
            query_for_record(nrecords)=nrequests
          end do
          end do
          end do
        end do
        end do
      end do
      call fetch_mpi_amr_tree_cells_3d(tree,owners,alpha,nv-5,requests(:,:nrequests),remote,comm,success)
      if(.not.success) return
      found=.true.
      do r=1,nrecords
        p=records(1,r);axis=records(2,r);face=records(3:5,r);n=tree%nodes(p)%shape
        inside=face;side=1
        if(face(axis)==0) then
          inside(axis)=1;side=0
        end if
        neighbor_alpha=1
        if(query_for_record(r)>0) neighbor_alpha=remote(:,query_for_record(r))
        select case(axis)
        case(1)
          delta=candidate_flux(p)%x(:,face(1),face(2),face(3))-raw_flux(p)%x(:,face(1),face(2),face(3))
          if(present(internal_energy_rate)) delta_aux=auxiliary_flux(p)%x(1,face(1),face(2),face(3))- &
            auxiliary_raw(p)%x(1,face(1),face(2),face(3))
          face_area=geometry(p)%dy*geometry(p)%dz*geometry(p)%x_face_fraction(face(1),face(2),face(3))
        case(2)
          delta=candidate_flux(p)%y(:,face(1),face(2),face(3))-raw_flux(p)%y(:,face(1),face(2),face(3))
          if(present(internal_energy_rate)) delta_aux=auxiliary_flux(p)%y(1,face(1),face(2),face(3))- &
            auxiliary_raw(p)%y(1,face(1),face(2),face(3))
          face_area=geometry(p)%dx*geometry(p)%dz*geometry(p)%y_face_fraction(face(1),face(2),face(3))
        case(3)
          delta=candidate_flux(p)%z(:,face(1),face(2),face(3))-raw_flux(p)%z(:,face(1),face(2),face(3))
          if(present(internal_energy_rate)) delta_aux=auxiliary_flux(p)%z(1,face(1),face(2),face(3))- &
            auxiliary_raw(p)%z(1,face(1),face(2),face(3))
          face_area=geometry(p)%dx*geometry(p)%dy*geometry(p)%z_face_fraction(face(1),face(2),face(3))
        end select
        theta=1
        do q=1,nv-5
          if(delta(q+5)>0) then
            if(side==1) theta=min(theta,alpha%fields(p)%values(q,inside(1),inside(2),inside(3)))
            if(side==0) theta=min(theta,neighbor_alpha(q))
          else if(delta(q+5)<0) then
            if(side==0) theta=min(theta,alpha%fields(p)%values(q,inside(1),inside(2),inside(3)))
            if(side==1) theta=min(theta,neighbor_alpha(q))
          end if
        end do
        if(theta>=1) cycle
        new_flux=theta*delta
        fluid_volume=geometry(p)%volume_fraction(inside(1),inside(2),inside(3))* &
          geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
        if(.not.ieee_is_finite(fluid_volume).or.fluid_volume<=0) then
          found=.false.;exit
        end if
        sign_face=merge(1._dp,-1._dp,side==0)
        rates%fields(p)%values(:,inside(1),inside(2),inside(3))= &
          rates%fields(p)%values(:,inside(1),inside(2),inside(3))+ &
          sign_face*face_area*(new_flux-delta)/(dt*fluid_volume)
        if(present(internal_energy_rate)) auxiliary_rates%fields(p)%values(1,inside(1),inside(2),inside(3))= &
          auxiliary_rates%fields(p)%values(1,inside(1),inside(2),inside(3))+ &
          sign_face*face_area*(theta-1._dp)*delta_aux/(dt*fluid_volume)
        select case(axis)
        case(1)
          candidate_flux(p)%x(:,face(1),face(2),face(3))=raw_flux(p)%x(:,face(1),face(2),face(3))+new_flux
          if(present(internal_energy_rate)) auxiliary_flux(p)%x(1,face(1),face(2),face(3))= &
            auxiliary_raw(p)%x(1,face(1),face(2),face(3))+theta*delta_aux
        case(2)
          candidate_flux(p)%y(:,face(1),face(2),face(3))=raw_flux(p)%y(:,face(1),face(2),face(3))+new_flux
          if(present(internal_energy_rate)) auxiliary_flux(p)%y(1,face(1),face(2),face(3))= &
            auxiliary_raw(p)%y(1,face(1),face(2),face(3))+theta*delta_aux
        case(3)
          candidate_flux(p)%z(:,face(1),face(2),face(3))=raw_flux(p)%z(:,face(1),face(2),face(3))+new_flux
          if(present(internal_energy_rate)) auxiliary_flux(p)%z(1,face(1),face(2),face(3))= &
            auxiliary_raw(p)%z(1,face(1),face(2),face(3))+theta*delta_aux
        end select
      end do
      if(found) then
        do p=1,size(owners)
          if(owners(p)/=rank.or.(depth>=0.and.tree%nodes(p)%level/=depth)) cycle
          if(.not.all(ieee_is_finite(rates%fields(p)%values))) found=.false.
          if(present(internal_energy_rate)) then
            if(.not.all(ieee_is_finite(auxiliary_rates%fields(p)%values))) found=.false.
          endif
        end do
      end if
      call agree(found,success)
    end subroutine limit_boundary_faces
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine compute_mol_mpi_eb_tree_hydro_3d
end module mpi_eb_tree_hydro_3d_mod
