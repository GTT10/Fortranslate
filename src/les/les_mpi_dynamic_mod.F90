module les_mpi_dynamic_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use les_filter_3d_mod, only: les_filter_weights,les_filter_code
  use les_mpi_filter_mod, only: filter_mpi_les_field_3d
  use les_dynamic_smagorinsky_mod, only: dynamic_smagorinsky_quantities,dynamic_smagorinsky_coefficients, &
    compute_dynamic_les_rhs_3d
  implicit none
  private
  public :: compute_mpi_dynamic_les_coefficients_3d,compute_mpi_dynamic_les_rhs_3d
  public :: advance_mpi_dynamic_les_ssprk2_3d
contains
  subroutine advance_mpi_dynamic_les_ssprk2_3d(comm,species,state,temperature,spacing,cells,test_model, &
      test_ratio,dt,ok,grid_ratio,harmonic_mean)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    integer, intent(in) :: cells(3),test_ratio
    character(len=*), intent(in) :: test_model
    logical, optional, intent(in) :: harmonic_mean
    logical, intent(out) :: ok
    integer, optional, intent(in) :: grid_ratio
    real(dp), allocatable :: rhs(:,:,:,:),stage(:,:,:,:),candidate(:,:,:,:),t1(:,:,:),t2(:,:,:)
    real(dp) :: diffusion,dt_min,dt_max
    integer :: n(3),nv,ios,ierr
    logical :: valid,global_valid
    ok=.false.;n=shape(temperature);nv=size(state,1)
    valid=ieee_is_finite(dt)
    call agree()
    if (.not.global_valid) return
    call MPI_Allreduce(dt,dt_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,dt_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (dt_min/=dt_max.or.dt<0.or.dt>1e6_dp) return
    allocate(rhs(nv,size(state,2),size(state,3),size(state,4)),stage(nv,n(1),n(2),n(3)), &
      candidate(nv,n(1),n(2),n(3)),t1(n(1),n(2),n(3)),t2(n(1),n(2),n(3)),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    call compute_mpi_dynamic_les_rhs_3d(comm,species,state,temperature,spacing,cells,test_model,test_ratio, &
      rhs,valid,grid_ratio,diffusion,harmonic_mean)
    call agree()
    if (.not.global_valid) return
    if (dt*diffusion*sum(1/spacing**2)>.2_dp) return
    if (dt==0) then
      ok=.true.;return
    end if
    stage=state+dt*rhs
    call recover_reactive_temperatures_3d(species,stage,temperature,n(1),n(2),n(3),t1,valid)
    call agree()
    if (.not.global_valid) return
    call compute_mpi_dynamic_les_rhs_3d(comm,species,stage,t1,spacing,cells,test_model,test_ratio,rhs,valid, &
      grid_ratio,diffusion,harmonic_mean)
    call agree()
    if (.not.global_valid) return
    if (dt*diffusion*sum(1/spacing**2)>.2_dp) return
    candidate=.5_dp*state+.5_dp*stage+.5_dp*dt*rhs
    call recover_reactive_temperatures_3d(species,candidate,t1,n(1),n(2),n(3),t2,valid)
    call agree()
    if (.not.global_valid) return
    state=candidate;temperature=t2;ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
  end subroutine

  subroutine compute_mpi_dynamic_les_rhs_3d(comm,species,state,temperature,spacing,cells,test_model,test_ratio, &
      rhs,ok,grid_ratio,maximum_diffusivity,harmonic_mean)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3)
    integer, intent(in) :: cells(3),test_ratio
    character(len=*), intent(in) :: test_model
    real(dp), intent(inout) :: rhs(:,:,:,:)
    logical, optional, intent(in) :: harmonic_mean
    logical, intent(out) :: ok
    integer, optional, intent(in) :: grid_ratio
    real(dp), optional, intent(out) :: maximum_diffusivity
    real(dp), allocatable :: coeff(:,:,:,:),packed(:,:,:,:),halo(:,:,:,:),left_stream(:,:,:,:),right_stream(:,:,:,:)
    real(dp), allocatable :: send_left(:,:,:),send_right(:,:,:),recv_left(:,:,:),recv_right(:,:,:),halo_rhs(:,:,:,:)
    real(dp) :: diffusion,global_diffusion
    integer :: n(3),nv,nc,rank,nranks,ierr,ios,count,left,right,s
    integer :: mean_code,mean_min,mean_max
    logical :: valid,global_valid,harmonic
    ok=.false.;n=shape(temperature);nv=size(state,1);nc=nv+4
    if (present(maximum_diffusivity)) maximum_diffusivity=0
    harmonic=.false.
    if (present(harmonic_mean)) harmonic=harmonic_mean
    mean_code=merge(1,0,harmonic)
    call MPI_Allreduce(mean_code,mean_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(mean_code,mean_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.mean_min/=mean_max) return
    valid=all(shape(rhs)==shape(state)).and.all(n>=1)
    call agree()
    if (.not.global_valid) return
    allocate(coeff(3,n(1),n(2),n(3)),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    call compute_mpi_dynamic_les_coefficients_3d(comm,species,state,temperature,spacing,cells,test_model, &
      test_ratio,coeff,valid,grid_ratio)
    call agree()
    if (.not.global_valid) return
    valid=int(nc,int64)*n(2)*n(3)<=huge(count)
    call agree()
    if (.not.global_valid) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    allocate(packed(nc,n(1),n(2),n(3)),left_stream(nc,n(1),n(2),n(3)),right_stream(nc,n(1),n(2),n(3)), &
      halo(nc,n(1)+4,n(2),n(3)),halo_rhs(nv,n(1)+4,n(2),n(3)), &
      send_left(nc,n(2),n(3)),send_right(nc,n(2),n(3)),recv_left(nc,n(2),n(3)),recv_right(nc,n(2),n(3)),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    packed(:nv,:,:,:)=state;packed(nv+1,:,:,:)=temperature;packed(nv+2:,:,:,:)=coeff
    left_stream=packed;right_stream=packed;halo(:,3:n(1)+2,:,:)=packed
    left=modulo(rank-1,nranks);right=modulo(rank+1,nranks);count=nc*n(2)*n(3)
    do s=1,2
      send_left=right_stream(:,1,:,:);send_right=left_stream(:,n(1),:,:)
      call MPI_Sendrecv(send_left,count,MPI_DOUBLE_PRECISION,left,491,recv_right,count,MPI_DOUBLE_PRECISION, &
        right,491,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Sendrecv(send_right,count,MPI_DOUBLE_PRECISION,right,492,recv_left,count,MPI_DOUBLE_PRECISION, &
        left,492,comm,MPI_STATUS_IGNORE,ierr)
      if (ierr/=MPI_SUCCESS) return
      halo(:,3-s,:,:)=recv_left;halo(:,n(1)+2+s,:,:)=recv_right
      left_stream(:,2:n(1),:,:)=left_stream(:,1:n(1)-1,:,:);left_stream(:,1,:,:)=recv_left
      right_stream(:,1:n(1)-1,:,:)=right_stream(:,2:n(1),:,:);right_stream(:,n(1),:,:)=recv_right
    end do
    call compute_dynamic_les_rhs_3d(species,halo(:nv,:,:,:),halo(nv+1,:,:,:),spacing,test_model,test_ratio, &
      halo_rhs,valid,grid_ratio,diffusion,halo(nv+2:,:,:,:),[3,n(1)+2],harmonic)
    call agree()
    if (.not.global_valid) return
    call MPI_Allreduce(diffusion,global_diffusion,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    rhs=halo_rhs(:,3:n(1)+2,:,:)
    if (present(maximum_diffusivity)) maximum_diffusivity=global_diffusion
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
  end subroutine

  subroutine compute_mpi_dynamic_les_coefficients_3d(comm,species,state,temperature,spacing,cells,test_model, &
      test_ratio,coefficients,ok,grid_ratio)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3)
    integer, intent(in) :: cells(3),test_ratio
    character(len=*), intent(in) :: test_model
    real(dp), intent(inout) :: coefficients(:,:,:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: grid_ratio
    real(dp), allocatable :: weights(:),smooth(:),filtered_state(:,:,:,:),q(:,:,:,:),qt(:,:,:,:)
    real(dp), allocatable :: quantities(:,:,:,:),filtered_quantities(:,:,:,:),test_quantities(:,:,:,:),raw(:,:,:,:)
    real(dp) :: low_spacing(3),high_spacing(3)
    integer :: n(3),ratio,model_id,config(7),low(7),high(7),rank,nranks,ierr,ios
    logical :: valid,global_valid
    ok=.false.;n=shape(temperature);ratio=1;model_id=0
    if (present(grid_ratio)) ratio=grid_ratio
    model_id=les_filter_code(test_model)
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    config=[cells,size(species),model_id,test_ratio,ratio]
    call MPI_Allreduce(config,low,7,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,high,7,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(low==high).and.model_id>0.and.ratio>=1.and.ratio<=64.and.all(cells>=2).and.all(n>=1)
    valid=valid.and.size(species)>=1.and.size(state,1)==reactive_nvar(size(species))
    valid=valid.and.all([size(state,2),size(state,3),size(state,4)]==n).and.all(shape(coefficients)==[3,n])
    valid=valid.and.all(ieee_is_finite(spacing))
    call agree()
    if (.not.global_valid) return
    call MPI_Allreduce(spacing,low_spacing,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(spacing,high_spacing,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (any(low_spacing/=high_spacing).or.any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
    call les_filter_weights(test_model,test_ratio,weights,valid)
    call agree()
    if (.not.global_valid) return
    call les_filter_weights('box',6,smooth,valid)
    call agree()
    if (.not.global_valid) return
    allocate(filtered_state(size(state,1),n(1),n(2),n(3)),q(5,n(1),n(2),n(3)),qt(5,n(1),n(2),n(3)), &
      quantities(25,n(1),n(2),n(3)),filtered_quantities(25,n(1),n(2),n(3)), &
      test_quantities(25,n(1),n(2),n(3)),raw(3,n(1),n(2),n(3)),stat=ios)
    valid=ios==0;call agree()
    if (.not.global_valid) return
    call filter_mpi_les_field_3d(comm,state,cells,weights,filtered_state,valid)
    call agree()
    if (.not.global_valid) return
    call recover_fields()
    call agree()
    if (.not.global_valid) return
    call form_quantities(q,ratio,quantities)
    call agree()
    if (.not.global_valid) return
    call form_quantities(qt,test_ratio,test_quantities)
    call agree()
    if (.not.global_valid) return
    call filter_mpi_les_field_3d(comm,quantities,cells,weights,filtered_quantities,valid)
    call agree()
    if (.not.global_valid) return
    call contract_fields()
    call agree()
    if (.not.global_valid) return
    call filter_mpi_les_field_3d(comm,raw,cells,smooth,coefficients,ok)
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine

    subroutine recover_fields()
      real(dp) :: primitive(reactive_nprim(size(species))),t,sound
      integer :: i,j,k
      valid=.false.
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),primitive,t,sound,valid)
            if (.not.valid) return
            valid=maxval(abs(primitive(2:4)))<=1e30_dp.and.primitive(1)<=1e8_dp.and.t<=1e8_dp
            if (.not.valid) return
            q(:,i,j,k)=[primitive(1:4),t]
            call reactive_conserved_to_primitive(species,filtered_state(:,i,j,k),temperature(i,j,k),primitive,t,sound,valid)
            if (.not.valid) return
            valid=maxval(abs(primitive(2:4)))<=1e30_dp.and.primitive(1)<=1e8_dp.and.t<=1e8_dp
            if (.not.valid) return
            qt(:,i,j,k)=[primitive(1:4),t]
          end do
        end do
      end do
    end subroutine

    subroutine form_quantities(primitive,width,output)
      real(dp), intent(in) :: primitive(:,:,:,:)
      integer, intent(in) :: width
      real(dp), intent(out) :: output(:,:,:,:)
      real(dp), allocatable :: halo(:,:,:,:),send_left(:,:,:),send_right(:,:,:),recv_left(:,:,:),recv_right(:,:,:)
      real(dp) :: g(3,3),gt(3),kij(3,3),rut(3),aij(3,3),a,ft(3)
      integer :: i,j,k,d,lo(3),hi(3),cell(3),left,right,count
      allocate(halo(5,n(1)+2,n(2),n(3)),send_left(5,n(2),n(3)),send_right(5,n(2),n(3)), &
        recv_left(5,n(2),n(3)),recv_right(5,n(2),n(3)),stat=ios)
      valid=ios==0;call agree()
      if (.not.global_valid) then
        valid=.false.;return
      end if
      left=modulo(rank-1,nranks);right=modulo(rank+1,nranks);count=5*n(2)*n(3)
      send_left=primitive(:,1,:,:);send_right=primitive(:,n(1),:,:)
      call MPI_Sendrecv(send_left,count,MPI_DOUBLE_PRECISION,left,481,recv_right,count,MPI_DOUBLE_PRECISION, &
        right,481,comm,MPI_STATUS_IGNORE,ierr)
      valid=ierr==MPI_SUCCESS
      if (.not.valid) return
      call MPI_Sendrecv(send_right,count,MPI_DOUBLE_PRECISION,right,482,recv_left,count,MPI_DOUBLE_PRECISION, &
        left,482,comm,MPI_STATUS_IGNORE,ierr)
      valid=ierr==MPI_SUCCESS
      if (.not.valid) return
      halo(:,2:n(1)+1,:,:)=primitive;halo(:,1,:,:)=recv_left;halo(:,n(1)+2,:,:)=recv_right
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            cell=[i+1,j,k]
            do d=1,3
              lo=cell;hi=cell
              if (d==1) then
                lo(d)=cell(d)-1;hi(d)=cell(d)+1
              else
                lo(d)=modulo(cell(d)-2,n(d))+1;hi(d)=modulo(cell(d),n(d))+1
              end if
              g(:,d)=(halo(2:4,hi(1),hi(2),hi(3))-halo(2:4,lo(1),lo(2),lo(3)))/(2*spacing(d))
              gt(d)=(halo(5,hi(1),hi(2),hi(3))-halo(5,lo(1),lo(2),lo(3)))/(2*spacing(d))
            end do
            call dynamic_smagorinsky_quantities(g,gt,primitive(1,i,j,k),primitive(2:4,i,j,k),primitive(5,i,j,k), &
              spacing,width,kij,rut,aij,a,ft,valid)
            if (.not.valid) return
            output(:,i,j,k)=[reshape(kij,[9]),rut,reshape(aij,[9]),a,ft]
          end do
        end do
      end do
    end subroutine

    subroutine contract_fields()
      integer :: i,j,k
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            call dynamic_smagorinsky_coefficients(qt(1,i,j,k),qt(2:4,i,j,k),qt(5,i,j,k), &
              reshape(filtered_quantities(1:9,i,j,k),[3,3]),filtered_quantities(10:12,i,j,k), &
              reshape(filtered_quantities(13:21,i,j,k),[3,3]),filtered_quantities(22,i,j,k), &
              filtered_quantities(23:25,i,j,k),reshape(test_quantities(13:21,i,j,k),[3,3]), &
              test_quantities(22,i,j,k),test_quantities(23:25,i,j,k),raw(:,i,j,k),valid)
            if (.not.valid) return
          end do
        end do
      end do
    end subroutine
  end subroutine
end module
