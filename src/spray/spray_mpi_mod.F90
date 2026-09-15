module spray_mpi_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use spray_parcel_mod
  use spray_mixture_mod, only: valid_component_inventory,liquid_mixture_properties
  use spray_breakup_mod, only: split_spray_parcels
  use spray_partition_mod
  implicit none
  private
  public :: migrate_spray_parcels,reserve_spray_ids,split_mpi_spray_parcels,inject_mpi_cone_parcels
contains
  subroutine inject_mpi_cone_parcels(comm,liquid,total_mass,number,center,axis,speed,half_angle,diameter,temperature, &
      maximum_parcels,last_id,parcels,injected_totals,ok,components,fractions,component_mass, &
      diameter_shape,diameter_min,diameter_max,diameter_weighting,diameter_distribution,diameter_std,nozzle_radius,nozzle_model)
    type(MPI_Comm), intent(in) :: comm
    type(spray_liquid), intent(in) :: liquid
    real(dp), intent(in) :: total_mass,center(3),axis(3),speed,half_angle,diameter,temperature
    integer, intent(in) :: number,maximum_parcels
    integer(int64), intent(inout) :: last_id
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(out) :: injected_totals(:)
    logical, intent(out) :: ok
    type(spray_liquid), optional, intent(in) :: components(:)
    real(dp), optional, intent(in) :: fractions(:)
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    real(dp), optional, intent(in) :: diameter_shape,diameter_min,diameter_max,diameter_std,nozzle_radius
    character(len=*), optional, intent(in) :: diameter_weighting,diameter_distribution
    type(spray_nozzle_model), optional, intent(in) :: nozzle_model
    type(spray_parcel), allocatable :: candidate(:),born(:)
    type(spray_liquid) :: material
    real(dp), allocatable :: masses(:,:)
    real(dp) :: rho,cp,e,ledger(size(injected_totals)),inventory(5)
    integer :: p,c,nc,lo,hi,old_count,ierr
    integer(int64) :: next_id,first_id
    logical :: valid,agreed,mixed
    ok=.false.;injected_totals=0;nc=0;mixed=present(components)
    valid=allocated(parcels).and.maximum_parcels>=1.and.maximum_parcels<=1000000
    valid=valid.and.number>=0.and.number<=1000000.and.last_id>=0.and.size(injected_totals)>=6
    valid=valid.and.(mixed.eqv.present(fractions)).and.(mixed.eqv.present(component_mass))
    if (present(component_mass)) valid=valid.and.allocated(component_mass)
    if (.not.ieee_is_finite(total_mass)) valid=.false.
    if (number==0.and.total_mass/=0) valid=.false.
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    old_count=size(parcels);valid=number<=maximum_parcels-old_count
    material=liquid
    if (mixed) then
      nc=size(components)
      valid=valid.and.nc>=1.and.nc<=16.and.size(fractions)==nc
      valid=valid.and.all(shape(component_mass)==[nc,old_count])
      if (valid) then
        call liquid_mixture_properties(components,fractions,temperature,rho,cp,e,valid)
        if (valid) then
          material=components(1);material%density=rho;material%heat_capacity=cp
          material%reference_temperature=temperature;material%reference_energy=e
          material%density_model='constant';material%density_coefficients=0;material%density_temperature_range=0
          material%heat_capacity_model='constant';material%heat_capacity_coefficients=0
          material%heat_capacity_temperature_range=0;material%boiling_model='reject';material%critical_temperature=0
        end if
      end if
    else
      valid=valid.and.valid_spray_liquid(material)
    end if
    if (valid) then
      do p=1,old_count
        valid=valid.and.valid_spray_parcel(parcels(p)).and.parcels(p)%id<=last_id
        if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
        if (mixed) valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      end do
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(nc,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.lo/=hi) return
    next_id=last_id;first_id=0
    call reserve_spray_ids(comm,number,next_id,first_id,valid)
    if (.not.valid) return
    allocate(born(0));ledger=0
    if (number>0) call inject_cone_parcels(material,total_mass,number,first_id,center,axis,speed,half_angle,diameter, &
      temperature,born,valid,diameter_shape,diameter_min,diameter_max,diameter_weighting,diameter_distribution, &
      diameter_std,nozzle_radius,nozzle_model)
    if (valid) then
      do p=1,size(born)
        inventory=parcel_totals(born(p),material);ledger(1:5)=ledger(1:5)+inventory
        if (mixed) then
          do c=1,nc
            if (components(c)%vapor_index>size(ledger)-5) then
              valid=.false.;exit
            end if
            ledger(5+components(c)%vapor_index)=ledger(5+components(c)%vapor_index)+inventory(1)*fractions(c)
          end do
        else
          if (material%vapor_index>size(ledger)-5) then
            valid=.false.;exit
          end if
          ledger(5+material%vapor_index)=ledger(5+material%vapor_index)+inventory(1)
        end if
      end do
    end if
    valid=valid.and.all(ieee_is_finite(ledger))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    candidate=[parcels,born]
    if (mixed) then
      allocate(masses(nc,size(candidate)));masses(:,:old_count)=component_mass
      do p=1,size(born)
        masses(:,old_count+p)=born(p)%mass*fractions
      end do
      call move_alloc(masses,component_mass)
    end if
    call move_alloc(candidate,parcels);last_id=next_id;injected_totals=ledger;ok=.true.
  end subroutine

  subroutine split_mpi_spray_parcels(comm,parcels,maximum_multiplicity,maximum_parcels,last_id,ok,component_mass)
    type(MPI_Comm), intent(in) :: comm
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: maximum_multiplicity
    integer, intent(in) :: maximum_parcels
    integer(int64), intent(inout) :: last_id
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: candidate(:)
    real(dp), allocatable :: masses(:,:)
    integer :: p,parent,old_count,nc,lo,hi,ierr
    integer(int64) :: next_id,first_id
    logical :: valid,agreed
    ok=.false.;valid=allocated(parcels).and.last_id>=0;nc=0
    if (present(component_mass)) then
      valid=valid.and.allocated(component_mass)
      if (allocated(component_mass)) nc=size(component_mass,1)
      valid=valid.and.nc>=1.and.nc<=16
    end if
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(nc,lo,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(nc,hi,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.lo/=hi) return
    old_count=size(parcels);candidate=parcels
    if (present(component_mass)) valid=size(component_mass,2)==old_count
    if (valid) then
      do p=1,old_count
        valid=valid.and.valid_spray_parcel(parcels(p)).and.parcels(p)%id<=last_id
        if (p>1) valid=valid.and.parcels(p)%id>parcels(p-1)%id
        if (present(component_mass)) valid=valid.and.valid_component_inventory(parcels(p),component_mass(:,p))
      end do
    end if
    if (valid) call split_spray_parcels(candidate,maximum_multiplicity,maximum_parcels,valid)
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    if (present(component_mass)) then
      allocate(masses(nc,size(candidate)));masses(:,:old_count)=component_mass;parent=1
      do p=old_count+1,size(candidate)
        do while(parent<old_count)
          if (parcels(parent)%id>=candidate(p)%parent_id) exit
          parent=parent+1
        end do
        masses(:,p)=component_mass(:,parent)
      end do
    end if
    next_id=last_id;first_id=0
    call reserve_spray_ids(comm,size(candidate)-old_count,next_id,first_id,valid)
    if (.not.valid) return
    do p=old_count+1,size(candidate)
      candidate(p)%id=first_id+int(p-old_count-1,int64)
    end do
    call move_alloc(candidate,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    last_id=next_id;ok=.true.
  end subroutine

  ! All ranks retain the same high-water mark; rank order breaks allocation ties.
  subroutine reserve_spray_ids(comm,local_count,last_id,first_id,ok)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: local_count
    integer(int64), intent(inout) :: last_id,first_id
    logical, intent(out) :: ok
    integer(int64) :: minimum,maximum,requested,total,prefix,next_first
    integer :: rank,ierr
    logical :: valid,agreed
    ok=.false.;valid=local_count>=0.and.local_count<=1000000.and.last_id>=0
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
    if (ierr/=MPI_SUCCESS.or..not.agreed) return
    call MPI_Allreduce(last_id,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(last_id,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS.or.minimum/=maximum) return
    requested=int(local_count,int64)
    call MPI_Allreduce(requested,total,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (total>huge(last_id)-last_id) return
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Exscan(requested,prefix,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (rank==0) prefix=0
    next_first=0
    if (local_count>0) next_first=last_id+prefix+1_int64
    first_id=next_first;last_id=last_id+total;ok=.true.
  end subroutine

  subroutine migrate_spray_parcels(comm,parcels,lower,spacing,cells,parts,periodic,maximum_parcels,ok,component_mass)
    type(MPI_Comm), intent(in) :: comm
    type(spray_parcel), allocatable, intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3)
    integer, intent(in) :: cells(3),parts(3),maximum_parcels
    logical, intent(in) :: periodic(3)
    logical, intent(out) :: ok
    real(dp), allocatable, optional, intent(inout) :: component_mass(:,:)
    type(spray_parcel), allocatable :: outgoing(:),incoming(:)
    real(dp), allocatable :: masses(:,:),send_real(:,:),recv_real(:,:)
    integer(int64), allocatable :: send_int(:,:),recv_int(:,:)
    integer, allocatable :: counts(:),offsets(:),received(:),displacements(:)
    integer :: rank,nranks,ierr,nc,nreal,nrecv,p,r,config(10),lo(10),hi(10)
    real(dp) :: bounds(6),minbounds(6),maxbounds(6)
    logical :: valid,global_valid
    ok=.false.
    call MPI_Comm_rank(comm,rank,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Comm_size(comm,nranks,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=allocated(parcels).and.maximum_parcels>=0.and.maximum_parcels<=1000000;nc=0
    if (allocated(parcels)) valid=valid.and.size(parcels)<=1000000
    if (present(component_mass)) then
      valid=valid.and.allocated(component_mass)
      if (allocated(component_mass)) nc=size(component_mass,1)
      valid=valid.and.nc>=1.and.nc<=16
    end if
    valid=valid.and.all(ieee_is_finite([lower,spacing]))
    call agree()
    if (.not.global_valid) return
    config=[cells,parts,merge(1,0,periodic),nc];bounds=[lower,spacing]
    call MPI_Allreduce(config,lo,10,MPI_INTEGER,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(config,hi,10,MPI_INTEGER,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,minbounds,6,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Allreduce(bounds,maxbounds,6,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=all(lo==hi).and.all(minbounds==maxbounds)
    valid=valid.and.spray_cartesian_owner(lower,lower,spacing,cells,parts,periodic)>=0
    call agree()
    if (.not.global_valid) return
    valid=product(parts)==nranks
    outgoing=parcels
    if (valid) then
      if (present(component_mass)) then
        masses=component_mass
        call partition_spray_parcels(outgoing,lower,spacing,cells,parts,periodic,counts,offsets,valid,masses)
      else
        call partition_spray_parcels(outgoing,lower,spacing,cells,parts,periodic,counts,offsets,valid)
      end if
    end if
    call agree()
    if (.not.global_valid) return
    allocate(received(nranks),displacements(nranks))
    call MPI_Alltoall(counts,1,MPI_INTEGER,received,1,MPI_INTEGER,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    valid=sum(int(received,int64))<=int(maximum_parcels,int64)
    call agree()
    if (.not.global_valid) return
    nrecv=sum(received);displacements(1)=0
    do r=2,nranks
      displacements(r)=displacements(r-1)+received(r-1)
    end do
    nreal=14+nc
    allocate(send_int(3,size(outgoing)),send_real(nreal,size(outgoing)),recv_int(3,nrecv),recv_real(nreal,nrecv))
    do p=1,size(outgoing)
      send_int(:,p)=[outgoing(p)%id,outgoing(p)%parent_id,int(outgoing(p)%material_id,int64)]
      send_real(1:14,p)=[outgoing(p)%position,outgoing(p)%velocity,outgoing(p)%mass,outgoing(p)%temperature, &
        outgoing(p)%multiplicity,outgoing(p)%deformation,outgoing(p)%deformation_rate,outgoing(p)%rt_time, &
        outgoing(p)%khrt_initial_count,outgoing(p)%khrt_shed_volume]
      if (nc>0) send_real(15:,p)=masses(:,p)
    end do
    call MPI_Alltoallv(send_int,3*counts,3*offsets,MPI_INTEGER8,recv_int,3*received,3*displacements,MPI_INTEGER8,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    call MPI_Alltoallv(send_real,nreal*counts,nreal*offsets,MPI_DOUBLE_PRECISION, &
      recv_real,nreal*received,nreal*displacements,MPI_DOUBLE_PRECISION,comm,ierr)
    if (ierr/=MPI_SUCCESS) return
    if (nrecv>1) call sort_records(1,nrecv)
    allocate(incoming(nrecv))
    if (allocated(masses)) deallocate(masses)
    if (nc>0) allocate(masses(nc,nrecv))
    valid=.true.
    do p=1,nrecv
      incoming(p)=spray_parcel(recv_int(1,p),recv_real(1:3,p),recv_real(4:6,p),recv_real(7,p), &
        recv_real(8,p),recv_real(9,p),recv_real(10,p),recv_real(11,p),recv_int(2,p),recv_real(12,p), &
        int(recv_int(3,p)),recv_real(13,p),recv_real(14,p))
      valid=valid.and.valid_spray_parcel(incoming(p))
      if (p>1) valid=valid.and.incoming(p)%id>incoming(p-1)%id
      valid=valid.and.spray_cartesian_owner(incoming(p)%position,lower,spacing,cells,parts,periodic)==rank
      if (nc>0) then
        masses(:,p)=recv_real(15:,p)
        valid=valid.and.valid_component_inventory(incoming(p),masses(:,p))
      end if
    end do
    call agree()
    if (.not.global_valid) return
    call move_alloc(incoming,parcels)
    if (present(component_mass)) call move_alloc(masses,component_mass)
    ok=.true.
  contains
    subroutine agree()
      call MPI_Allreduce(valid,global_valid,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) global_valid=.false.
    end subroutine
    recursive subroutine sort_records(first,last)
      integer, intent(in) :: first,last
      integer :: left,right
      integer(int64) :: pivot,tmp_int(3)
      real(dp) :: tmp_real(nreal)
      left=first;right=last;pivot=recv_int(1,first+(last-first)/2)
      do
        do while(recv_int(1,left)<pivot)
          left=left+1
        end do
        do while(recv_int(1,right)>pivot)
          right=right-1
        end do
        if (left<=right) then
          tmp_int=recv_int(:,left);recv_int(:,left)=recv_int(:,right);recv_int(:,right)=tmp_int
          tmp_real=recv_real(:,left);recv_real(:,left)=recv_real(:,right);recv_real(:,right)=tmp_real
          left=left+1;right=right-1
        end if
        if (left>right) exit
      end do
      if (first<right) call sort_records(first,right)
      if (left<last) call sort_records(left,last)
    end subroutine
  end subroutine
end module
