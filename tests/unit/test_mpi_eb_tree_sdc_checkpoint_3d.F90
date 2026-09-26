program test_mpi_eb_tree_sdc_checkpoint_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d,react_sdc_eb_tree_3d
  use eb_reactive_patch_tree_3d_mod, only: advance_reactive_eb_tree_sdc_3d,build_eb_tree_extended_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d,react_sdc_mpi_eb_tree_3d,advance_mpi_eb_tree_sdc_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,reactive_eb_tree_checkpoint_record_3d
  use mpi_eb_tree_checkpoint_3d_mod, only: write_mpi_eb_tree_checkpoint_3d
  use mpi_eb_tree_checkpoint_read_3d_mod, only: read_mpi_eb_tree_checkpoint_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:),wall(:),chemical(:)
  real(dp) :: residual,theta,speed
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,mode,integrator_index,rank,nranks,status
  logical :: ok,chemistry_on,transport_on
  character(len=8), parameter :: integrators(2)=[character(len=8)::'explicit','implicit']
    type(eb_tree_checkpoint_3d) :: record,restored,serial_record
    type(amr_tree_state_3d) :: local,lt
    type(eb_geometry_3d) :: geometry(5)
    integer :: owners(5),p,stream,checkpoint_mode,layout
    character(len=512) :: message
    character(len=*), parameter :: path='mpi_eb_sdc_checkpoint.dat',serial_path='serial_eb_sdc_checkpoint.dat'
    logical :: good
    integer :: source_count,line_index,io_code
    integer(int64) :: row_position
    character(len=4096) :: line
    real(dp), allocatable :: bad_row(:)
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  allocate(wall(nv),chemical(nv))
  do axis=3,3
  do count=5,5
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
    tree%nodes(node)%state(5,:,:,:)=1.1_dp*tree%nodes(node)%state(5,:,:,:)
    do k=1,4
    do j=1,4
    do i=1,4
      if (tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
      call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k),t(i,j,k),primitive, &
        tree%nodes(node)%temperature(i,j,k),speed,ok)
      call require(ok,'temperature perturbation EOS')
    end do
    end do
    end do
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  call check_checkpoint()
  deallocate(u,t)
  end do
  end do
  if(rank==0) print *, "MPI EB SDC checkpoint: PASS"
  call MPI_Finalize(status)
contains
  subroutine check_checkpoint()
    config%time_integrator='sdc';config%external_forcing=[1._dp,2._dp,3._dp]
    config%redistribution='none';config%plane_axis='z';config%plane_position=.31_dp
    record%forcing=1.e-8_dp+0*before
    record%initial=before-record%forcing;record%scale=abs(before)
    record%boundary=0*before;record%wall=0*before;record%chemical=0*before
    record%time=1.e-8_dp;record%steps=1;record%minimum_dt=1.e-8_dp;record%minimum_theta=1
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
    end do
    do checkpoint_mode=1,2
      config%chemistry_enabled=checkpoint_mode==1
      owners=[(mod(p-1,nranks),p=1,count)]
      call local_fields()
      if(checkpoint_mode==1) then
        allocate(record%reaction_source%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          allocate(record%reaction_source%fields(p)%values(nv-4,4,4,4),source=-real(p,dp))
        end do
      else
        deallocate(record%reaction_source%fields)
      end if
      open(newunit=stream,status='scratch',form='formatted')
      call write_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
        geometry,local,lt,record,MPI_COMM_WORLD,good)
      close(stream)
      call require(good,'MPI V4 write')
      do layout=1,2
        if(layout==2) owners=0
        call local_fields()
        do p=1,count
          if(owners(p)==rank) local%fields(p)%values=.9_dp*local%fields(p)%values
        end do
        open(newunit=stream,file=path,access='stream',form='formatted',status='old')
        call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
          geometry,local,lt,restored,MPI_COMM_WORLD,good)
        close(stream)
        call require(good,'MPI V4 read with changed ownership')
        call check_fields()
      end do
      tree=saved
      open(newunit=stream,file=path,access='stream',form='formatted',status='old')
      call reactive_eb_tree_checkpoint_record_3d(stream,.true.,species,reactions,transport,config, &
        tree,serial_record,good,message)
      close(stream)
      call require(good,'serial reads V4: '//trim(message))
      call require(all(serial_record%forcing==record%forcing),'serial forcing exact')
      do p=1,count
        call require(all(tree%nodes(p)%state==saved%nodes(p)%state),'serial V4 gas exact')
        if(checkpoint_mode==1) call require(all(serial_record%reaction_source%fields(p)%values==-real(p,dp)),'serial V4 IR exact')
      end do
      ! Each process reads the two-rank file using a one-rank communicator.
      owners=0
      if(allocated(local%fields)) deallocate(local%fields,lt%fields)
      allocate(local%fields(count),lt%fields(count))
      do p=1,count
        local%fields(p)%values=saved%nodes(p)%state
        allocate(lt%fields(p)%values(1,4,4,4))
        lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
      end do
      open(newunit=stream,file=path,access='stream',form='formatted',status='old')
      call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
        geometry,local,lt,restored,MPI_COMM_SELF,good)
      close(stream)
      call require(good,'two-rank file read with one rank')
      do p=1,count
        call require(all(local%fields(p)%values==saved%nodes(p)%state),'one-rank gas exact')
        if(checkpoint_mode==1) call require(all(restored%reaction_source%fields(p)%values==-real(p,dp)),'one-rank IR exact')
      end do
      if(rank==0) then
        open(newunit=stream,file=serial_path,access='stream',form='formatted',status='replace')
        call reactive_eb_tree_checkpoint_record_3d(stream,.false.,species,reactions,transport,config, &
          tree,serial_record,good,message)
        close(stream)
        call require(good,'serial V3 write')
      end if
      call MPI_Barrier(MPI_COMM_WORLD,status)
      owners=[(mod(p-1,nranks),p=1,count)];call local_fields()
      open(newunit=stream,file=serial_path,access='stream',form='formatted',status='old')
      call read_mpi_eb_tree_checkpoint_3d(serial_path,stream,species,reactions,transport,config,topology,owners, &
        geometry,local,lt,restored,MPI_COMM_WORLD,good)
      close(stream)
      call require(good,'MPI reads serial V3')
      call check_fields()
      restored%time=17
      config%external_forcing(1)=2
      open(newunit=stream,file=path,access='stream',form='formatted',status='old')
      call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
        geometry,local,lt,restored,MPI_COMM_WORLD,good)
      close(stream)
      config%external_forcing(1)=1
      call require(.not.good.and.restored%time==17,'mismatched force rejects without publication')
      call check_fields()
      if(checkpoint_mode==1) then
        if(rank==0) then
          open(newunit=stream,file=path,access='stream',form='formatted',status='old',action='readwrite')
          source_count=0
          do
            read(stream,'(a)',iostat=io_code) line
            call require(io_code==0,'locate final source record')
            if(trim(line)/='SDC_REACTION_SOURCE') cycle
            source_count=source_count+1
            if(source_count==count) exit
          end do
          do line_index=1,63
            read(stream,'(a)') line
          end do
          inquire(stream,pos=row_position)
          allocate(bad_row(nv-4),source=ieee_value(0._dp,ieee_quiet_nan))
          write(stream,'(*(es27.18e3,1x))',pos=row_position) bad_row
          write(stream,'(a)') 'END_EB_TREE_CHECKPOINT'
          deallocate(bad_row)
          close(stream)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,status)
        do p=1,count
          if(owners(p)==rank) local%fields(p)%values=.9_dp*local%fields(p)%values
        end do
        open(newunit=stream,file=path,access='stream',form='formatted',status='old')
        call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
          geometry,local,lt,restored,MPI_COMM_WORLD,good)
        close(stream)
        call require(.not.good.and.restored%time==17,'late corrupt source rejects without publishing record')
        do p=1,count
          if(owners(p)/=rank) cycle
          call require(all(local%fields(p)%values==.9_dp*saved%nodes(p)%state),'late failure gas rollback')
          call require(all(lt%fields(p)%values(1,:,:,:)==saved%nodes(p)%temperature),'late failure T rollback')
          call require(all(restored%reaction_source%fields(p)%values==-real(p,dp)),'late failure IR rollback')
        end do
      end if
      call MPI_Barrier(MPI_COMM_WORLD,status)
    end do
  end subroutine check_checkpoint
    subroutine local_fields()
      integer :: q
      if(allocated(local%fields)) deallocate(local%fields,lt%fields)
      allocate(local%fields(count),lt%fields(count))
      do q=1,count
        if(owners(q)/=rank) cycle
        local%fields(q)%values=saved%nodes(q)%state
        allocate(lt%fields(q)%values(1,4,4,4))
        lt%fields(q)%values(1,:,:,:)=saved%nodes(q)%temperature
      end do
    end subroutine
    subroutine check_fields()
      integer :: q
      call require(all(restored%forcing==record%forcing),'MPI forcing exact')
      do q=1,count
        call require(allocated(local%fields(q)%values).eqv.(owners(q)==rank),'owner-local restart gas')
        if(checkpoint_mode==1) call require(allocated(restored%reaction_source%fields(q)%values).eqv.(owners(q)==rank), &
          'owner-local restart IR')
        if(owners(q)/=rank) cycle
        call require(all(local%fields(q)%values==saved%nodes(q)%state),'MPI gas exact')
        call require(all(lt%fields(q)%values(1,:,:,:)==saved%nodes(q)%temperature),'MPI T exact')
        if(checkpoint_mode==1) call require(all(restored%reaction_source%fields(q)%values==-real(q,dp)),'MPI IR exact')
      end do
      if(checkpoint_mode==2) call require(.not.allocated(restored%reaction_source%fields),'chemistry off no IR')
    end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if(condition) return
    print *, 'FAIL: ',label
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
end program test_mpi_eb_tree_sdc_checkpoint_3d
