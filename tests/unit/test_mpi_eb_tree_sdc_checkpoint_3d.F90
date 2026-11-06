program test_mpi_eb_tree_sdc_checkpoint_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use dual_energy_mod, only: initialize_nasa7_internal_energy_patch
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
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
    character(len=32) :: checkpoint_action
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call get_command_argument(1,checkpoint_action)
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
  if(len_trim(checkpoint_action)==0) call check_checkpoint()
  call check_held_checkpoint()
  deallocate(u,t)
  end do
  end do
  if(rank==0) print *, "MPI EB SDC checkpoint: PASS"
  call MPI_Finalize(status)
contains
  subroutine check_held_checkpoint()
    type(eb_reactive_tree_3d) :: seed,serial_tree
    type(amr_tree_state_3d) :: seed_h,gas,temp,kept_gas,kept_temp
    type(eb_tree_checkpoint_3d) :: written,loaded,serial_record,kept_record
    type(rotating_frame_3d), allocatable :: rotation
    type(reactive_eb_3d_config) :: cfg
    real(dp) :: inventory(nv),rho,phi,position(3)
    integer :: condition,mode,layout,owner(count),p,x,y,z,stream,code,scalar_row
    integer(int64) :: final_scalar,final_metadata
    character(len=128) :: mpi_path,serial_path
    character(len=4096) :: text
    logical :: good,write_only,read_only
    write_only=checkpoint_action=='held_write';read_only=checkpoint_action=='held_read'
    call require(len_trim(checkpoint_action)==0.or.write_only.or.read_only,'known held checkpoint diagnostic phase')
    do condition=1,2
      if(condition==2) then
        allocate(rotation);rotation%axis=2;rotation%angular_speed=1.e10_dp;rotation%axis_point=[10._dp,10._dp,10._dp]
      endif
      seed=saved
      allocate(seed_h%fields(count))
      do p=1,count
        allocate(seed_h%fields(p)%values(1,4,4,4),source=27._dp)
        where(seed%nodes(p)%geometry%volume_fraction>0) seed%nodes(p)%temperature=merge(300._dp,700._dp,condition==1)
        call initialize_nasa7_internal_energy_patch(species,seed%nodes(p)%state,seed%nodes(p)%temperature, &
          seed%nodes(p)%geometry%volume_fraction,seed_h%fields(p)%values(1,:,:,:),good)
        call require(good,'held tree checkpoint independent native thermal seed')
        do z=1,4
        do y=1,4
        do x=1,4
          if(seed%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          rho=seed%nodes(p)%state(1,x,y,z)
          seed%nodes(p)%state(2:4,x,y,z)=[rho*1.e10_dp,0._dp,0._dp];seed%nodes(p)%state(5,x,y,z)=.5_dp*rho*1.e20_dp
          if(allocated(rotation)) then
            associate(g=>seed%nodes(p)%geometry)
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_potential_3d(rotation,position,phi,good)
            call require(good,'held tree checkpoint rotating gas representation')
            seed%nodes(p)%state(5,x,y,z)=seed%nodes(p)%state(5,x,y,z)-rho*phi
          endif
        enddo
        enddo
        enddo
        geometry(p)=seed%nodes(p)%geometry
      enddo
      call average_down_reactive_eb_tree_3d(species,seed,good,rotation,seed_h)
      call require(good,'paired held checkpoint coarse consistency')
      call integrate_reactive_eb_tree_3d(seed,inventory,good)
      call require(good,'independent held checkpoint composite inventory')
      do mode=1,3
        cfg=config;cfg%time_integrator='split';cfg%external_forcing=0;cfg%chemistry_enabled=.false.
        cfg%redistribution='none';cfg%plane_axis='z';cfg%plane_position=.31_dp
        if(mode>=2) then
          cfg%time_integrator='sdc';cfg%external_forcing=[1._dp,2._dp,3._dp];cfg%chemistry_enabled=.true.
        endif
        write(mpi_path,'(a,i0,a,i0,a)') 'held_mpi_c',condition,'_m',mode,'.dat'
        write(serial_path,'(a,i0,a,i0,a)') 'held_serial_c',condition,'_m',mode,'.dat'
        written=eb_tree_checkpoint_3d();written%energy_correction=real((-1)**condition,dp)*1.e-3_dp*abs(inventory(5))
        written%energy_policy%eta2=2;written%energy_policy%allow_negative_energy=condition==1
        written%initial=inventory;written%initial(5)=written%initial(5)-written%energy_correction
        written%scale=abs(inventory);written%boundary=0*inventory;written%wall=0*inventory;written%chemical=0*inventory
        written%time=1.e-8_dp;written%steps=1;written%minimum_dt=1.e-8_dp;written%minimum_theta=1
        if(mode>=2) then
          written%forcing=1.e-8_dp+0*inventory;written%initial=written%initial-written%forcing
        endif
        owner=[(mod(p-1,nranks),p=1,count)]
        allocate(gas%fields(count),temp%fields(count),written%internal_energy%fields(count))
        if(mode>=2) allocate(written%reaction_source%fields(count))
        if(mode==3) allocate(written%heat_release%fields(count))
        do p=1,count
          if(owner(p)/=rank) cycle
          gas%fields(p)%values=seed%nodes(p)%state
          allocate(temp%fields(p)%values(1,4,4,4));temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          written%internal_energy%fields(p)=seed_h%fields(p)
          if(mode>=2) allocate(written%reaction_source%fields(p)%values(nv-4,4,4,4),source=-real(p,dp))
          if(mode==3) allocate(written%heat_release%fields(p)%values(1,4,4,4),source=(-1._dp)**p*real(p,dp)*1.e9_dp)
        enddo
        if(.not.read_only) then
          open(newunit=stream,status='scratch',form='formatted')
          call write_mpi_eb_tree_checkpoint_3d(trim(mpi_path),stream,species,reactions,transport,cfg,topology,owner, &
            geometry,gas,temp,written,MPI_COMM_WORLD,good,frame=rotation)
          close(stream);call require(good,'MPI indexed held tree checkpoint write with one composite prefix')
        endif
        deallocate(gas%fields,temp%fields)
        if(write_only) cycle
        do layout=1,2
          owner=[(mod(p+nranks-2,nranks),p=1,count)];if(layout==2) owner=0
          allocate(gas%fields(count),temp%fields(count))
          loaded=eb_tree_checkpoint_3d();loaded%energy_policy=written%energy_policy
          allocate(loaded%internal_energy%fields(count))
          loaded%time=17;loaded%energy_correction=19
          do p=1,count
            if(owner(p)/=rank) cycle
            gas%fields(p)%values=.9_dp*seed%nodes(p)%state
            allocate(temp%fields(p)%values(1,4,4,4));temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
            loaded%internal_energy%fields(p)%values=.9_dp*seed_h%fields(p)%values
          enddo
          open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old')
          call read_mpi_eb_tree_checkpoint_3d(trim(mpi_path),stream,species,reactions,transport,cfg,topology,owner, &
            geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
          close(stream);call require(good,'MPI held read across ownership and rank count')
          call require(loaded%energy_correction==written%energy_correction.and.loaded%time==written%time, &
            'signed composite held prefix and time recovered exactly once')
          do p=1,count
            if(owner(p)/=rank) then
              call require(.not.allocated(loaded%internal_energy%fields(p)%values),'idle restart rank has no UEINT payload')
              cycle
            endif
            call require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
              all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature).and. &
              all(loaded%internal_energy%fields(p)%values==seed_h%fields(p)%values),'exact changed-owner gas/T/UEINT restart')
            if(mode>=2) call require(all(loaded%reaction_source%fields(p)%values==-real(p,dp)), &
              'held restart retains owner-local reaction estimates')
            if(mode==3) call require(all(loaded%heat_release%fields(p)%values==(-1._dp)**p*real(p,dp)*1.e9_dp), &
              'held restart retains signed owner-local heat')
          enddo
          ! An idle rank's expected reset policy remains part of the restart contract.
          kept_gas=gas;kept_temp=temp;kept_record=loaded
          if(rank==nranks-1) loaded%energy_policy%eta2=3
          open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old')
          call read_mpi_eb_tree_checkpoint_3d(trim(mpi_path),stream,species,reactions,transport,cfg,topology,owner, &
            geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
          close(stream);call require(.not.good,'idle or owner reset-policy mismatch refuses collectively')
          do p=1,count
            if(owner(p)/=rank) cycle
            call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
              all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
              all(loaded%internal_energy%fields(p)%values==kept_record%internal_energy%fields(p)%values), &
              'failed held policy read preserves all prior owner fields')
          enddo
          call require(loaded%energy_correction==kept_record%energy_correction.and.loaded%time==kept_record%time, &
            'failed held policy read retains composite metadata')
          loaded=kept_record
          deallocate(gas%fields,temp%fields)
        enddo
        serial_tree=seed;serial_record=eb_tree_checkpoint_3d()
        serial_record%internal_energy=seed_h;serial_record%energy_policy=written%energy_policy
        open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old')
        call reactive_eb_tree_checkpoint_record_3d(stream,.true.,species,reactions,transport,cfg, &
          serial_tree,serial_record,good,message,frame=rotation)
        close(stream);if(.not.good) print *,trim(message)
        call require(good.and.serial_record%energy_correction==written%energy_correction,'serial reads indexed MPI held format')
        do p=1,count
          call require(all(serial_tree%nodes(p)%state==seed%nodes(p)%state).and. &
            all(serial_tree%nodes(p)%temperature==seed%nodes(p)%temperature).and. &
            all(serial_record%internal_energy%fields(p)%values==seed_h%fields(p)%values),'exact MPI/serial held interchange')
        enddo
        if(rank==0) then
          open(newunit=stream,file=trim(serial_path),form='formatted',status='replace')
          call reactive_eb_tree_checkpoint_record_3d(stream,.false.,species,reactions,transport,cfg, &
            serial_tree,serial_record,good,message,frame=rotation)
          close(stream);if(.not.good) print *,trim(message)
          call require(good,'serial held checkpoint writes unindexed envelope')
        endif
        call MPI_Barrier(MPI_COMM_WORLD,status)
        owner=0;allocate(gas%fields(count),temp%fields(count))
        loaded=eb_tree_checkpoint_3d();loaded%energy_policy=written%energy_policy
        allocate(loaded%internal_energy%fields(count))
        do p=1,count
          if(rank/=0) cycle
          gas%fields(p)%values=seed%nodes(p)%state
          allocate(temp%fields(p)%values(1,4,4,4));temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          loaded%internal_energy%fields(p)=seed_h%fields(p)
        enddo
        open(newunit=stream,file=trim(serial_path),access='stream',form='formatted',status='old')
        call read_mpi_eb_tree_checkpoint_3d(trim(serial_path),stream,species,reactions,transport,cfg,topology,owner, &
          geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
        close(stream);call require(good.and.loaded%energy_correction==written%energy_correction,'MPI reads serial held format')
        do p=1,count
          if(rank/=0) cycle
          call require(all(gas%fields(p)%values==seed%nodes(p)%state).and. &
            all(temp%fields(p)%values(1,:,:,:)==seed%nodes(p)%temperature).and. &
            all(loaded%internal_energy%fields(p)%values==seed_h%fields(p)%values),'exact unindexed serial/MPI held fields')
        enddo
        ! Corrupt the last node's finite scalar at the same byte width. Earlier
        ! patch records and metadata parse before the native EOS rejects it.
        if(rank==0) then
          open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old',action='readwrite')
          final_scalar=0;final_metadata=0
          do
            read(stream,'(a)',iostat=code) text
            if(code/=0) exit
            if(text/='HELD_INTERNAL_ENERGY_V1') cycle
            inquire(stream,pos=final_metadata)
            read(stream,'(a)') text
            do scalar_row=1,64
              inquire(stream,pos=final_scalar);read(stream,'(a)') text
            enddo
          enddo
          call require(final_scalar>0.and.final_metadata>0,'locate actual last held patch payload')
          write(stream,'(es27.18e3)',pos=final_scalar) 1.e30_dp
          close(stream)
        endif
        call MPI_Barrier(MPI_COMM_WORLD,status)
        kept_gas=gas;kept_temp=temp;kept_record=loaded
        open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old')
        call read_mpi_eb_tree_checkpoint_3d(trim(mpi_path),stream,species,reactions,transport,cfg,topology,owner, &
          geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
        close(stream);call require(.not.good,'late finite owner UEINT EOS refuses entire MPI checkpoint')
        do p=1,count
          if(rank/=0) cycle
          call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
            all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
            all(loaded%internal_energy%fields(p)%values==kept_record%internal_energy%fields(p)%values), &
            'late held read discards all earlier private patch restores')
        enddo
        call require(loaded%energy_correction==kept_record%energy_correction.and.loaded%time==kept_record%time, &
          'late held checkpoint failure preserves global prefix/time')
        if(rank==0) then
          open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old',action='readwrite')
          write(stream,'(es27.18e3)',pos=final_scalar) seed_h%fields(count)%values(1,4,4,4)
          write(stream,'(3(i0,1x),2(es27.18e3,1x))',pos=final_metadata) &
            merge(1,0,written%energy_policy%allow_small_energy),merge(1,0,written%energy_policy%allow_negative_energy), &
            merge(1,0,written%energy_policy%update_total_from_internal),written%energy_policy%eta2,2*written%energy_correction
          close(stream)
        endif
        call MPI_Barrier(MPI_COMM_WORLD,status)
        open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old')
        call read_mpi_eb_tree_checkpoint_3d(trim(mpi_path),stream,species,reactions,transport,cfg,topology,owner, &
          geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
        close(stream);call require(.not.good,'embedded patch prefix cannot differ from the single composite correction')
        do p=1,count
          if(rank/=0) cycle
          call require(all(gas%fields(p)%values==kept_gas%fields(p)%values).and. &
            all(temp%fields(p)%values==kept_temp%fields(p)%values).and. &
            all(loaded%internal_energy%fields(p)%values==kept_record%internal_energy%fields(p)%values), &
            'embedded prefix mismatch discards every private patch restore')
        enddo
        call require(loaded%energy_correction==kept_record%energy_correction.and.loaded%time==kept_record%time, &
          'embedded prefix mismatch preserves prior global metadata')
        if(nranks>1) then
          if(rank==nranks-1) deallocate(loaded%internal_energy%fields)
          open(newunit=stream,file=trim(mpi_path),access='stream',form='formatted',status='old')
          call read_mpi_eb_tree_checkpoint_3d(trim(mpi_path),stream,species,reactions,transport,cfg,topology,owner, &
            geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
          close(stream);call require(.not.good,'rank-disagreeing held restart presence refuses before payload reads')
          call require(loaded%energy_correction==kept_record%energy_correction,'presence refusal preserves numeric prefix')
          loaded=kept_record
          if(rank==nranks-1) loaded%energy_policy%eta2=3
          open(newunit=stream,status='scratch',form='formatted')
          call write_mpi_eb_tree_checkpoint_3d(trim(serial_path),stream,species,reactions,transport,cfg,topology,owner, &
            geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
          close(stream);call require(.not.good,'idle-rank reset policy mismatch refuses checkpoint publication')
          loaded=kept_record
          open(newunit=stream,file=trim(serial_path),access='stream',form='formatted',status='old')
          call read_mpi_eb_tree_checkpoint_3d(trim(serial_path),stream,species,reactions,transport,cfg,topology,owner, &
            geometry,gas,temp,loaded,MPI_COMM_WORLD,good,frame=rotation)
          close(stream);call require(good.and.loaded%energy_correction==written%energy_correction, &
            'refused collective writer leaves the prior accepted checkpoint readable')
        endif
        deallocate(gas%fields,temp%fields)
      enddo
      deallocate(seed_h%fields)
    enddo
    if(rank==0) print *, 'Held checkpoint phase/ranks: ',trim(checkpoint_action),nranks
  end subroutine check_held_checkpoint
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
    do checkpoint_mode=1,4
      config%chemistry_enabled=mod(checkpoint_mode,2)==1
      owners=[(mod(p-1,nranks),p=1,count)]
      call local_fields()
      if(config%chemistry_enabled) then
        allocate(record%reaction_source%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          allocate(record%reaction_source%fields(p)%values(nv-4,4,4,4),source=-real(p,dp))
        end do
      else
        deallocate(record%reaction_source%fields)
      end if
      if(allocated(record%heat_release%fields)) deallocate(record%heat_release%fields)
      if(checkpoint_mode>=3) then
        allocate(record%heat_release%fields(count))
        do p=1,count
          if(owners(p)/=rank) cycle
          allocate(record%heat_release%fields(p)%values(1,4,4,4),source=(-1._dp)**p*real(p,dp)*1.e9_dp)
        enddo
      endif
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
        if(layout==2) then
          restored%time=17
          if(rank==nranks-1) config%state_redist_max_order=2
          open(newunit=stream,file=path,access='stream',form='formatted',status='old')
          call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
            geometry,local,lt,restored,MPI_COMM_WORLD,good)
          close(stream)
          call require(.not.good.and.restored%time==17,'idle-rank reconstruction mismatch rejects restart')
          call check_fields()
          if(.not.config%chemistry_enabled) then
            open(newunit=stream,status='scratch',form='formatted')
            call write_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
              geometry,local,lt,record,MPI_COMM_WORLD,good)
            close(stream)
            call require(.not.good,'idle-rank reconstruction mismatch rejects checkpoint write')
          endif
          config%state_redist_max_order=0
        endif
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
        if(config%chemistry_enabled) call require(all(serial_record%reaction_source%fields(p)%values==-real(p,dp)), &
          'serial IR exact')
        if(checkpoint_mode>=3) call require(all(serial_record%heat_release%fields(p)%values== &
          (-1._dp)**p*real(p,dp)*1.e9_dp),'serial signed heat exact')
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
        if(config%chemistry_enabled) call require(all(restored%reaction_source%fields(p)%values==-real(p,dp)),'one-rank IR exact')
        if(checkpoint_mode>=3) call require(all(restored%heat_release%fields(p)%values== &
          (-1._dp)**p*real(p,dp)*1.e9_dp),'changed-rank heat exact')
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
      config%state_redist_max_order=2
      open(newunit=stream,file=path,access='stream',form='formatted',status='old')
      call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
        geometry,local,lt,restored,MPI_COMM_WORLD,good)
      close(stream)
      config%state_redist_max_order=0
      call require(.not.good.and.restored%time==17,'changed restart reconstruction order rejects without publication')
      call check_fields()
      config%external_forcing(1)=2
      open(newunit=stream,file=path,access='stream',form='formatted',status='old')
      call read_mpi_eb_tree_checkpoint_3d(path,stream,species,reactions,transport,config,topology,owners, &
        geometry,local,lt,restored,MPI_COMM_WORLD,good)
      close(stream)
      config%external_forcing(1)=1
      call require(.not.good.and.restored%time==17,'mismatched force rejects without publication')
      call check_fields()
      if(config%chemistry_enabled) then
        if(rank==0) then
          open(newunit=stream,file=path,access='stream',form='formatted',status='old',action='readwrite')
          source_count=0
          do
            read(stream,'(a)',iostat=io_code) line
            call require(io_code==0,'locate final source record')
            if(checkpoint_mode<3.and.trim(line)/='SDC_REACTION_SOURCE') cycle
            if(checkpoint_mode>=3.and.trim(line)/='SDC_HEAT_RELEASE') cycle
            source_count=source_count+1
            if(source_count==count) exit
          end do
          do line_index=1,63
            read(stream,'(a)') line
          end do
          inquire(stream,pos=row_position)
          allocate(bad_row(merge(1,nv-4,checkpoint_mode>=3)),source=ieee_value(0._dp,ieee_quiet_nan))
          if(checkpoint_mode>=3) then
            write(stream,'(es27.18e3)',pos=row_position) bad_row
          else
            write(stream,'(*(es27.18e3,1x))',pos=row_position) bad_row
          endif
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
          if(checkpoint_mode>=3) call require(all(restored%heat_release%fields(p)%values== &
            (-1._dp)**p*real(p,dp)*1.e9_dp),'late heat failure rollback')
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
      call require(allocated(restored%heat_release%fields).eqv.(checkpoint_mode>=3),'MPI heat presence')
      do q=1,count
        call require(allocated(local%fields(q)%values).eqv.(owners(q)==rank),'owner-local restart gas')
        if(config%chemistry_enabled) call require(allocated(restored%reaction_source%fields(q)%values).eqv.(owners(q)==rank), &
          'owner-local restart IR')
        if(checkpoint_mode>=3) call require(allocated(restored%heat_release%fields(q)%values).eqv.(owners(q)==rank), &
          'owner-local restart heat')
        if(owners(q)/=rank) cycle
        if(checkpoint_mode>=3) call require(all(restored%heat_release%fields(q)%values== &
          (-1._dp)**q*real(q,dp)*1.e9_dp),'MPI signed heat exact')
        call require(all(local%fields(q)%values==saved%nodes(q)%state),'MPI gas exact')
        call require(all(lt%fields(q)%values(1,:,:,:)==saved%nodes(q)%temperature),'MPI T exact')
        if(config%chemistry_enabled) call require(all(restored%reaction_source%fields(q)%values==-real(q,dp)),'MPI IR exact')
      end do
      if(.not.config%chemistry_enabled) call require(.not.allocated(restored%reaction_source%fields),'chemistry off no IR')
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
