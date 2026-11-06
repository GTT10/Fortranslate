program test_mpi_eb_tree_remap_3d
  use mpi_f08
  use precision_mod, only: dp
  use dual_energy_mod, only: initialize_nasa7_internal_energy_patch
  use amr_patch_tree_remap_3d_mod, only: covered
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,remap_reactive_eb_tree_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_remap_3d_mod, only: remap_mpi_eb_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(rotating_frame_3d), allocatable :: frame
  type(rotating_frame_3d) :: mismatched_frame
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology,moved,root_only,current_tree
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  type(amr_tree_state_3d) :: local,lt,original,original_t,source,serial_source,original_source
  type(amr_tree_state_3d) :: heat,serial_heat,original_heat
  real(dp) :: heat_before(1),heat_after(1)
  type(eb_geometry_3d), allocatable :: geometry(:),original_geometry(:)
  integer, allocatable :: owners(:),next_owners(:)
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),source_before(:),source_after(:)
  real(dp) :: residual,position(3),potential
  real(dp), allocatable :: primitive(:)
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k,rank,nranks,status,stage,rotation
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),primitive(reactive_nprim(size(species))))
  do rotation=0,1
  if(rotation==1) then
    allocate(frame)
    frame%angular_speed=100._dp;frame%axis_point=[10._dp,10._dp,10._dp]
  end if
  do axis=1,3
  if(allocated(frame)) frame%axis=axis
  do count=4,5
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
  if(allocated(frame)) then
    do k=1,4
    do j=1,4
    do i=1,4
      if(coarse%volume_fraction(i,j,k)<=0) cycle
      position=(real([i,j,k],dp)-.5_dp)/4
      call rotating_potential_3d(frame,position,potential,ok);call require(ok,'root potential')
      u(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
    end do
    end do
    end do
  end if
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok,frame)
  call require(ok,'tree initialization')
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  do i=1,4
    tree%nodes(4)%state(:,i,:,:)=tree%nodes(4)%state(:,i,:,:)*(1._dp+.01_dp*real(i,dp))
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok,frame)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree;current_tree=topology;stage=0
  allocate(geometry(count),owners(count),local%fields(count),lt%fields(count))
  owners=[(mod(node-1,nranks),node=1,count)]
  do node=1,count
    geometry(node)=tree%nodes(node)%geometry
    if (owners(node)/=rank) cycle
    local%fields(node)%values=tree%nodes(node)%state
    allocate(lt%fields(node)%values(1,4,4,4))
    lt%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
  end do
  allocate(source%fields(count),serial_source%fields(count),source_before(nv-4),source_after(nv-4))
  do node=1,count
    allocate(serial_source%fields(node)%values(nv-4,4,4,4),source=-real(node,dp))
    do i=1,4
      serial_source%fields(node)%values(:,i,:,:)=serial_source%fields(node)%values(:,i,:,:)+.1_dp*i
    end do
    if(owners(node)==rank) source%fields(node)%values=serial_source%fields(node)%values
  end do
  allocate(heat%fields(count),serial_heat%fields(count))
  do node=1,count
    allocate(serial_heat%fields(node)%values(1,4,4,4))
    serial_heat%fields(node)%values=serial_source%fields(node)%values(:1,:,:,:)
    if(owners(node)==rank) heat%fields(node)%values=serial_heat%fields(node)%values
  end do
  call integrate_mpi_eb_tree_3d(current_tree,owners,geometry,heat,heat_before,MPI_COMM_WORLD,ok)
  call require(ok,'initial signed heat inventory')
  call integrate_mpi_eb_tree_3d(current_tree,owners,geometry,source,source_before,MPI_COMM_WORLD,ok)
  call require(ok,'initial source inventory')
  call initialize_amr_patch_3d(4,4,4,2,3,3,4,2,3,2,boxes(3),ok)
  call require(ok,'shifted grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),moved,ok)
  call require(ok,'new topology')
  call remap_reactive_eb_tree_3d(species,moved,finest,tree,ok,serial_source,frame,serial_heat)
  call require(ok,'move grandchild and remove other refinement')
  call check_inventory()
  call check_mpi_remap(tree%topology)
  call require(all(tree%nodes(4)%state(:,3:4,:,:)==saved%nodes(4)%state(:,1:2,:,:)), &
    'exact same-level overlap states')
  call require(all(tree%nodes(4)%temperature(3:4,:,:)==saved%nodes(4)%temperature(1:2,:,:)), &
    'exact same-level overlap temperatures')
  ! Root-only target coarsens all old descendants; the input finest geometry
  ! still describes the deeper of the source and target trees.
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],boxes(:0),root_only,ok)
  call require(ok,'root-only topology')
  call remap_reactive_eb_tree_3d(species,root_only,finest,tree,ok,serial_source,frame,serial_heat)
  call require(ok.and.size(tree%nodes)==1,'complete coarsening')
  call check_inventory()
  call check_mpi_remap(tree%topology)
  call remap_reactive_eb_tree_3d(species,topology,finest,tree,ok,serial_source,frame,serial_heat)
  call require(ok,'refine root-only state to three levels')
  call check_inventory()
  call check_mpi_remap(tree%topology)
  original=local;original_t=lt;original_geometry=geometry;original_source=source;original_heat=heat
  if(allocated(frame).and.nranks>1) then
    mismatched_frame=frame
    if(rank==1) mismatched_frame%angular_speed=mismatched_frame%angular_speed+1._dp
    call average_down_mpi_eb_tree_3d(species,current_tree,owners,geometry,local,lt,MPI_COMM_WORLD,ok,mismatched_frame)
    call require(.not.ok,'rank-mismatched frame restriction rejected')
    call unchanged_mpi()
    next_owners=[(mod(node,nranks),node=1,size(moved%nodes))]
    call remap_mpi_eb_tree_3d(species,current_tree,moved,owners,next_owners,finest,geometry,local,lt, &
      MPI_COMM_WORLD,ok,source,mismatched_frame)
    call require(.not.ok,'rank-mismatched frame rejected')
    call unchanged_mpi()
  end if
  if (owners(count)==rank) local%fields(count)%values(1,4,4,4)=-1._dp
  original=local
  next_owners=[(mod(node,nranks),node=1,size(moved%nodes))]
  call remap_mpi_eb_tree_3d(species,current_tree,moved,owners,next_owners,finest,geometry,local,lt, &
    MPI_COMM_WORLD,ok,source,frame,heat)
  call require(.not.ok,'MPI bad source rejected')
  call unchanged_mpi()
  if (owners(count)==rank) local%fields(count)%values(1,4,4,4)=tree%nodes(count)%state(1,4,4,4)
  original=local
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    axes(axis),.32_dp,finest,ok)
  call require(ok,'changed geometry fixture')
  call remap_mpi_eb_tree_3d(species,current_tree,moved,owners,next_owners,finest,geometry,local,lt, &
    MPI_COMM_WORLD,ok,source,frame,heat)
  call require(.not.ok,'MPI geometry change rejected')
  call unchanged_mpi()
  deallocate(geometry,owners,local%fields,lt%fields,source%fields,serial_source%fields,source_before,source_after)
  deallocate(heat%fields,serial_heat%fields)
  deallocate(u,t)
  end do
  end do
  end do
  call check_held_remap()
  if (rank==0) print *, 'test_mpi_eb_tree_remap_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_held_remap()
    type(eb_geometry_3d) :: fg
    type(eb_geometry_3d), allocatable :: geom(:),g_saved(:)
    type(eb_reactive_tree_3d) :: seed,reference
    type(amr_patch_tree_3d) :: cur,target
    type(rotating_frame_3d), allocatable :: spin
    type(amr_tree_state_3d) :: rh,rir,rheat,gas,temp,held,ir,hfield,gas_saved,temp_saved,held_saved,ir_saved,heat_saved
    integer, allocatable :: old_owner(:),new_owner(:)
    real(dp) :: baseline(2*nv-2),inventory(2*nv-2),rho,position(3),phi,recovered,cs,q(nv)
    integer :: p,x,y,z,condition,step,nodes
    logical :: good
    call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,fg,good)
    call require(good,'native held remap immutable geometry')
    do condition=1,2
      if(condition==2) then
        allocate(spin);spin%axis=2;spin%angular_speed=1.e10_dp;spin%axis_point=[10._dp,10._dp,10._dp]
      endif
      seed=tree;nodes=size(seed%nodes)
      allocate(rh%fields(nodes),rir%fields(nodes),rheat%fields(nodes))
      do p=1,nodes
        allocate(rh%fields(p)%values(1,4,4,4),source=0._dp)
        allocate(rir%fields(p)%values(nv-4,4,4,4),source=-real(p,dp))
        allocate(rheat%fields(p)%values(1,4,4,4),source=real(p,dp)*.5_dp)
        call initialize_nasa7_internal_energy_patch(species,seed%nodes(p)%state,seed%nodes(p)%temperature, &
          seed%nodes(p)%geometry%volume_fraction,rh%fields(p)%values(1,:,:,:),good)
        call require(good,'held remap native initial rho*e independent of total subtraction')
        do z=1,4
        do y=1,4
        do x=1,4
          if(seed%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          rho=seed%nodes(p)%state(1,x,y,z)
          seed%nodes(p)%state(2:4,x,y,z)=[rho*1.e10_dp,0._dp,0._dp]
          seed%nodes(p)%state(5,x,y,z)=.5_dp*rho*1.e20_dp
          if(allocated(spin)) then
            associate(g=>seed%nodes(p)%geometry)
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_potential_3d(spin,position,phi,good)
            call require(good,'held remap native modified energy convention')
            seed%nodes(p)%state(5,x,y,z)=seed%nodes(p)%state(5,x,y,z)-rho*phi
          endif
        enddo
        enddo
        enddo
      enddo
      reference=seed;cur=seed%topology
      allocate(old_owner(nodes),geom(nodes),gas%fields(nodes),temp%fields(nodes),held%fields(nodes), &
        ir%fields(nodes),hfield%fields(nodes))
      old_owner=[(mod(p-1,nranks),p=1,nodes)]
      do p=1,nodes
        geom(p)=seed%nodes(p)%geometry
        if(old_owner(p)/=rank) cycle
        gas%fields(p)%values=seed%nodes(p)%state
        allocate(temp%fields(p)%values(1,4,4,4))
        temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
        held%fields(p)=rh%fields(p);ir%fields(p)=rir%fields(p);hfield%fields(p)=rheat%fields(p)
      enddo
      call held_remap_inventory(reference,rh,rir,rheat,baseline)
      do step=1,3
        if(step==1) target=moved
        if(step==2) target=root_only
        if(step==3) target=topology
        new_owner=[(mod(p+step,nranks),p=1,size(target%nodes))]
        call remap_reactive_eb_tree_3d(species,target,fg,reference,good,rir,spin,rheat,rh)
        call require(good,'native serial held remap moves/coarsens/refines without total-derived root recovery')
        call held_remap_inventory(reference,rh,rir,rheat,inventory)
        call require(maxval(abs(inventory-baseline)/max(1._dp,abs(baseline)))<3.e-14_dp, &
          'independent gas/UEINT/IR/heat composite inventories survive native held remap')
        call remap_mpi_eb_tree_3d(species,cur,target,old_owner,new_owner,fg,geom,gas,temp, &
          MPI_COMM_WORLD,good,ir,spin,hfield,held)
        call require(good,'native MPI held remap changes owners through move/coarsening/refinement')
        do p=1,size(target%nodes)
          call require(allocated(held%fields(p)%values).eqv.(new_owner(p)==rank), &
            'remapped UEINT has only its new owner-local payload')
          call require(matching_eb_metrics_3d(geom(p),reference%nodes(p)%geometry),'paired remap geometry matches serial')
          if(new_owner(p)/=rank) cycle
          call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
            all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
            all(held%fields(p)%values==rh%fields(p)%values).and. &
            all(ir%fields(p)%values==rir%fields(p)%values).and.all(hfield%fields(p)%values==rheat%fields(p)%values), &
            'native serial/MPI remapped gas/T/UEINT/IR/heat agree exactly')
          do z=1,4
          do y=1,4
          do x=1,4
            if(geom(p)%volume_fraction(x,y,z)<=0) cycle
            call recover_reactive_eb_cell_3d(species,geom(p),[x,y,z],gas%fields(p)%values(:,x,y,z), &
              temp%fields(p)%values(1,x,y,z),q,recovered,cs,good,spin,held%fields(p)%values(1,x,y,z))
            call require(good.and.abs(recovered/temp%fields(p)%values(1,x,y,z)-1)<1.e-12_dp, &
              'remapped native temperature is consistent with directly retained UEINT')
          enddo
          enddo
          enddo
        enddo
        cur=target;old_owner=new_owner
      enddo
      if(old_owner(2)==rank) held%fields(2)%values(1,4,4,4)=1.e30_dp
      gas_saved=gas;temp_saved=temp;held_saved=held;ir_saved=ir;heat_saved=hfield;g_saved=geom
      call remap_mpi_eb_tree_3d(species,cur,moved,old_owner, &
        [(mod(p,nranks),p=1,size(moved%nodes))],fg,geom,gas,temp,MPI_COMM_WORLD,good,ir,spin,hfield,held)
      call require(.not.good,'native late held EOS failure refuses remap collectively after private source staging')
      do p=1,size(cur%nodes)
        call require(matching_eb_metrics_3d(geom(p),g_saved(p)),'failed paired remap preserves all geometry')
        if(old_owner(p)/=rank) cycle
        call require(all(gas%fields(p)%values==gas_saved%fields(p)%values).and. &
          all(temp%fields(p)%values==temp_saved%fields(p)%values).and. &
          all(held%fields(p)%values==held_saved%fields(p)%values).and. &
          all(ir%fields(p)%values==ir_saved%fields(p)%values).and. &
          all(hfield%fields(p)%values==heat_saved%fields(p)%values), &
          'late EOS failed remap preserves all prior gas/T/UEINT/IR/heat payloads')
      enddo
      if(nranks>1) then
        if(rank==nranks-1) then
          call remap_mpi_eb_tree_3d(species,cur,moved,old_owner, &
            [(mod(p,nranks),p=1,size(moved%nodes))],fg,geom,gas,temp,MPI_COMM_WORLD,good,ir,spin,hfield)
        else
          call remap_mpi_eb_tree_3d(species,cur,moved,old_owner, &
            [(mod(p,nranks),p=1,size(moved%nodes))],fg,geom,gas,temp,MPI_COMM_WORLD,good,ir,spin,hfield,held)
        endif
        call require(.not.good,'rank-disagreeing remap UEINT presence refuses collectively')
        do p=1,size(cur%nodes)
          if(old_owner(p)/=rank) cycle
          call require(all(gas%fields(p)%values==gas_saved%fields(p)%values).and. &
            all(temp%fields(p)%values==temp_saved%fields(p)%values).and. &
            all(held%fields(p)%values==held_saved%fields(p)%values).and. &
            all(ir%fields(p)%values==ir_saved%fields(p)%values).and. &
            all(hfield%fields(p)%values==heat_saved%fields(p)%values), &
            'failed paired remap preserves all prior gas/T/UEINT/IR/heat payloads')
        enddo
      endif
      deallocate(rh%fields,rir%fields,rheat%fields,old_owner,geom,gas%fields,temp%fields,held%fields, &
        ir%fields,hfield%fields)
    enddo
  end subroutine check_held_remap

  subroutine held_remap_inventory(layout,stored,ir,heat,inventory)
    type(eb_reactive_tree_3d), intent(in) :: layout
    type(amr_tree_state_3d), intent(in) :: stored,ir,heat
    real(dp), intent(out) :: inventory(:)
    real(dp) :: volume
    integer :: p,x,y,z
    inventory=0
    do p=1,size(layout%nodes)
      associate(g=>layout%nodes(p)%geometry)
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0.or.covered(layout%topology,p,[x,y,z])) cycle
        volume=g%dx*g%dy*g%dz*g%volume_fraction(x,y,z)
        inventory=inventory+volume*[layout%nodes(p)%state(:,x,y,z),stored%fields(p)%values(1,x,y,z), &
          ir%fields(p)%values(:,x,y,z),heat%fields(p)%values(1,x,y,z)]
      enddo
      enddo
      enddo
      end associate
    enddo
  end subroutine held_remap_inventory

  subroutine check_inventory()
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok,'remapped inventory')
    residual=maxval(abs(after-before)/max(1._dp,abs(before)))
    print *, 'axis, source nodes, residual: ',axes(axis),count,residual
    call require(residual<3.e-14_dp,'remap composite conservation')
  end subroutine check_inventory
  subroutine check_mpi_remap(target)
    type(amr_patch_tree_3d), intent(in) :: target
    real(dp) :: inventory(nv),error,recovered,speed
    integer :: p
    stage=stage+1
    next_owners=[(mod(p+stage,nranks),p=1,size(target%nodes))]
    call remap_mpi_eb_tree_3d(species,current_tree,target,owners,next_owners,finest,geometry,local,lt, &
      MPI_COMM_WORLD,ok,source,frame,heat)
    call require(ok,'MPI remap')
    call integrate_mpi_eb_tree_3d(target,next_owners,geometry,heat,heat_after,MPI_COMM_WORLD,ok)
    call require(ok.and.all(abs(heat_after-heat_before)<3.e-14_dp*max(1._dp,abs(heat_before))), &
      'MPI signed heat remap conservation')
    current_tree=target;owners=next_owners
    do p=1,size(target%nodes)
      call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'new owner-local allocation')
      call require(allocated(source%fields(p)%values).eqv.(owners(p)==rank),'owner-local remapped source')
      call require(allocated(heat%fields(p)%values).eqv.(owners(p)==rank),'owner-local remapped heat')
      call require(matching_eb_metrics_3d(geometry(p),tree%nodes(p)%geometry),'serial target geometry')
      if (owners(p)/=rank) cycle
      call require(all(source%fields(p)%values==serial_source%fields(p)%values),'exact serial remapped IR')
      call require(all(heat%fields(p)%values==serial_heat%fields(p)%values),'exact serial remapped stage heat')
      call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial remapped U')
      call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial remapped T')
      if(allocated(frame)) then
        do k=1,4
        do j=1,4
        do i=1,4
          if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
          position=[geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower]+ &
            (real([i,j,k],dp)-.5_dp)*[geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]
          call rotating_potential_3d(frame,position,potential,ok);call require(ok,'MPI receiving potential')
          call reactive_conserved_to_primitive(species,local%fields(p)%values(:,i,j,k), &
            lt%fields(p)%values(1,i,j,k),primitive,recovered,speed,ok,potential)
          call require(ok,'MPI receiving EOS')
          call require(abs(recovered/lt%fields(p)%values(1,i,j,k)-1)<1.e-12_dp,'MPI receiving temperature')
        end do
        end do
        end do
      end if
    end do
    call integrate_mpi_eb_tree_3d(target,owners,geometry,local,inventory,MPI_COMM_WORLD,ok)
    call require(ok,'MPI inventory')
    call integrate_mpi_eb_tree_3d(target,owners,geometry,source,source_after,MPI_COMM_WORLD,ok)
    call require(ok,'remapped source inventory')
    call require(maxval(abs(source_after-source_before)/max(1._dp,abs(source_before)))<3.e-14_dp, &
      'signed source conservation')
    error=maxval(abs(inventory-before)/max(1._dp,abs(before)))
    call require(error<3.e-14_dp,'MPI remap conservation')
    if (rank==0) print *, 'MPI axis, nodes, stage, residual: ',axes(axis),count,stage,error
  end subroutine
  subroutine unchanged_mpi()
    integer :: p
    do p=1,size(owners)
      call require(matching_eb_metrics_3d(geometry(p),original_geometry(p)),'geometry rollback')
      if (owners(p)/=rank) cycle
      call require(all(source%fields(p)%values==original_source%fields(p)%values),'IR rollback')
      call require(all(heat%fields(p)%values==original_heat%fields(p)%values),'stage heat rollback')
      call require(all(local%fields(p)%values==original%fields(p)%values),'U rollback')
      call require(all(lt%fields(p)%values==original_t%fields(p)%values),'T rollback')
    end do
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_remap_3d
