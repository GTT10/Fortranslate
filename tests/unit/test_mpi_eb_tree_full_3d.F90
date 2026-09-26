program test_mpi_eb_tree_full_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_full_3d,react_sdc_eb_tree_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_full_3d,react_sdc_mpi_eb_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
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
  call check_sdc_reaction()
  do integrator_index=1,2
  do mode=0,3
    if (integrator_index==2.and.mode/=3) cycle
    chemistry_on=btest(mode,0);transport_on=btest(mode,1)
    tree=saved
    call advance_reactive_eb_tree_full_3d(species,reactions,transport,'hllc',1.e-8_dp,chemistry_on, &
      2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true.,tree,boundary,wall,chemical,theta,ok, &
      integrators(integrator_index))
    call require(ok,'combined step')
    call integrate_reactive_eb_tree_3d(tree,after,ok)
    call require(ok,'final inventory')
    residual=maxval(abs(after-before-boundary-wall-chemical)/max(1._dp,abs(before)))
    print *, 'integrator, mode, residual: ',integrators(integrator_index),mode,residual
    call require(residual<3.e-14_dp,'combined conservation')
    call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
    if (.not.transport_on) call require(theta==1._dp,'disabled transport limiter')
    if (.not.chemistry_on) call require(all(chemical==0._dp),'disabled chemistry ledger')
    if (chemistry_on) then
      call require(any(chemical(6:)/=0._dp),'nonzero chemical species exchange')
      call require(all(abs(chemical(:5))<=3.e-14_dp*max(1._dp,abs(before(:5)))),'chemical invariants')
    end if
    call check_mpi_full()
  end do
  end do
  deallocate(u,t)
  end do
  end do
  if (rank==0) print *, 'test_mpi_eb_tree_full_3d: PASS'
  call MPI_Finalize(status)
contains
  subroutine check_sdc_reaction()
    type(amr_tree_state_3d) :: rates,ref_state,ref_temp,ref_ir,local,lt,lr,ir,old,ot,kept,kt,ki
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: tolerance
    integer :: p,layout,owners(count)
    logical :: good
    allocate(rates%fields(count),ref_ir%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      allocate(rates%fields(p)%values(nv,4,4,4),source=0._dp)
      allocate(ref_ir%fields(p)%values(nv-4,4,4,4),source=7._dp)
      rates%fields(p)%values(2,:,:,:)=2.e4_dp
      rates%fields(p)%values(5,:,:,:)=2.e4_dp*saved%nodes(p)%state(2,:,:,:)/saved%nodes(p)%state(1,:,:,:)
    end do
    call react_sdc_eb_tree_3d(species,reactions,saved,rates,1.e-9_dp,2.e-7_dp,1.e-12_dp, &
      ref_state,ref_temp,ref_ir,good)
    call require(good,'serial SDC reaction reference')
    do layout=1,2
      owners=[(mod(p-1,nranks),p=1,count)]
      if(layout==2) owners=0
      allocate(old%fields(count),ot%fields(count),lr%fields(count),ir%fields(count))
      do p=1,count
        if(owners(p)/=rank) cycle
        old%fields(p)%values=saved%nodes(p)%state
        allocate(ot%fields(p)%values(1,4,4,4))
        ot%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
        lr%fields(p)%values=rates%fields(p)%values
        allocate(ir%fields(p)%values(nv-4,4,4,4),source=7._dp)
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(good,'MPI SDC reaction')
      do p=1,count
        call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'SDC owner-local allocation')
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==ref_state%fields(p)%values).and. &
          all(lt%fields(p)%values==ref_temp%fields(p)%values).and.all(ir%fields(p)%values==ref_ir%fields(p)%values), &
          'MPI SDC exact serial state temperature and reaction')
        local%fields(p)%values=-17;lt%fields(p)%values=-19;ir%fields(p)%values=7
      end do
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good,target_depth=2)
      call require(good,'MPI selected level SDC reaction')
      do p=1,count
        if(owners(p)/=rank) cycle
        if(topology%nodes(p)%level==2) then
          call require(all(local%fields(p)%values==ref_state%fields(p)%values).and. &
            all(ir%fields(p)%values==ref_ir%fields(p)%values),'MPI selected SDC exact reaction')
        else
          call require(all(local%fields(p)%values==-17).and.all(lt%fields(p)%values==-19).and. &
            all(ir%fields(p)%values==7),'MPI SDC inactive level retained')
        end if
      end do
      kept=local;kt=lt;ki=ir
      if(owners(count)==rank) old%fields(count)%values(1,4,4,4)=-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,2.e-7_dp,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(.not.good,'MPI failing owner rejects globally')
      if(owners(count)==rank) old%fields(count)%values=saved%nodes(count)%state
      tolerance=2.e-7_dp
      if(rank==nranks-1) tolerance=-1
      call react_sdc_mpi_eb_tree_3d(species,reactions,topology,owners,geometry,old,ot,lr, &
        1.e-9_dp,tolerance,1.e-12_dp,local,lt,ir,MPI_COMM_WORLD,good)
      call require(.not.good,'MPI invalid control including idle rank rejects')
      do p=1,count
        if(owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==kept%fields(p)%values).and. &
          all(lt%fields(p)%values==kt%fields(p)%values).and.all(ir%fields(p)%values==ki%fields(p)%values), &
          'MPI SDC collective rollback')
      end do
      deallocate(old%fields,ot%fields,lr%fields,ir%fields)
    end do
  end subroutine check_sdc_reaction

  subroutine check_mpi_full()
    type(amr_tree_state_3d) :: local,lt,initial,it
    type(eb_geometry_3d) :: geometry(count)
    real(dp) :: b(nv),w(nv),chem(nv),limiter,inventory(nv),error
    integer :: owners(count),p
    logical :: good,enabled
    owners=[(mod(p-1,nranks),p=1,count)]
    if (integrator_index==2) owners=nranks-1
    allocate(local%fields(count),lt%fields(count))
    do p=1,count
      geometry(p)=saved%nodes(p)%geometry
      if (owners(p)/=rank) cycle
      local%fields(p)%values=saved%nodes(p)%state
      allocate(lt%fields(p)%values(1,4,4,4))
      lt%fields(p)%values(1,:,:,:)=saved%nodes(p)%temperature
    end do
    initial=local;it=lt
    call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'hllc',1.e-8_dp, &
      chemistry_on,2.e-7_dp,1.e-12_dp,transport_on,.true.,.true.,.true.,.true., &
      local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,integrators(integrator_index))
    call require(good,'MPI full step')
    do p=1,count
      call require(allocated(local%fields(p)%values).eqv.(owners(p)==rank),'owner-local full U')
      if (owners(p)/=rank) cycle
      call require(all(local%fields(p)%values==tree%nodes(p)%state),'exact serial full U')
      call require(all(lt%fields(p)%values(1,:,:,:)==tree%nodes(p)%temperature),'exact serial full T')
    end do
    call require(all(b==boundary).and.limiter==theta,'exact boundary and limiter')
    call require(maxval(abs(w-wall)/max(1._dp,abs(before)))<3.e-14_dp,'wall parity')
    call require(maxval(abs(chem-chemical)/max(1._dp,abs(before)))<3.e-14_dp,'chemical ledger parity')
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,inventory,MPI_COMM_WORLD,good)
    call require(good,'MPI inventory')
    error=maxval(abs(inventory-before-b-w-chem)/max(1._dp,abs(before)))
    call require(error<3.e-14_dp,'MPI full conservation')
    if (rank==0) print *, 'MPI integrator, mode, residual: ',integrators(integrator_index),mode,error
    if (integrator_index==1.and.mode==3) then
      local=initial;lt=it;b=17._dp;w=19._dp;chem=23._dp;limiter=29._dp
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'invalid',1.e-8_dp, &
        .true.,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit')
      call require(.not.good,'late hydro rejection after reaction and transport')
      do p=1,count
        if (owners(p)/=rank) cycle
        call require(all(local%fields(p)%values==initial%fields(p)%values),'full U rollback')
        call require(all(lt%fields(p)%values==it%fields(p)%values),'full T rollback')
      end do
      call require(all(b==17._dp).and.all(w==19._dp).and.all(chem==23._dp).and.limiter==29._dp,'full ledger rollback')
      enabled=rank==0
      call advance_mpi_eb_tree_full_3d(species,reactions,transport,topology,owners,geometry,'hllc',1.e-8_dp, &
        enabled,2.e-7_dp,1.e-12_dp,.true.,.true.,.true.,.true.,.true., &
        local,lt,b,w,chem,limiter,MPI_COMM_WORLD,good,'explicit')
      call require(.not.good,'different chemistry flags rejected')
    end if
  end subroutine check_mpi_full
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (.not.condition) then
      print *, 'FAIL: ',label
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine require
end program test_mpi_eb_tree_full_3d
