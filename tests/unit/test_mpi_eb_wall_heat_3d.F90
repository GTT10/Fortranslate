program test_mpi_eb_wall_heat_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_mass_fraction_component,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_heat_source_3d
  use mpi_eb_wall_heat_3d_mod, only: compute_mpi_eb_wall_heat_3d
  use mpi_eb_tree_transport_3d_mod, only: compute_sdc_mpi_eb_tree_transport_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use elementary_kinetics_mod, only: elementary_reaction
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(2)
  type(amr_patch_tree_3d) :: topology
  type(eb_reactive_tree_3d) :: serial
  type(eb_geometry_3d) :: padded,finest,geometry(3),padded_six
  type(eb_geometry_3d), allocatable :: extended(:),extended_six(:)
  type(amr_tree_state_3d) :: local,lt,rate,stage,baseline
  type(amr_tree_state_3d) :: evolved,evolved_t
  type(eb_integrated_flux_3d), allocatable :: flux(:)
  real(dp), allocatable :: q(:),y(:),u(:),root(:,:,:,:),rt(:,:,:),wide(:,:,:,:),wt(:,:,:),heat(:,:,:)
  real(dp) :: mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,t,sound,expected,ref_frequency,power,frequency
  real(dp) :: theta,wall_energy,stage_frequency,control
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),chemical(:),external(:)
  real(dp) :: residual
  character(len=16) :: faces(6)
  integer :: status,rank,ranks,nv,s,node,owners(3),n(3),mode,scheme
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,ranks,status)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(reactions(0),before(nv),after(nv),boundary(nv),wall(nv),chemical(nv),external(nv))
  faces='outflow'
  allocate(q(reactive_nprim(size(species))),y(size(species)),u(nv),root(nv,4,4,4),rt(4,4,4))
  y=0
  do s=1,size(species)
    if(species(s)%name=='N2') y(s)=1
  end do
  call mixture_mass_properties(species,y,500._dp,mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,ok)
  call require(ok,'properties')
  q=0;q(1)=1;q(5)=rgas*500
  do s=1,size(species)
    q(reactive_mass_fraction_component(s))=y(s)
  end do
  call reactive_primitive_to_conserved(species,q,u,t,sound,ok);call require(ok,'state')
  do s=1,nv
    root(s,:,:,:)=u(s)
  end do
  rt=500
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),ok);call require(ok,'left patch')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,boxes(2),ok);call require(ok,'right patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1],boxes,topology,ok);call require(ok,'tree')
  call build_axis_plane_eb_geometry_3d(16,16,16,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp, &
    'x',.43_dp,padded,ok);call require(ok,'padded geometry')
  call extract_eb_geometry_3d(padded,[5,5,5],[12,12,12],finest,ok);call require(ok,'finest')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,rt,serial,ok);call require(ok,'serial tree')
  do node=1,3
    geometry(node)=serial%nodes(node)%geometry
  end do
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,ok,halo_width=2)
  call require(ok,'extended tree geometry')
  call build_axis_plane_eb_geometry_3d(32,32,32,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
    'x',.43_dp,padded_six,ok);call require(ok,'hydro padded geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded_six,extended_six,ok)
  call require(ok,'hydro extended tree geometry')
  allocate(wide(nv,16,16,16),wt(16,16,16),heat(8,8,8))
  do s=1,nv
    wide(s,:,:,:)=u(s)
  end do
  wt=500
  call compute_eb_wall_heat_source_3d(species,transport,wide,wt,padded,[5,5,5],[12,12,12], &
    400._dp,heat,expected,ref_frequency,ok);call require(ok,'whole fine reference')
  do mode=1,2
    owners=[0,0,mod(1,ranks)]
    if(mode==2) owners=0
    allocate(local%fields(3),lt%fields(3))
    do node=1,3
      if(owners(node)/=rank) cycle
      n=topology%nodes(node)%shape
      local%fields(node)%values=serial%nodes(node)%state
      allocate(lt%fields(node)%values(1,n(1),n(2),n(3)))
      lt%fields(node)%values(1,:,:,:)=serial%nodes(node)%temperature
    end do
    call compute_mpi_eb_wall_heat_3d(species,transport,topology,owners,geometry,extended,local,lt,400._dp, &
      rate,power,frequency,MPI_COMM_WORLD,ok)
    call require(ok,'distributed heat')
    call require(abs(power-expected)<1e-12_dp*abs(expected),'composite heat excludes covered coarse wall')
    call require(frequency>0,'frequency')
    call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1e-7_dp, &
      .false.,.true.,.false.,.false.,local,lt,baseline,flux,theta,MPI_COMM_WORLD,ok)
    call require(ok,'baseline MPI stage')
    call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1e-7_dp, &
      .false.,.true.,.false.,.false.,local,lt,stage,flux,theta,MPI_COMM_WORLD,ok, &
      embedded_wall_temperature=400._dp,extended_geometry=extended,embedded_wall_energy=wall_energy, &
      embedded_frequency=stage_frequency)
    call require(ok,'heated MPI stage')
    call require(abs(wall_energy-1e-7_dp*power)<1e-20_dp,'integrated stage ledger')
    do node=1,3
      if(owners(node)/=rank) cycle
      call require(maxval(abs(stage%fields(node)%values-baseline%fields(node)%values-rate%fields(node)%values)) &
        <1e-10_dp,'source routed once')
    end do
    control=400
    if(rank==ranks-1) control=-1
    power=123;frequency=456
    call compute_mpi_eb_wall_heat_3d(species,transport,topology,owners,geometry,extended,local,lt,control, &
      rate,power,frequency,MPI_COMM_WORLD,ok)
    call require(.not.ok.and.power==123.and.frequency==456,'collective failure rollback')
    if(mode==1) then
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,before,MPI_COMM_WORLD,ok)
      call require(ok,'initial inventory')
      do scheme=1,2
        evolved=local;evolved_t=lt
        boundary=0;wall=0;chemical=0;external=0;theta=1
        if(scheme==1) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
            extended_six,1e-7_dp,2,.false.,[.false.,.true.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
            embedded_wall_temperature=400._dp)
        else
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
            1e-7_dp,2,.false.,[.false.,.true.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,extended_six,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
            embedded_wall_temperature=400._dp)
        end if
        call require(ok,'evolved isothermal wall')
        call integrate_mpi_eb_tree_3d(topology,owners,geometry,evolved,after,MPI_COMM_WORLD,ok)
        call require(ok,'final inventory')
        residual=maxval(abs(after-before-boundary-wall-chemical-external)/max(1._dp,abs(before)))
        call require(residual<1e-12_dp,'time-integrated conservation')
        call require(wall(5)<0.and.abs(wall(5)-1e-7_dp*expected)<1e-5_dp*abs(1e-7_dp*expected), &
          'time-integrated cooling ledger')
        if(rank==0) print *, 'Wall scheme, heat, budget residual:',scheme,wall(5),residual
      end do
    end if
    deallocate(local%fields,lt%fields)
  end do
  if(rank==0) print *, 'PASS: cross-owner donors, composite wall heat, SDC source and collective rejection',expected
  call MPI_Finalize(status)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: code
    if(condition) return
    print *, 'FAIL: ',message
    call MPI_Abort(MPI_COMM_WORLD,1,code)
  end subroutine
end program
