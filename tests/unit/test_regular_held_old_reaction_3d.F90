program test_regular_held_old_reaction_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use mixture_thermo_mod, only: mixture_mass_properties
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_sdc_subcycle_3d
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_sdc_subcycle_3d
  use dual_energy_mod, only: dual_energy_policy
  implicit none
  type :: trial
    type(amr_tree_state_3d) :: gas,temp,held,ir,heat
    real(dp), allocatable :: boundary(:),source(:)
    real(dp) :: prefix,correction,theta,work_change
  end type
  type(trial) :: serial,distributed,saved
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: boxes(1),empty(0)
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction) :: reactions(0)
  type(dual_energy_policy) :: policy
  integer, allocatable :: owners(:)
  real(dp), allocatable :: y(:)
  real(dp) :: spacing(3),dt,mw,rg,cp,cv,gamma,h,e700,e1000,entropy,energy_gap
  integer :: rank,nrank,status,nv,kind,node,n(3),c,calls,initial_calls,cycle
  character(len=16) :: faces(6)
  logical :: good,reject_reactor
  call MPI_Init(status);call MPI_Comm_rank(MPI_COMM_WORLD,rank,status);call MPI_Comm_size(MPI_COMM_WORLD,nrank,status)
  call require(nrank==3,'two owners and idle')
  call load_h2o2_full_thermo(species,good);call require(good,'thermo')
  call load_h2o2_full_transport(transport,good);call require(good,'transport')
  nv=5+size(species);allocate(y(size(species)));y=0
  do c=1,size(species)
    if(species(c)%name=='N2') y(c)=1
  enddo
  call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,e700,entropy,good)
  call require(good,'actualH thermo')
  call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,e1000,entropy,good)
  call require(good,'distinct bulk thermo')
  spacing=.001_dp;dt=2e-8_dp;faces='periodic'
  policy%eta2=2;policy%update_total_from_internal=.false.;policy%allow_negative_energy=.false.
  do kind=1,2
    if(kind==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer ::],empty,tree,good);owners=[0]
    else
      call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),good);call require(good,'fine patch')
      call initialize_amr_patch_tree_3d([4,4,4],[1],boxes,tree,good);owners=[0,1]
    endif
    call require(good,'tree')
    do cycle=1,2
      call initialize(serial,.false.);call initialize(distributed,.true.)
      calls=0;reject_reactor=.false.;energy_gap=0
      call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,serial%gas,serial%temp, &
        spacing,dt,faces,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
        serial%boundary,serial%source,serial%correction,serial%theta,good,reaction_source=serial%ir, &
        evaluate_reaction=reactor,mol_mode=.true.,internal_energy=serial%held,energy_correction=serial%prefix, &
        energy_policy=policy,heat_release=serial%heat,subcycling=cycle==2,internal_energy_pressure_work_change=serial%work_change)
      call require(good.and.calls>0,'serial native actual oldH reactor packing')
      calls=0;energy_gap=0
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,distributed%gas,distributed%temp, &
        spacing,dt,faces,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
        distributed%boundary,distributed%source,distributed%correction,distributed%theta,MPI_COMM_WORLD,good, &
        reaction_source=distributed%ir,evaluate_reaction=reactor,mol_mode=.true.,internal_energy=distributed%held, &
        energy_correction=distributed%prefix,energy_policy=policy,heat_release=distributed%heat,subcycling=cycle==2, &
        internal_energy_pressure_work_change=distributed%work_change)
      call require(good,'MPI native actual oldH reactor packing')
      call require(serial%work_change==0.and.distributed%work_change==0,'native MOL has zero pressure-work source')
      call require((rank>=size(owners)).eqv.(calls==0),'idle rank does not call reactor')
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        call require(all(serial%gas%fields(node)%values==distributed%gas%fields(node)%values).and. &
          all(serial%temp%fields(node)%values==distributed%temp%fields(node)%values).and. &
          all(serial%held%fields(node)%values==distributed%held%fields(node)%values), 'serial MPI exact gas/T/held')
        call require(all(serial%ir%fields(node)%values==distributed%ir%fields(node)%values).and. &
          all(serial%heat%fields(node)%values==distributed%heat%fields(node)%values),'serial MPI exact IR/heat')
      enddo
      call require(serial%prefix==distributed%prefix,'serial MPI reset history')
      saved=distributed;calls=0;reject_reactor=.true.
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,distributed%gas,distributed%temp, &
        spacing,dt,faces,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
        distributed%boundary,distributed%source,distributed%correction,distributed%theta,MPI_COMM_WORLD,good, &
        reaction_source=distributed%ir,evaluate_reaction=reactor,mol_mode=.true.,internal_energy=distributed%held, &
        energy_correction=distributed%prefix,energy_policy=policy,heat_release=distributed%heat,subcycling=cycle==2, &
        internal_energy_pressure_work_change=distributed%work_change)
      call require(.not.good,'reactor failure refuses collectively')
      do node=1,size(owners)
        if(owners(node)/=rank) cycle
        call require(all(saved%gas%fields(node)%values==distributed%gas%fields(node)%values).and. &
          all(saved%temp%fields(node)%values==distributed%temp%fields(node)%values).and. &
          all(saved%held%fields(node)%values==distributed%held%fields(node)%values),'gas/T/H rollback')
        call require(all(saved%ir%fields(node)%values==distributed%ir%fields(node)%values).and. &
          all(saved%heat%fields(node)%values==distributed%heat%fields(node)%values),'IR/heat rollback')
      enddo
      call require(saved%work_change==distributed%work_change,'late reactor refusal preserves pressure-work output')
      call require(saved%prefix==distributed%prefix.and.all(saved%boundary==distributed%boundary).and. &
        all(saved%source==distributed%source),'all history rollback')
    enddo
  enddo
  if(rank==0) print *, 'PASS: regular native serial/MPI actual oldH reactor packing and rollback'
  call MPI_Finalize(status)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(condition) return
    print *, 'FAIL rank',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
  subroutine initialize(t,local)
    type(trial), intent(inout) :: t
    logical, intent(in) :: local
    integer :: peer,dims(3)
    if(allocated(t%gas%fields)) deallocate(t%gas%fields,t%temp%fields,t%held%fields,t%ir%fields,t%heat%fields, &
      t%boundary,t%source)
    allocate(t%gas%fields(size(owners)),t%temp%fields(size(owners)),t%held%fields(size(owners)), &
      t%ir%fields(size(owners)),t%heat%fields(size(owners)),t%boundary(nv),t%source(nv))
    t%boundary=3;t%source=4;t%prefix=5;t%correction=6;t%theta=7;t%work_change=19
    do peer=1,size(owners)
      if(local.and.owners(peer)/=rank) cycle
      dims=tree%nodes(peer)%shape
      allocate(t%gas%fields(peer)%values(nv,dims(1),dims(2),dims(3)),t%temp%fields(peer)%values(1,dims(1),dims(2),dims(3)), &
        t%held%fields(peer)%values(1,dims(1),dims(2),dims(3)),t%ir%fields(peer)%values(nv-4,dims(1),dims(2),dims(3)), &
        t%heat%fields(peer)%values(1,dims(1),dims(2),dims(3)))
      t%gas%fields(peer)%values=0;t%gas%fields(peer)%values(1,:,:,:)=.4_dp
      t%gas%fields(peer)%values(5,:,:,:)=.4_dp*e1000
      do c=1,size(species)
        t%gas%fields(peer)%values(5+c,:,:,:)=.4_dp*y(c)
      enddo
      t%temp%fields(peer)%values=700;t%held%fields(peer)%values=.4_dp*e700
      t%ir%fields(peer)%values=0;t%heat%fields(peer)%values=53
    enddo
  end subroutine
  subroutine reactor(sp,rx,initial_energy,external_y,external_e,duration,rtol,atol,densities,temp,count,accepted,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: initial_energy,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: accepted
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: rho,expected,molecular,gas_r,cpm,cvm,gm,hm,em,sm
    logical :: valid
    calls=calls+1;rho=sum(densities)
    call mixture_mass_properties(sp,densities/rho,temp,molecular,gas_r,cpm,cvm,gm,hm,em,sm,valid)
    call require(valid,'old EOS in reactor observer')
    expected=rho*em
    energy_gap=abs(initial_energy-expected)
    call require(energy_gap<=128*epsilon(1._dp)*max(1._dp,abs(expected)), 'reactor receives actual old H rather than bulk E')
    if(reject_reactor.and.calls==2) then
      densities=-1;temp=-1;count=0;accepted=.false.;return
    endif
    call advance_forced_reactor_interval(sp,rx,initial_energy,external_y,external_e,duration,rtol,atol,densities,temp, &
      count,accepted,maximum_steps=maximum_steps)
  end subroutine
end program
