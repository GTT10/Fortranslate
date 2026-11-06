program test_amr_tree_initial_reaction
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use lidryer_mechanism_mod, only: load_lidryer_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use mixture_thermo_mod, only: mixture_mass_properties,temperature_from_internal_energy
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_chemistry_3d_mod, only: initialize_amr_tree_reaction_source_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: initialize_mpi_amr_tree_reaction_source_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patch(1)
  type(amr_tree_state_3d) :: gas,temp,held,ir,local_gas,local_temp,local_held,local_ir,saved
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction) :: reactions(0)
  real(dp) :: y(9),mw,r,cp,cv,gamma,h,e,s,dt,expected,sum_dt
  integer :: rank,nranks,status,node,n(3),cycle,calls,owners(2),reject_after
  logical :: good
  call MPI_Init(status);call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  call require(nranks==3,'two owners and one idle rank')
  call load_lidryer_thermo(species,good);call require(good,'thermo')
  y=0;y(1)=.1_dp;y(9)=.9_dp
  call mixture_mass_properties(species,y,2000._dp,mw,r,cp,cv,gamma,h,e,s,good)
  call require(good,'initial energy')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),good)
  call require(good,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,good);call require(good,'tree')
  owners=[0,1];dt=2e-8_dp
  allocate(gas%fields(2),temp%fields(2),held%fields(2),ir%fields(2))
  allocate(local_gas%fields(2),local_temp%fields(2),local_held%fields(2),local_ir%fields(2))
  do node=1,2
    n=tree%nodes(node)%shape
    allocate(gas%fields(node)%values(14,n(1),n(2),n(3)),source=0._dp)
    allocate(temp%fields(node)%values(1,n(1),n(2),n(3)),source=2000._dp)
    allocate(held%fields(node)%values(1,n(1),n(2),n(3)),source=.4_dp*e)
    allocate(ir%fields(node)%values(10,n(1),n(2),n(3)),source=9._dp)
    gas%fields(node)%values(1,:,:,:)=.4_dp;gas%fields(node)%values(5,:,:,:)=.4_dp*e
    gas%fields(node)%values(6,:,:,:)=.04_dp;gas%fields(node)%values(14,:,:,:)=.36_dp
    if(owners(node)/=rank) cycle
    local_gas%fields(node)=gas%fields(node);local_temp%fields(node)=temp%fields(node)
    local_held%fields(node)=held%fields(node);local_ir%fields(node)=ir%fields(node)
  enddo
  do cycle=1,2
    calls=0;sum_dt=0;reject_after=0
    call initialize_amr_tree_reaction_source_3d(tree,species,reactions,gas,temp,dt,cycle==2,1e-8_dp,1e-14_dp,ir,good, &
      evaluate_reaction=reactor,internal_energy=held)
    call require(good.and.calls==320,'all valid cells initialized once')
    expected=64*dt+256*dt/real(cycle,dp)
    call require(abs(sum_dt-expected)<=1e-13_dp*expected,'level-specific initial reaction duration')
    call require(all(ir%fields(1)%values(:,1:2,:,:)==0),'pinned root restriction precedes fine initialization')
    call require(maxval(abs(ir%fields(2)%values(1,:,:,:)))>0,'nonzero chemical estimate')
    calls=0;sum_dt=0
    call initialize_mpi_amr_tree_reaction_source_3d(tree,owners,species,reactions,local_gas,local_temp,dt,cycle==2, &
      1e-8_dp,1e-14_dp,local_ir,MPI_COMM_WORLD,good,evaluate_reaction=reactor,internal_energy=local_held)
    call require(good,'collective initialization')
    call require((rank==2).eqv.(calls==0),'idle rank has no reactor calls')
    do node=1,2
      call require(all(gas%fields(node)%values(1,:,:,:)==.4_dp).and. &
        all(gas%fields(node)%values(6,:,:,:)==.04_dp).and.all(temp%fields(node)%values==2000), 'public state unchanged')
      if(owners(node)/=rank) cycle
      call require(all(local_ir%fields(node)%values==ir%fields(node)%values),'exact serial/MPI history')
      call require(all(local_gas%fields(node)%values==gas%fields(node)%values).and. &
        all(local_temp%fields(node)%values==temp%fields(node)%values).and. &
        all(local_held%fields(node)%values==held%fields(node)%values),'gas/T/UEINT unchanged')
    enddo
    saved=ir;calls=0;reject_after=65
    call initialize_amr_tree_reaction_source_3d(tree,species,reactions,gas,temp,dt,cycle==2,1e-8_dp,1e-14_dp,ir,good, &
      evaluate_reaction=reactor,internal_energy=held)
    call require(.not.good,'late fine failure refused')
    do node=1,2
      call require(all(ir%fields(node)%values==saved%fields(node)%values),'whole-tree source rollback')
    enddo
    saved=local_ir;calls=0;reject_after=0
    if(rank==1) reject_after=1
    call initialize_mpi_amr_tree_reaction_source_3d(tree,owners,species,reactions,local_gas,local_temp,dt,cycle==2, &
      1e-8_dp,1e-14_dp,local_ir,MPI_COMM_WORLD,good,evaluate_reaction=reactor,internal_energy=local_held)
    call require(.not.good,'remote fine failure refused collectively')
    do node=1,2
      if(owners(node)/=rank) cycle
      call require(all(local_ir%fields(node)%values==saved%fields(node)%values),'collective source rollback')
    enddo
  enddo
  if(rank==0) print *, 'PASS: initial MOL level clocks, nonzero source, immutable gas, idle rank and late rollback'
  call MPI_Finalize(status)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(condition) return
    print *, 'FAIL rank',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
  subroutine reactor(sp,rx,energy,external_y,external_e,duration,rtol,atol,densities,t,count,ok,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: energy,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),t
    integer, intent(out) :: count
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    real(dp) :: rho,delta,recovered
    calls=calls+1;sum_dt=sum_dt+duration
    call require(all(external_y==0).and.external_e==0,'initialization has no spatial forcing')
    call require(energy==.4_dp*e,'actual internal energy supplied')
    if(reject_after>0.and.calls==reject_after) then
      densities=-1;t=-1;count=0;ok=.false.;return
    endif
    ! Synthetic H2 -> 2H transfer conserves mass/elements; its dt dependence
    ! makes an incorrect fine-level initialization interval observable.
    rho=sum(densities);delta=1e6_dp*duration**2
    densities(1)=densities(1)-delta;densities(4)=densities(4)+delta
    call temperature_from_internal_energy(sp,densities/rho,energy/rho,t,recovered,ok)
    if(ok) t=recovered
    count=1
  end subroutine
end program test_amr_tree_initial_reaction
