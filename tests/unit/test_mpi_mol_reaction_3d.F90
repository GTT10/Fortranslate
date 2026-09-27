program test_mpi_mol_reaction_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: mixture_mass_properties
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,extract_eb_geometry_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(0)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: geometry(1),extended(1)
  type(amr_tree_state_3d) :: initial,initial_t,state,temp,saved,saved_t
  real(dp), allocatable :: y(:),q(:),u(:),reference(:),b(:),w(:),ch(:),source(:)
  real(dp) :: mw,r,cp,cv,gamma,h,e,entropy,t,c,reference_t,theta,error
  integer :: rank,status,nv,s,iteration,steps,i,j,k,calls,global_calls
  logical :: ok,fail_late
  character(len=16) :: faces(6)
  call MPI_Init(status);call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok);call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(y(size(species)),q(reactive_nprim(size(species))),u(nv),reference(size(species)), &
    b(nv),w(nv),ch(nv),source(nv))
  y=0
  do s=1,size(species)
    select case(species(s)%name)
    case('H2');y(s)=.1_dp
    case('O2');y(s)=.2_dp
    case('N2');y(s)=.699_dp
    case('H');y(s)=.001_dp
    end select
  end do
  call mixture_mass_properties(species,y,1200._dp,mw,r,cp,cv,gamma,h,e,entropy,ok)
  call require(ok,'properties')
  q=0;q(1)=1;q(5)=r*1200;q(6:)=y
  call reactive_primitive_to_conserved(species,q,u,t,c,ok);call require(ok,'state')
  reference=u(6:);reference_t=t
  call advance_forced_reactor_interval(species,reactions,u(5),0*y,0._dp,1e-10_dp,1e-8_dp,1e-12_dp, &
    reference,reference_t,steps,ok)
  call require(ok.and.maxval(abs(reference-u(6:)))>1e-10_dp,'nontrivial homogeneous reactor')
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],boxes,tree,ok);call require(ok,'tree')
  call build_axis_plane_eb_geometry_3d(14,14,14,-3._dp,4._dp,-3._dp,4._dp,-3._dp,4._dp, &
    'x',-10._dp,extended(1),ok);call require(ok,'extended regular geometry')
  call extract_eb_geometry_3d(extended(1),[7,7,7],[8,8,8],geometry(1),ok);call require(ok,'geometry')
  allocate(initial%fields(1),initial_t%fields(1))
  if(rank==0) then
    allocate(initial%fields(1)%values(nv,2,2,2),initial_t%fields(1)%values(1,2,2,2))
    do k=1,2
    do j=1,2
    do i=1,2
      initial%fields(1)%values(:,i,j,k)=u
    end do
    end do
    end do
    initial_t%fields(1)%values=t
  end if
  faces='periodic';fail_late=.false.
  do iteration=1,2
    state=initial;temp=initial_t;calls=0
    call step(iteration,ok);call require(ok,'reacting MOL step')
    call MPI_Allreduce(calls,global_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(global_calls==8*iteration,'MOL reaction iteration count')
    error=0
    if(rank==0) then
      do k=1,2
      do j=1,2
      do i=1,2
        error=max(error,maxval(abs(state%fields(1)%values(6:,i,j,k)-reference)))
        error=max(error,abs(temp%fields(1)%values(1,i,j,k)-reference_t)/reference_t)
      end do
      end do
      end do
      call require(error<1e-11_dp,'homogeneous reactor agreement')
      call require(maxval(abs(state%fields(1)%values(:,1,1,1)-u-b-w-ch-source))<1e-9_dp,'chemical ledger')
      print *, 'MOL iterations, reactor calls, maximum error:',iteration,global_calls,error
    end if
  end do
  saved=state;saved_t=temp;calls=0;fail_late=.true.;b=17;w=19;ch=23;source=29
  call step(2,ok);call require(.not.ok,'late reactor failure')
  call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(source==29),'ledger rollback')
  if(rank==0) then
    call require(all(state%fields(1)%values==saved%fields(1)%values),'state rollback')
    call require(all(temp%fields(1)%values==saved_t%fields(1)%values),'temperature rollback')
    print *, 'PASS: reacting MOL and late-stage transactional rollback'
  end if
  call MPI_Finalize(status)
contains
  subroutine step(iterations,good)
    integer, intent(in) :: iterations
    logical, intent(out) :: good
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
      1e-10_dp,iterations,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
      1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,good,evaluate_reaction=reactor,mol_mode=.true.)
  end subroutine
  subroutine reactor(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    calls=calls+1
    if(fail_late.and.calls==9) then
      species_density=-1;temperature=-1;completed_steps=0;ok=.false.;return
    end if
    call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps=maximum_steps)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: code
    if(condition) return
    print *, 'FAIL ',label
    call MPI_Abort(MPI_COMM_WORLD,1,code)
  end subroutine
end program
