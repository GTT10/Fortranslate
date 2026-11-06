program test_regular_held_diffusion
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use elementary_kinetics_mod, only: elementary_reaction
  use mixture_thermo_mod, only: mixture_mass_properties
  use dual_energy_mod, only: dual_energy_policy
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_sdc_3d,advance_mpi_amr_tree_sdc_subcycle_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: compute_sdc_mpi_amr_tree_transport_3d, &
    compute_mpi_amr_tree_transport_timestep_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d
  use amr_patch_tree_subcycle_3d_mod, only: compute_amr_patch_tree_subcycle_cfl_3d
  use mpi_amr_patch_tree_subcycle_3d_mod, only: compute_mpi_amr_tree_subcycle_cfl_3d
  use amr_patch_tree_transport_3d_mod, only: compute_amr_patch_tree_transport_timestep_3d
  use amr_patch_tree_transport_subcycle_3d_mod, only: compute_amr_patch_tree_subcycle_transport_timestep_3d
  use mpi_amr_patch_tree_transport_subcycle_3d_mod, only: compute_mpi_amr_tree_subcycle_transport_timestep_3d
  use les_mpi_amr_tree_mod, only: compute_mpi_amr_tree_les_rhs_3d,compute_mpi_amr_tree_les_timestep_3d
  implicit none
  type :: trial
    type(amr_tree_state_3d) :: gas,temp,held
    real(dp), allocatable :: boundary(:),source(:)
    real(dp) :: prefix,correction,theta
  end type
  type(trial) :: run(2),saved
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: boxes(1),empty(0)
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction) :: reactions(0)
  type(dual_energy_policy) :: policy
  type(amr_tree_state_3d) :: rates(2),first_callback
  type(amr_tree_flux_3d), allocatable :: flux_a(:),flux_b(:)
  integer, allocatable :: owners(:)
  real(dp), allocatable :: y(:)
  real(dp) :: spacing(3),dt,bound(2),diffusion(2),theta(2),local_change,sample_change
  integer :: rank,nrank,status,nv,kind,a,node,n(3),i,j,k,c,model_case,calls
  character(len=16) :: faces(6),model
  logical :: ok,reject_callback
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,status)
  call require(nrank==3,'two owners and idle')
  call load_h2o2_full_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_full_transport(transport,ok);call require(ok,'transport')
  nv=5+size(species);allocate(y(size(species)));y=0
  do c=1,size(species)
    if(species(c)%name=='N2') y(c)=1
  enddo
  faces='periodic';spacing=.001_dp;dt=2e-8_dp
  policy%allow_negative_energy=.false.;policy%eta2=2;policy%update_total_from_internal=.false.
  do kind=1,2
    if(kind==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer ::],empty,tree,ok)
      owners=[0]
    else
      call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),ok)
      call require(ok,'fine patch')
      call initialize_amr_patch_tree_3d([4,4,4],[1],boxes,tree,ok)
      owners=[0,1]
    endif
    call require(ok,'tree')
    do a=1,2
      call initialize(run(a),a==1)
    enddo
    do a=1,2
      call compute_mpi_amr_tree_subcycle_cfl_3d(tree,owners,species,run(a)%gas,run(a)%temp,spacing,.2_dp, &
        bound(a),MPI_COMM_WORLD,ok,internal_energy=run(a)%held)
      call require(ok,'recursive actualH hydro bound')
    enddo
    call require(bound(1)==bound(2),'recursive hydro bound independent of Eden')
    do a=1,2
      call compute_mpi_amr_tree_subcycle_transport_timestep_3d(tree,owners,species,transport,run(a)%gas,run(a)%temp, &
        spacing,.1_dp,.true.,.true.,.true.,bound(a),diffusion(a),MPI_COMM_WORLD,ok, &
        faces=faces,internal_energy=run(a)%held)
      call require(ok,'recursive actualH transport bound')
    enddo
    call require(bound(1)==bound(2).and.diffusion(1)==diffusion(2),'recursive transport bound independent of Eden')
    if(kind==1.and.rank==0) then
      do model_case=1,4
        do a=1,2
          select case(model_case)
          case(1)
            call compute_amr_patch_tree_cfl_3d(tree,species,run(a)%gas,run(a)%temp,spacing,.2_dp,bound(a),ok,run(a)%held)
          case(2)
            call compute_amr_patch_tree_subcycle_cfl_3d(tree,species,run(a)%gas,run(a)%temp, &
              spacing,.2_dp,bound(a),ok,run(a)%held)
          case(3)
            call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,run(a)%gas,run(a)%temp,spacing, &
              .1_dp,.true.,.true.,.true.,bound(a),diffusion(a),ok,faces=faces,internal_energy=run(a)%held)
          case(4)
            call compute_amr_patch_tree_subcycle_transport_timestep_3d(tree,species,transport,run(a)%gas,run(a)%temp, &
              spacing,.1_dp,.true.,.true.,.true.,bound(a),diffusion(a),ok,faces=faces,internal_energy=run(a)%held)
          end select
          call require(ok,'serial actualH timestep')
        enddo
        call require(bound(1)==bound(2),'serial common/recursive hydro/transport bound independent of Eden')
      enddo
    endif
    if(kind==1) then
      ! Operators must see the actual nonuniform 700K family, independently of
      ! the distinct bulk energy. Dynamic must filter actual H too.
      do model_case=1,4
        select case(model_case)
        case(1);model='smagorinsky'
        case(2);model='wale'
        case(3);model='vreman'
        case(4);model='dynamic'
        end select
        do a=1,2
          if(a==1) then
            call les(run(a)%gas,run(a)%temp,-1,rates(a),flux_a,diffusion(a),ok,run(a)%held)
          else
            call les(run(a)%gas,run(a)%temp,-1,rates(a),flux_b,diffusion(a),ok,run(a)%held)
          endif
          call require(ok,'LES actualH operator')
          call compute_mpi_amr_tree_les_timestep_3d(tree,owners,species,run(a)%gas,run(a)%temp,spacing,.1_dp,faces, &
            .17_dp,.7_dp,.7_dp,bound(a),diffusion(a),MPI_COMM_WORLD,ok,model=model,ci=.1_dp, &
            internal_energy=run(a)%held)
          call require(ok,'LES actualH timestep')
        enddo
        call require(bound(1)==bound(2).and.diffusion(1)==diffusion(2),'LES independent actualH bound')
        call compare_rates()
        call require_any_nonzero('LES spatial family reached')
      enddo
      do a=1,2
        theta(a)=1
        if(a==1) then
          call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,run(a)%gas,run(a)%temp,spacing,dt, &
            faces,[.true.,.true.,.true.,.false.],rates(a),flux_a,theta(a),MPI_COMM_WORLD,ok, &
            internal_energy=run(a)%held)
        else
          call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,run(a)%gas,run(a)%temp,spacing,dt, &
            faces,[.true.,.true.,.true.,.false.],rates(a),flux_b,theta(a),MPI_COMM_WORLD,ok, &
            internal_energy=run(a)%held)
        endif
        call require(ok,'transport actualH operator')
        call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,run(a)%gas,run(a)%temp,spacing, &
          .1_dp,.true.,.true.,.true.,bound(a),diffusion(a),MPI_COMM_WORLD,ok,internal_energy=run(a)%held)
        call require(ok,'transport actualH timestep')
      enddo
      call require(bound(1)==bound(2).and.diffusion(1)==diffusion(2),'transport independent actualH bound')
      call compare_rates();call require_any_nonzero('transport spatial family reached')
    endif
    model='vreman';reject_callback=.false.
    do a=1,2
      calls=0;sample_change=0
      call advance(run(a),ok)
      call require(ok,'native actualH molecular plus LES interval')
      call require(calls>=2.and.sample_change>1e-8_dp,'callback receives distinct native old/new H')
    enddo
    local_change=0
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      call require(all(abs(run(1)%held%fields(node)%values-run(2)%held%fields(node)%values)<= &
        128*epsilon(1._dp)*max(1._dp,abs(run(2)%held%fields(node)%values))),'native actualH paired endpoint')
      call require(all(abs(run(1)%temp%fields(node)%values-run(2)%temp%fields(node)%values)<= &
        128*epsilon(1._dp)*max(1._dp,abs(run(2)%temp%fields(node)%values))),'native actualH paired EOS')
    enddo
    call require(run(1)%prefix==3.5_dp.and.run(2)%prefix==3.5_dp,'native numerical prefix unchanged')
    saved=run(1);reject_callback=.true.;calls=0
    call advance(run(1),ok)
    call require(.not.ok.and.calls==2,'late diffusion callback refusal')
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      call require(all(run(1)%gas%fields(node)%values==saved%gas%fields(node)%values).and. &
        all(run(1)%temp%fields(node)%values==saved%temp%fields(node)%values).and. &
        all(run(1)%held%fields(node)%values==saved%held%fields(node)%values),'late callback field rollback')
    enddo
    call require(all(run(1)%boundary==saved%boundary).and.all(run(1)%source==saved%source).and. &
      run(1)%prefix==saved%prefix,'late callback history rollback')
    if(rank==0) print *, 'PASS: actualH molecular/LES operators, native stage chronology and rollback, case',kind
  enddo
  call MPI_Finalize(status)
contains
  subroutine require(ready,message)
    logical, intent(in) :: ready
    character(len=*), intent(in) :: message
    if(ready) return
    print *, 'FAIL rank',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
  subroutine initialize(t,bulk_hot)
    type(trial), intent(inout) :: t
    logical, intent(in) :: bulk_hot
    real(dp) :: x(3),v(3),thermal,mw,rg,cp,cv,gamma,h,e,entropy,bulk_e,dx
    if(allocated(t%gas%fields)) deallocate(t%gas%fields,t%temp%fields,t%held%fields,t%boundary,t%source)
    allocate(t%gas%fields(size(owners)),t%temp%fields(size(owners)),t%held%fields(size(owners)))
    allocate(t%boundary(nv),t%source(nv));t%boundary=0;t%source=0;t%prefix=3.5_dp;t%theta=1;t%correction=0
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,bulk_e,entropy,ok)
    call require(ok,'bulk thermo')
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;dx=.001_dp/real(2**tree%nodes(node)%level,dp)
      allocate(t%gas%fields(node)%values(nv,n(1),n(2),n(3)),t%temp%fields(node)%values(1,n(1),n(2),n(3)), &
        t%held%fields(node)%values(1,n(1),n(2),n(3)))
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        x=(real([i,j,k],dp)-.5_dp)*dx
        thermal=700+25*sin(2*acos(-1._dp)*x(1)/.004_dp)
        v=[8*sin(2*acos(-1._dp)*x(1)/.004_dp),2*sin(2*acos(-1._dp)*x(2)/.004_dp),0._dp]
        call mixture_mass_properties(species,y,thermal,mw,rg,cp,cv,gamma,h,e,entropy,ok)
        call require(ok,'actual thermo')
        t%gas%fields(node)%values(:,i,j,k)=0
        t%gas%fields(node)%values(1,i,j,k)=.4_dp;t%gas%fields(node)%values(2:4,i,j,k)=.4_dp*v
        t%gas%fields(node)%values(5,i,j,k)=.4_dp*merge(bulk_e,e,bulk_hot)+.2_dp*sum(v*v)
        t%gas%fields(node)%values(6:,i,j,k)=.4_dp*y
        t%temp%fields(node)%values(1,i,j,k)=thermal;t%held%fields(node)%values(1,i,j,k)=.4_dp*e
      enddo
      enddo
      enddo
    enddo
  end subroutine
  subroutine les(u,t,depth,rates,flux,diffusion,good,held,fractions,coarse_start,coarse_end)
    type(amr_tree_state_3d), intent(in) :: u,t,held
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    integer, intent(in) :: depth
    real(dp), optional, intent(in) :: fractions(:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(out) :: diffusion
    logical, intent(out) :: good
    diffusion=0
    call compute_mpi_amr_tree_les_rhs_3d(tree,owners,species,u,t,spacing,faces,.17_dp,.7_dp,.7_dp, &
      rates,flux,diffusion,MPI_COMM_WORLD,good,model=model,ci=.1_dp,target_depth=depth, &
      fractions=fractions,coarse_start=coarse_start,coarse_end=coarse_end,internal_energy=held)
  end subroutine
  subroutine callback(u,t,depth,fractions,coarse_start,coarse_end,rates,flux,good,internal_energy)
    type(amr_tree_state_3d), intent(in) :: u,t
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,coarse_start,coarse_end
    integer, intent(in) :: depth
    real(dp), optional, intent(in) :: fractions(:)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: good
    real(dp) :: diffusion,local
    integer :: peer
    calls=calls+1
    call require(present(internal_energy),'typed callback actualH presence')
    if(calls==1) first_callback=internal_energy
    local=0
    do peer=1,size(owners)
      if(owners(peer)/=rank) cycle
      if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
      local=max(local,maxval(abs(internal_energy%fields(peer)%values-first_callback%fields(peer)%values)))
      if(present(coarse_start)) then
        if(tree%nodes(peer)%level<depth) call require(size(coarse_start%fields(peer)%values,1)==nv+2, &
          'typed callback full actualH ancestor packet')
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
    sample_change=max(sample_change,local)
    call les(u,t,depth,rates,flux,diffusion,good,internal_energy,fractions,coarse_start,coarse_end)
    if(reject_callback.and.calls==2) good=.false.
  end subroutine
  subroutine advance(t,good)
    type(trial), intent(inout) :: t
    logical, intent(out) :: good
    if(kind==1) then
      call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,t%gas,t%temp,spacing,dt,faces, &
        2,.false.,[.true.,.true.,.true.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
        t%boundary,t%source,t%correction,t%theta,MPI_COMM_WORLD,good,evaluate_diffusion=callback,mol_mode=.true., &
        internal_energy=t%held,energy_correction=t%prefix,energy_policy=policy)
    else
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,t%gas,t%temp,spacing,dt,faces, &
        2,.false.,[.true.,.true.,.true.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
        t%boundary,t%source,t%correction,t%theta,MPI_COMM_WORLD,good,evaluate_diffusion=callback,mol_mode=.true., &
        internal_energy=t%held,energy_correction=t%prefix,energy_policy=policy)
    endif
  end subroutine
  subroutine compare_rates()
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      call require(all(rates(1)%fields(node)%values==rates(2)%fields(node)%values),'actualH physical rates independent of Eden')
      call require(all(flux_a(node)%x==flux_b(node)%x).and.all(flux_a(node)%y==flux_b(node)%y).and. &
        all(flux_a(node)%z==flux_b(node)%z),'actualH physical faces independent of Eden')
    enddo
  end subroutine
  subroutine require_any_nonzero(message)
    character(len=*), intent(in) :: message
    local_change=0
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      local_change=max(local_change,maxval(abs(rates(1)%fields(node)%values)))
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,local_change,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
    call require(status==MPI_SUCCESS.and.local_change>1e-6_dp,message)
  end subroutine
end program
