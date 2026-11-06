! Matched-state probe for the two-level periodic reacting TG control.
! This exercises the production tree integrator with imported SI cell states.
! LiDryer CVODE chemistry; optional pinned Simple transport, no clipping.
program amr_tree_reacting_probe
#ifdef PELEF_REFERENCE_MPI
  use mpi_f08
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_sdc_3d,advance_mpi_amr_tree_sdc_subcycle_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: initialize_mpi_amr_tree_reaction_source_3d
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,configure_nasa7_eos_policy
  use lidryer_mechanism_mod, only: load_lidryer_thermo,load_lidryer_mechanism,lidryer_nspecies,load_lidryer_transport
  use sundials_forced_reactor_mod, only: advance_forced_reactor_cvode
  use dual_energy_mod, only: dual_energy_policy
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species,load_gas_transport_data,set_simple_transport_fits
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_sdc_3d,advance_amr_patch_tree_sdc_subcycle_3d
  use amr_patch_tree_chemistry_3d_mod, only: initialize_amr_tree_reaction_source_3d
  implicit none
#ifdef PELEF_REFERENCE_MPI
  integer :: rank,nranks,ierr,owners(2)
  type(amr_tree_state_3d) :: output_state,output_temp,output_held,output_history
#endif
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,seen,held,history
  type(dual_energy_policy) :: policy
  integer, parameter :: nv=5+lidryer_nspecies
  real(dp) :: spacing(3),dt,budget(nv),source(nv),correction,theta,values(2*nv-2),force(3),energy_correction
  real(dp) :: boundary_integral(nv),source_integral(nv)
  character(len=8),parameter :: budget_names(5)=[character(len=8)::"density","xmom","ymom","zmom","rho_E"]
  integer :: unit,ios,steps,iterations,subcycle,node,i,j,k,records,row,n(3),step,c,root_cells,transport_flag
  real(dp) :: transport_values(5,lidryer_nspecies),fits(4,2*lidryer_nspecies+lidryer_nspecies*(lidryer_nspecies+1)/2)
  integer :: geometries(lidryer_nspecies)
  character(len=24) :: transport_names(lidryer_nspecies)
  logical :: isothermal_slip_faces(6)
  real(dp), allocatable :: wall_temperature(:)
  real(dp) :: parsed_wall_temperature
  logical :: ok,transport_on,reaction_on,slip_x,no_slip_x,outflow_x
  character(len=1024) :: input,prefix,path
  character(len=64) :: magic
  character(len=32) :: faces(6)
  if(command_argument_count()<2.or.command_argument_count()>6) &
    error stop 'Usage: probe input.dat prefix [--initialize-history] [--reaction-off] [--slip-x|--no-slip-x|--outflow-x]'
#ifdef PELEF_REFERENCE_MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  owners=mod([0,1],nranks)
#endif
  reaction_on=.true.;slip_x=.false.;no_slip_x=.false.;outflow_x=.false.
  do c=3,command_argument_count()
    call get_command_argument(c,magic)
    if(index(trim(magic),'--wall-temperature=')==1) then
      if(allocated(wall_temperature)) error stop 'Duplicate wall temperature'
      read(magic(20:),*,iostat=ios) parsed_wall_temperature
      if(ios/=0) error stop 'Invalid wall temperature'
      if(.not.ieee_is_finite(parsed_wall_temperature).or.parsed_wall_temperature<=0) &
        error stop 'Invalid wall temperature'
      allocate(wall_temperature(6));wall_temperature=-1
      wall_temperature(1:2)=parsed_wall_temperature
      cycle
    endif
    if(trim(magic)=='--reaction-off') reaction_on=.false.
    if(trim(magic)=='--slip-x') slip_x=.true.
    if(trim(magic)=='--no-slip-x') no_slip_x=.true.
    if(trim(magic)=='--outflow-x') outflow_x=.true.
    if(trim(magic)/='--reaction-off'.and.trim(magic)/='--initialize-history'.and.trim(magic)/='--slip-x'.and. &
      trim(magic)/='--no-slip-x'.and.trim(magic)/='--outflow-x') &
      error stop 'Unknown probe mode'
  enddo
  call get_command_argument(1,input);call get_command_argument(2,prefix)
  open(newunit=unit,file=trim(input),status='old',action='read',iostat=ios)
  if(ios/=0) error stop 'Cannot open probe input'
  read(unit,*,iostat=ios) magic
  if(ios/=0.or.(magic/='PELEF_REACTING_AMR_REFERENCE_V2'.and. &
    magic/='PELEF_REACTING_AMR_REFERENCE_V3'.and.magic/='PELEF_REACTING_AMR_REFERENCE_V4')) &
    error stop 'Invalid probe input magic'
  transport_on=magic=='PELEF_REACTING_AMR_REFERENCE_V3'
  root_cells=4
  if(magic=='PELEF_REACTING_AMR_REFERENCE_V4') then
    read(unit,*,iostat=ios) root_cells,transport_flag
    if(ios/=0) error stop 'Invalid mesh controls'
    if(.not.any(root_cells==[4,8,16,32,64]).or.(transport_flag/=0.and.transport_flag/=1)) &
      error stop 'Unsupported mesh/transport controls'
    transport_on=transport_flag==1
  endif
  read(unit,*,iostat=ios) steps,iterations,subcycle,dt
  if(ios/=0) error stop 'Invalid stepping row'
  if(steps<0.or.steps>100.or.iterations<1.or.(subcycle/=0.and.subcycle/=1).or.dt<=0) error stop 'Invalid stepping controls'
  read(unit,*,iostat=ios) spacing
  if(ios/=0.or.any(spacing<=0)) error stop 'Invalid probe spacing'
  if(transport_on) then
    do row=1,size(fits,2)
      read(unit,*,iostat=ios) fits(:,row)
      if(ios/=0) error stop 'Invalid Simple transport fits'
    enddo
  endif
  read(unit,*,iostat=ios) records
  if(ios/=0.or.records/=5*root_cells**3) error stop 'Probe requires all coarse and fine cells'
  call load_lidryer_thermo(species,ok)
  if(.not.ok) error stop 'Cannot load LiDryer thermo'
  call load_lidryer_mechanism(reactions,ok)
  if(.not.ok) error stop "Cannot load LiDryer reactions"
  call configure_nasa7_eos_policy(species,.true.,ok)
  if(.not.ok) error stop "Invalid EOS policy"
  policy%allow_negative_energy=.false.;policy%eta2=2
  energy_correction=0;boundary_integral=0;source_integral=0
  if(transport_on) then
    call load_lidryer_transport(transport_names,geometries,transport_values,ok)
    if(.not.ok) error stop 'Cannot load LiDryer transport'
    call load_gas_transport_data(transport_names,geometries,transport_values,transport,ok)
    if(.not.ok) error stop 'Invalid LiDryer transport'
    call set_simple_transport_fits(transport,fits,ok)
    if(.not.ok) error stop 'Invalid pinned Simple transport'
  else
    allocate(transport(0))
  endif
  call initialize_amr_patch_3d(root_cells,root_cells,root_cells,1,root_cells/2,1,root_cells,1,root_cells,2,patch(1),ok)
  if(.not.ok) error stop 'Invalid reference patch'
  call initialize_amr_patch_tree_3d([root_cells,root_cells,root_cells],[1],patch,tree,ok)
  if(.not.ok) error stop 'Invalid reference tree'
  allocate(state%fields(2),temperature%fields(2),seen%fields(2),held%fields(2),history%fields(2))
  do node=1,2
    n=tree%nodes(node)%shape
    allocate(state%fields(node)%values(nv,n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)), &
      seen%fields(node)%values(1,n(1),n(2),n(3)))
    allocate(held%fields(node)%values(1,n(1),n(2),n(3)),history%fields(node)%values(nv-4,n(1),n(2),n(3)))
    history%fields(node)%values=0;held%fields(node)%values=0
    state%fields(node)%values=0;temperature%fields(node)%values=0;seen%fields(node)%values=0
  enddo
  do row=1,records
    read(unit,*,iostat=ios) node,i,j,k,values
    if(ios/=0.or.node<1.or.node>2) error stop 'Invalid imported cell'
    n=tree%nodes(node)%shape
    if(any([i,j,k]<1).or.any([i,j,k]>n)) error stop 'Imported cell index outside patch'
    if(seen%fields(node)%values(1,i,j,k)/=0) error stop 'Duplicate imported cell'
    seen%fields(node)%values(1,i,j,k)=1
    state%fields(node)%values(:,i,j,k)=values(:nv)
    temperature%fields(node)%values(1,i,j,k)=values(nv+2)
    held%fields(node)%values(1,i,j,k)=values(nv+1)
    ! Native post_init seeds the old reaction estimate before the first MOL
    ! predictor. It is part of the imported initial state, not future forcing.
    history%fields(node)%values(:,i,j,k)=values(nv+3:)
  enddo
  read(unit,*,iostat=ios) magic
  if(ios>=0) error stop 'Unexpected records after imported cells'
  close(unit)
  do node=1,2
    if(any(seen%fields(node)%values/=1)) error stop 'Missing imported cell'
  enddo
  if(allocated(wall_temperature)) then
    if((.not.no_slip_x.and..not.slip_x).or..not.transport_on) error stop 'Wall temperature requires a conducting wall'
  endif
  faces='periodic';force=0;isothermal_slip_faces=.false.
  if(count([slip_x,no_slip_x,outflow_x])>1) error stop 'Conflicting x boundary modes'
  if(slip_x) faces(1:2)='slip_wall'
  if(no_slip_x) faces(1:2)='no_slip_wall'
  if(allocated(wall_temperature)) then
    if(slip_x) then
      isothermal_slip_faces(1:2)=.true.
    else
      faces(1:2)='isothermal_wall'
    endif
  endif
  if(outflow_x) faces(1:2)='outflow'
  if(.not.reaction_on) then
    do node=1,2
      if(any(history%fields(node)%values/=0)) error stop 'Inactive chemistry requires zero imported history'
    enddo
  endif
#ifdef PELEF_REFERENCE_MPI
  output_state=state;output_temp=temperature;output_held=held;output_history=history
  do node=1,2
    if(owners(node)==rank) cycle
    deallocate(state%fields(node)%values,temperature%fields(node)%values, &
      held%fields(node)%values,history%fields(node)%values)
  enddo
#endif
  do c=3,command_argument_count()
    call get_command_argument(c,magic)
    if(trim(magic)/='--initialize-history') cycle
    if(.not.reaction_on) error stop 'Cannot initialize disabled chemistry'
#ifdef PELEF_REFERENCE_MPI
    call initialize_mpi_amr_tree_reaction_source_3d(tree,owners,species,reactions,state,temperature,dt,subcycle==1, &
      1e-10_dp,1e-14_dp,history,MPI_COMM_WORLD,ok,evaluate_reaction=advance_forced_reactor_cvode,internal_energy=held)
#else
    call initialize_amr_tree_reaction_source_3d(tree,species,reactions,state,temperature,dt,subcycle==1, &
      1e-10_dp,1e-14_dp,history,ok,evaluate_reaction=advance_forced_reactor_cvode,internal_energy=held)
#endif
    if(.not.ok) error stop 'Initial reaction estimate refused'
  enddo
  do step=0,steps
#ifdef PELEF_REFERENCE_MPI
    call collect(state,output_state)
    call collect(temperature,output_temp)
    call collect(held,output_held)
    call collect(history,output_history)
    if(rank==0) then
#endif
    write(path,'(a,"_",i3.3,".csv")') trim(prefix),step
    open(newunit=unit,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open probe output'
    write(unit,'(a)',advance='no') 'node,level,i,j,k,time_s,density,xmom,ymom,zmom,rho_E'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rho_'//trim(species(c)%name)
    enddo
    write(unit,'(a)',advance='no') ',rho_e,Temp'
    do c=1,size(species)
      write(unit,'(a)',advance='no') ',rho_omega_'//trim(species(c)%name)
    enddo
    write(unit,'(a)') ',rhoe_dot'
    do node=1,2
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        write(unit,'(5(i0,","),*(es25.16e3,:,","))') node,tree%nodes(node)%level,i,j,k,step*dt, &
#ifdef PELEF_REFERENCE_MPI
          output_state%fields(node)%values(:,i,j,k),output_held%fields(node)%values(1,i,j,k), &
          output_temp%fields(node)%values(1,i,j,k),output_history%fields(node)%values(:,i,j,k)
#else
          state%fields(node)%values(:,i,j,k),held%fields(node)%values(1,i,j,k),temperature%fields(node)%values(1,i,j,k), &
          history%fields(node)%values(:,i,j,k)
#endif
      enddo
      enddo
      enddo
    enddo
    close(unit)
    write(path,'(a,"_",i3.3,".budget.csv")') trim(prefix),step
    open(newunit=unit,file=trim(path),status="replace",action="write")
    write(unit,'(a)',advance='no') "step,time_s,numerical_energy_J"
    do c=1,5
      write(unit,'(a)',advance='no') ',boundary_'//trim(budget_names(c))//',source_'//trim(budget_names(c))
    enddo
    if(outflow_x) then
      do c=1,lidryer_nspecies
        write(unit,'(a)',advance='no') ',boundary_rho_'//trim(species(c)%name)//',source_rho_'//trim(species(c)%name)
      enddo
    endif
    write(unit,*)
    write(unit,'(i0,2(",",es25.16e3))',advance='no') step,step*dt,energy_correction
    do c=1,5
      write(unit,'(2(",",es25.16e3))',advance='no') boundary_integral(c),source_integral(c)
    enddo
    if(outflow_x) then
      do c=6,nv
        write(unit,'(2(",",es25.16e3))',advance='no') boundary_integral(c),source_integral(c)
      enddo
    endif
    write(unit,*)
    close(unit)
#ifdef PELEF_REFERENCE_MPI
    endif
#endif
    if(step==steps) exit
    if(.not.reaction_on.and.subcycle==1) then
#ifdef PELEF_REFERENCE_MPI
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,.false.,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#else
      call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,.false.,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,ok,mol_mode=.true.,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#endif
    else if(.not.reaction_on) then
#ifdef PELEF_REFERENCE_MPI
      call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,.false.,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#else
      call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,.false.,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,ok,mol_mode=.true.,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#endif
    else if(subcycle==1) then
#ifdef PELEF_REFERENCE_MPI
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,reaction_on,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,reaction_source=history, &
        evaluate_reaction=advance_forced_reactor_cvode,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#else
      call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,reaction_on,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,ok,mol_mode=.true.,reaction_source=history, &
        evaluate_reaction=advance_forced_reactor_cvode,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#endif
    else
#ifdef PELEF_REFERENCE_MPI
      call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,reaction_on,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,MPI_COMM_WORLD,ok,mol_mode=.true.,reaction_source=history, &
        evaluate_reaction=advance_forced_reactor_cvode,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#else
      call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,reaction_on,[transport_on,transport_on,transport_on,transport_on],force,1e-10_dp,1e-14_dp, &
        budget,source,correction,theta,ok,mol_mode=.true.,reaction_source=history, &
        evaluate_reaction=advance_forced_reactor_cvode,internal_energy=held, &
        energy_correction=energy_correction,energy_policy=policy,wall_temperature=wall_temperature, &
        isothermal_slip_faces=isothermal_slip_faces)
#endif
    endif
    if(.not.ok) error stop 'Reference-state MOL step rejected'
    boundary_integral=boundary_integral+budget;source_integral=source_integral+source
  enddo
#ifdef PELEF_REFERENCE_MPI
  call MPI_Finalize(ierr)
contains
  subroutine collect(local,global)
    type(amr_tree_state_3d),intent(in)::local
    type(amr_tree_state_3d),intent(inout)::global
    integer::peer,status
    do peer=1,2
      if(rank==owners(peer)) global%fields(peer)%values=local%fields(peer)%values
      call MPI_Bcast(global%fields(peer)%values,size(global%fields(peer)%values), &
        MPI_DOUBLE_PRECISION,owners(peer),MPI_COMM_WORLD,status)
      if(status/=MPI_SUCCESS) error stop 'Reference output gather failed'
    enddo
  end subroutine collect
#endif
end program amr_tree_reacting_probe
