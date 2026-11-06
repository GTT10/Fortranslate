! Matched-state probe for the frozen two-level periodic HIT control.
! This exercises the production tree integrator with imported SI cell states.
! No chemistry, transport, initialization, projection or clipping is applied.
program amr_tree_mol_probe
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_gamma_law_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_sdc_3d,advance_amr_patch_tree_sdc_subcycle_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patch(1)
  type(amr_patch_tree_3d) :: tree
  type(amr_tree_state_3d) :: state,temperature,seen
  real(dp) :: spacing(3),dt,budget(6),source(6),correction,theta,values(7),force(3)
  integer :: unit,ios,steps,iterations,subcycle,node,i,j,k,records,row,n(3),step
  logical :: ok
  character(len=1024) :: input,prefix,path
  character(len=32) :: magic,faces(6)
  if(command_argument_count()/=2) error stop 'Usage: probe input.dat output-prefix'
  call get_command_argument(1,input);call get_command_argument(2,prefix)
  open(newunit=unit,file=trim(input),status='old',action='read',iostat=ios)
  if(ios/=0) error stop 'Cannot open probe input'
  read(unit,*,iostat=ios) magic
  if(ios/=0.or.magic/='PELEF_AMR_MOL_REFERENCE_V1') error stop 'Invalid probe input magic'
  read(unit,*,iostat=ios) steps,iterations,subcycle,dt
  if(ios/=0) error stop 'Invalid stepping row'
  if(steps<0.or.steps>100.or.iterations<1.or.(subcycle/=0.and.subcycle/=1).or.dt<=0) error stop 'Invalid stepping controls'
  read(unit,*,iostat=ios) spacing
  if(ios/=0.or.any(spacing<=0)) error stop 'Invalid probe spacing'
  read(unit,*,iostat=ios) records
  if(ios/=0.or.records/=320) error stop 'Probe requires all 64 coarse and 256 fine cells'
  call load_gamma_law_thermo(1.4_dp,species,ok)
  if(.not.ok) error stop 'Cannot load gamma-law thermo'
  allocate(reactions(0),transport(0))
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok)
  if(.not.ok) error stop 'Invalid reference patch'
  call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
  if(.not.ok) error stop 'Invalid reference tree'
  allocate(state%fields(2),temperature%fields(2),seen%fields(2))
  do node=1,2
    n=tree%nodes(node)%shape
    allocate(state%fields(node)%values(6,n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)), &
      seen%fields(node)%values(1,n(1),n(2),n(3)))
    state%fields(node)%values=0;temperature%fields(node)%values=0;seen%fields(node)%values=0
  enddo
  do row=1,records
    read(unit,*,iostat=ios) node,i,j,k,values
    if(ios/=0.or.node<1.or.node>2) error stop 'Invalid imported cell'
    n=tree%nodes(node)%shape
    if(any([i,j,k]<1).or.any([i,j,k]>n)) error stop 'Imported cell index outside patch'
    if(seen%fields(node)%values(1,i,j,k)/=0) error stop 'Duplicate imported cell'
    seen%fields(node)%values(1,i,j,k)=1
    state%fields(node)%values(:,i,j,k)=values(:6)
    temperature%fields(node)%values(1,i,j,k)=values(7)
  enddo
  read(unit,*,iostat=ios) magic
  if(ios>=0) error stop 'Unexpected records after imported cells'
  close(unit)
  do node=1,2
    if(any(seen%fields(node)%values/=1)) error stop 'Missing imported cell'
  enddo
  faces='periodic';force=0
  do step=0,steps
    write(path,'(a,"_",i3.3,".csv")') trim(prefix),step
    open(newunit=unit,file=trim(path),status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'Cannot open probe output'
    write(unit,'(a)') 'node,level,i,j,k,time_s,density,xmom,ymom,zmom,rho_E,rho_AIR,Temp'
    do node=1,2
      n=tree%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        write(unit,'(5(i0,","),7(es25.16e3,","),es25.16e3)') node,tree%nodes(node)%level,i,j,k,step*dt, &
          state%fields(node)%values(:,i,j,k),temperature%fields(node)%values(1,i,j,k)
      enddo
      enddo
      enddo
    enddo
    close(unit)
    if(step==steps) exit
    if(subcycle==1) then
      call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,.false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp, &
        budget,source,correction,theta,ok,mol_mode=.true.)
    else
      call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,.false.,[.false.,.false.,.false.,.false.],force,1e-6_dp,1e-14_dp, &
        budget,source,correction,theta,ok,mol_mode=.true.)
    endif
    if(.not.ok) error stop 'Reference-state MOL step rejected'
  enddo
end program amr_tree_mol_probe
