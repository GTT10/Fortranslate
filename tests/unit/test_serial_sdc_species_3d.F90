program test_serial_sdc_species_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    advance_reactive_eb_tree_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(0)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: geometry,extended(1)
  type(eb_reactive_tree_3d) :: state,saved
  real(dp), allocatable :: q(:),u(:),root(:,:,:,:),temperature(:,:,:),b(:),w(:),ch(:),force(:),aux(:)
  real(dp) :: t,c,theta,source_sign
  integer :: nv,s,i,j,k,absent,calls,mode
  logical :: ok
  character(len=16) :: faces(6)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok);call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(q(reactive_nprim(size(species))),u(nv),root(nv,2,2,2),temperature(2,2,2), &
    b(nv),w(nv),ch(nv),force(nv),aux(nv))
  q=0;q(1)=1;q(5)=1e5_dp
  do s=1,size(species)
    if(species(s)%name=='N2') q(5+s)=1
  end do
  absent=5+findloc(q(6:),0._dp,dim=1)
  call reactive_primitive_to_conserved(species,q,u,t,c,ok);call require(ok,'state')
  do k=1,2
  do j=1,2
  do i=1,2
    root(:,i,j,k)=u
  end do
  end do
  end do
  temperature=t
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],boxes,topology,ok);call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(2,2,2,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    'x',-10._dp,geometry,ok);call require(ok,'geometry')
  call build_axis_plane_eb_geometry_3d(14,14,14,-3._dp,4._dp,-3._dp,4._dp,-3._dp,4._dp, &
    'x',-10._dp,extended(1),ok);call require(ok,'extended geometry')
  call initialize_reactive_eb_tree_3d(species,topology,geometry,root,temperature,saved,ok)
  call require(ok,'initialization')
  faces='periodic'
  do mode=1,2
    state=saved;b=17;w=19;ch=23;force=29;aux=31;theta=37;calls=0
    source_sign=real(2*mode-3,dp)
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-10_dp,1,.false., &
      [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp, &
      state,extended,faces,b,w,ch,force,theta,ok,auxiliary_source=source,auxiliary_change=aux)
    call require(calls==2,'both source stages reached')
    if(mode==1) then
      call require(.not.ok,'negative accepted species rejected')
      call require(all(state%nodes(1)%state==saved%nodes(1)%state),'state rollback')
      call require(all(state%nodes(1)%temperature==saved%nodes(1)%temperature),'temperature rollback')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(force==29).and.all(aux==31) &
        .and.theta==37,'diagnostic rollback')
    else
      call require(ok,'positive trace species accepted')
      call require(maxval(abs(state%nodes(1)%state(absent,:,:,:)-5e-19_dp))<1e-30_dp,'trace update')
      call require(abs(aux(absent)-5e-19_dp)<1e-30_dp,'trace inventory')
    end if
  end do
  print *, 'PASS: serial SDC species acceptance and rollback'
contains
  subroutine source(state,interval,stage,depth,rate,ok)
    type(eb_reactive_tree_3d), intent(in) :: state
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    calls=calls+1
    rate%fields(1)%values=0
    if(stage>0) then
      rate%fields(1)%values(1,:,:,:)=source_sign*1e-8_dp
      rate%fields(1)%values(absent,:,:,:)=source_sign*1e-8_dp
    end if
    ok=.true.
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(condition) return
    print *, 'FAIL: ',label
    error stop 1
  end subroutine
end program
