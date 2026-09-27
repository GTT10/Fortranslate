program test_les_amr_sdc
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    prolong_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use les_amr_sdc_mod, only: advance_amr_patch_tree_les_sdc_3d
  use les_sdc_3d_mod, only: advance_les_sdc_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction) :: reactions(0)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patch(1)
  type(amr_tree_state_3d) :: u,t,saved,saved_t
  real(dp) :: root(12,4,4,4),rt(1,4,4,4),q(12),sound,h(3),dt,bc(12),sc(12),correction,theta,diffusion
  real(dp) :: initial(12),final(12),scale(12),maximum,serial(12,4,4,4),st(4,4,4),sb(12),ss(12),sd
  integer :: i,j,k,c,geometry,mode,node
  character(len=16) :: faces(6),models(4)
  logical :: ok,flags(4)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  models=[character(len=16)::'smagorinsky','wale','vreman','dynamic']
  do k=1,4
  do j=1,4
  do i=1,4
    q=0;q(1)=.3_dp;q(2)=2*j;q(3)=3*i;q(4)=k;q(5)=101325+10*j;q(12)=1
    call reactive_primitive_to_conserved(species,q,root(:,i,j,k),rt(1,i,j,k),sound,ok)
    call require(ok,'initial state')
  end do
  end do
  end do
  h=.001_dp;dt=1e-9_dp;faces='periodic';maximum=0
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patch(1),ok);call require(ok,'patch')
  do geometry=1,2
    if(geometry==1) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer::],patch(:0),tree,ok)
    else
      call initialize_amr_patch_tree_3d([4,4,4],[1],patch,tree,ok)
    end if
    call require(ok,'tree')
    call prolong_amr_patch_tree_3d(tree,root,saved,ok);call require(ok,'gas')
    call prolong_amr_patch_tree_3d(tree,rt,saved_t,ok);call require(ok,'temperature')
    call composite_integrals_amr_patch_tree_3d(tree,saved,h,initial,ok);call require(ok,'initial inventory')
    do mode=1,4
      flags=.true.
      if(mode==3) flags=.false.
      u=saved;t=saved_t
      call advance_amr_patch_tree_les_sdc_3d(tree,species,reactions,transport,u,t,h,dt,faces,2,.false.,flags, &
        [0._dp,0._dp,0._dp],1e-6_dp,1e-12_dp,.17_dp,.7_dp,.7_dp,bc,sc,correction,theta,diffusion,ok,model=trim(models(mode)))
      call require(ok.and.diffusion>0,'coupled AMR LES SDC')
      call composite_integrals_amr_patch_tree_3d(tree,u,h,final,ok);call require(ok,'final inventory')
      scale=max(abs(initial),abs(final),1e-20_dp);scale(2:4)=max(scale(2:4),sqrt(initial(1)*initial(5)))
      maximum=max(maximum,maxval(abs(final-initial-bc-sc)/scale))
      call require(all(abs(final-initial-bc-sc)<=3e-12_dp*scale),'refluxed composite LES inventory')
      if(geometry==1) then
        serial=root;st=rt(1,:,:,:)
        call advance_les_sdc_3d(species,serial,st,h,dt,0._dp,2,.17_dp,.7_dp,.7_dp,sb,ss,sd,ok, &
          model=trim(models(mode)),transport=transport,transport_flags=flags)
        call require(ok,'regular LES reference')
        call require(maxval(abs(u%fields(1)%values-serial)/max(1._dp,abs(serial)))<1e-12_dp,'root LES SDC parity')
      end if
    end do
    u=saved;t=saved_t;bc=42;sc=43;diffusion=44
    call advance_amr_patch_tree_les_sdc_3d(tree,species,reactions,transport,u,t,h,dt,faces,2,.false.,flags, &
      [0._dp,0._dp,0._dp],1e-6_dp,1e-12_dp,.17_dp,.7_dp,.7_dp,bc,sc,correction,theta,diffusion,ok, &
      model='dynamic',test_ratio=-1)
    call require(.not.ok.and.all(bc==42).and.all(sc==43).and.diffusion==44,'AMR LES failure diagnostics')
    do node=1,size(tree%nodes)
      call require(all(u%fields(node)%values==saved%fields(node)%values).and. &
        all(t%fields(node)%values==saved_t%fields(node)%values),'AMR LES failure state rollback')
    end do
  end do
  print *, 'PASS: common AMR LES SDC four closures, root parity, reflux; residual ',maximum
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine
end program
