program test_eb_reaction_admission_contract_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use h2o2_full_mechanism_mod, only: load_h2o2_full_mechanism
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval, &
    reaction_admission_not_attempted,reaction_admission_accepted,reaction_admission_rejected, &
    reaction_admission_unknown,reaction_admission_provider_refused
  use amr_hierarchy_3d_mod, only: amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,advance_reactive_eb_tree_sdc_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use dual_energy_mod, only: dual_energy_policy
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: patches(0)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,padded,geometry(1)
  type(eb_geometry_3d), allocatable :: extended(:)
  type(eb_reactive_tree_3d) :: initial,work,authority
  type(amr_tree_state_3d) :: initial_h,held,authority_h
  type(eb_integrated_flux_3d), allocatable :: faces_out(:),scalar_out(:),saved_faces(:),saved_scalar(:)
  type(dual_energy_policy) :: policy
  real(dp) :: root(15,4,4,4),root_t(4,4,4),phi(0:4,0:4,0:4),padded_phi(0:16,0:16,0:16)
  real(dp) :: q(15),cell(15),t,sound,r,b(15),w(15),ch(15),ext(15),theta,prefix,pressure,projection,closure
  real(dp) :: authority_b(15),authority_w(15),authority_ch(15),authority_ext(15),authority_diagnostics(5)
  real(dp) :: last_packet(26)
  integer :: mode,actual_calls,admission_calls,legacy_calls,status,c,last_provider_preview
  logical :: ok,refused,refuse_actual,with_admission
  character(len=16) :: boundaries(6)='outflow'
  call load_h2o2_full_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_full_mechanism(reactions,ok);call require(ok,'mechanism')
  call load_h2o2_full_transport(transport,ok);call require(ok,'transport')
  q=0;q(1)=1;q(15)=1
  r=mixture_specific_gas_constant(species,q(6:),ok);call require(ok,'gas constant')
  q(5)=r*1000
  call reactive_primitive_to_conserved(species,q,cell,t,sound,ok);call require(ok,'uniform inactive gas')
  do c=1,15
    root(c,:,:,:)=cell(c)
  enddo
  root_t=t;phi=1;padded_phi=1
  call initialize_amr_patch_tree_3d([4,4,4],[integer::],patches,topology,ok)
  call require(ok,'root-only topology')
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],finest,ok)
  call require(ok,'full-fluid geometry')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,root_t,initial,ok,internal_energy=initial_h)
  call require(ok,'initial gas and held field')
  geometry(1)=initial%nodes(1)%geometry
  call build_level_set_eb_geometry_3d(padded_phi,[-1.5_dp,-1.5_dp,-1.5_dp],[2.5_dp,2.5_dp,2.5_dp],padded,ok)
  call require(ok,'six-cell padded geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,ok,halo_width=6)
  call require(ok,'extended geometry')
  policy=dual_energy_policy();policy%eta2=1;policy%update_total_from_internal=.false.
  do mode=1,9
    work=initial;held=initial_h
    actual_calls=0;admission_calls=0;last_provider_preview=0;last_packet=0
    b=17;w=19;ch=23;ext=29;theta=.75_dp;prefix=7;pressure=11;projection=13;closure=17
    status=999;refused=.true.;refuse_actual=mode==7.or.mode==8
    with_admission=mode/=1.and.mode/=8
    if(mode>=2) then
      faces_out=saved_faces;scalar_out=saved_scalar
    endif
    if(mode==9) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-8_dp,1,.true., &
        [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-8_dp,1.e-12_dp,work,extended,boundaries, &
        b,w,ch,ext,theta,ok,internal_energy=held,internal_energy_flux=scalar_out,integrated_flux=faces_out, &
        energy_policy=policy,energy_correction=prefix,internal_energy_pressure_work_change=pressure, &
        internal_energy_projection_change=projection,internal_energy_closure_change=closure, &
        evaluate_reaction_admission=admission,reaction_admission_status=status,provider_refused=refused)
    else if(with_admission) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-8_dp,1,.true., &
        [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-8_dp,1.e-12_dp,work,extended,boundaries, &
        b,w,ch,ext,theta,ok,internal_energy=held,internal_energy_flux=scalar_out,integrated_flux=faces_out, &
        energy_policy=policy,energy_correction=prefix,internal_energy_pressure_work_change=pressure, &
        internal_energy_projection_change=projection,internal_energy_closure_change=closure, &
        evaluate_reaction=provider,evaluate_reaction_admission=admission,reaction_admission_status=status, &
        provider_refused=refused)
    else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1.e-8_dp,1,.true., &
        [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1.e-8_dp,1.e-12_dp,work,extended,boundaries, &
        b,w,ch,ext,theta,ok,internal_energy=held,internal_energy_flux=scalar_out,integrated_flux=faces_out, &
        energy_policy=policy,energy_correction=prefix,internal_energy_pressure_work_change=pressure, &
        internal_energy_projection_change=projection,internal_energy_closure_change=closure, &
        evaluate_reaction=provider,reaction_admission_status=status,provider_refused=refused)
    endif
    print *, 'SERIAL_ADMISSION',mode,ok,status,actual_calls,admission_calls,refused
    if(mode<=2) then
      call require(ok.and..not.refused,'ordinary inactive interval accepted')
      if(mode==1) then
        call require(status==reaction_admission_not_attempted.and.admission_calls==0,'legacy capability absent')
        call require(actual_calls==64,'legacy ordinary provider called once per root cell')
        legacy_calls=actual_calls;authority=work;authority_h=held
        authority_b=b;authority_w=w;authority_ch=ch;authority_ext=ext
        authority_diagnostics=[theta,prefix,pressure,projection,closure]
        saved_faces=faces_out;saved_scalar=scalar_out
      else
        call require(status==reaction_admission_accepted,'successful final status accepted')
        call require(actual_calls==legacy_calls.and.admission_calls>=actual_calls, &
          'capability never duplicates evaluator side effects')
        call require(all(work%nodes(1)%state==authority%nodes(1)%state).and. &
          all(work%nodes(1)%temperature==authority%nodes(1)%temperature).and. &
          all(held%fields(1)%values==authority_h%fields(1)%values),'legacy and opt-in accepted state exact')
        call require(all(b==authority_b).and.all(w==authority_w).and.all(ch==authority_ch).and. &
          all(ext==authority_ext).and.all([theta,prefix,pressure,projection,closure]==authority_diagnostics), &
          'legacy and opt-in accepted ledgers exact')
        call check_faces()
      endif
    else
      call require(.not.ok,'nonacceptance refuses full caller interval')
      if(mode<=6.or.mode==9) call require(actual_calls==0,'admission failure never calls stateful evaluator')
      if(mode==3.or.mode==5.or.mode==9) call require(status==reaction_admission_unknown.and..not.refused, &
        'unknown is not provider refusal')
      if(mode==4) call require(status==reaction_admission_rejected.and..not.refused,'explicit certificate status retained')
      if(mode==6.or.mode==7) call require(status==reaction_admission_provider_refused.and.refused,'terminal refusal retained')
      if(mode==7.or.mode==8) call require(actual_calls==1.and.refused,'actual provider refusal is never retried')
      if(mode==8) call require(status==reaction_admission_not_attempted,'legacy refusal does not invent an admission result')
      call require(all(work%nodes(1)%state==initial%nodes(1)%state).and. &
        all(work%nodes(1)%temperature==initial%nodes(1)%temperature).and. &
        all(held%fields(1)%values==initial_h%fields(1)%values),'failure preserves all gas temperature and held values')
      call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(ext==29).and.theta==.75_dp.and. &
        prefix==7.and.pressure==11.and.projection==13.and.closure==17,'failure preserves every supplied ledger')
      call check_faces()
    endif
  enddo
  print *, 'PASS: serial provider call authority, exact packet, status distinctions and allocated-output rollback'
contains
  subroutine require(condition,label)
    logical,intent(in) :: condition
    character(len=*),intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL ',mode,': ',label
      error stop 1
    endif
  end subroutine
  subroutine check_faces()
    call require(allocated(faces_out).and.allocated(scalar_out),'retained allocated flux outputs')
    call require(all(faces_out(1)%x==saved_faces(1)%x).and.all(faces_out(1)%y==saved_faces(1)%y).and. &
      all(faces_out(1)%z==saved_faces(1)%z).and.all(scalar_out(1)%x==saved_scalar(1)%x).and. &
      all(scalar_out(1)%y==saved_scalar(1)%y).and.all(scalar_out(1)%z==saved_scalar(1)%z),'gas and scalar flux outputs exact')
  end subroutine
  subroutine admission(table,rates,energy,source_species,source_energy,duration,rtol,atol,density,temp,result)
    type(nasa7_species),intent(in) :: table(:)
    type(elementary_reaction),intent(in) :: rates(:)
    real(dp),intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol,density(:),temp
    integer,intent(out) :: result
    admission_calls=admission_calls+1
    last_packet=[energy,source_energy,duration,rtol,atol,temp,density,source_species]
    result=reaction_admission_accepted
    select case(mode)
    case(3);result=reaction_admission_unknown
    case(4);result=reaction_admission_rejected ! Manufactured status plumbing, not a physical certificate.
    case(5);result=999
    case(6);result=reaction_admission_provider_refused
    end select
  end subroutine
  subroutine provider(table,rates,energy,source_species,source_energy,duration,rtol,atol,density,temp,steps,good,maximum_steps)
    type(nasa7_species),intent(in) :: table(:)
    type(elementary_reaction),intent(in) :: rates(:)
    real(dp),intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
    real(dp),intent(inout) :: density(:),temp
    integer,intent(out) :: steps
    logical,intent(out) :: good
    integer,optional,intent(in) :: maximum_steps
    actual_calls=actual_calls+1
    if(with_admission) then
      call require(admission_calls>last_provider_preview,'every actual provider call has a fresh final admission')
      last_provider_preview=admission_calls
      call require(all(last_packet==[energy,source_energy,duration,rtol,atol,temp,density,source_species]), &
        'immediately certified full packet equals actual evaluator input')
    endif
    if(refuse_actual) then
      density=-1;temp=-1;steps=0;good=.false.;return
    endif
    call advance_forced_reactor_interval(table,rates,energy,source_species,source_energy,duration,rtol,atol, &
      density,temp,steps,good,maximum_steps=maximum_steps)
  end subroutine
end program test_eb_reaction_admission_contract_3d
