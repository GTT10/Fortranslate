program test_eb_native_donor_admission_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d
  use eb_reactive_reflux_3d_mod, only: supply_reactive_eb_gas_eos_donor_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: geometry
  real(dp), allocatable :: target(:),donor(:),debt(:),donor_debt(:),a(:),b(:),packet(:),y(:),pv(:)
  real(dp), allocatable :: saved_target(:),saved_donor(:),saved_debt(:),saved_donor_debt(:)
  real(dp) :: mw,rg,cp,cv,gamma,h,eg,eh,s,t,sound,bias(2),scale_error
  integer :: nv,sp,mode,calls,target_calls,donor_calls
  logical :: ok,complete
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(target(nv+1),donor(nv+1),debt(nv+1),donor_debt(nv+1),a(nv+1),b(nv+1),packet(nv+1), &
    y(size(species)),pv(5+size(species)))
  y=0
  do sp=1,size(species)
    if(species(sp)%name=='N2') y(sp)=1
  enddo
  call require(sum(y)==1,'pure inert composition')
  call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,h,eg,s,ok)
  call require(ok,'gas caloric state')
  call mixture_mass_properties(species,y,700._dp,mw,rg,cp,cv,gamma,h,eh,s,ok)
  call require(ok,'independent held caloric state')
  target=0;target(1)=1;target(5)=eg;target(6:nv)=y;target(nv+1)=eh
  donor=target;debt=.1_dp*target;donor_debt=.05_dp*donor;bias=eh-eg
  call require(target(5)/=target(nv+1),'held must differ from gas energy')
  saved_target=target;saved_donor=donor;saved_debt=debt;saved_donor_debt=donor_debt
  call build_axis_plane_eb_geometry_3d(1,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    'z',-1._dp,geometry,ok)
  call require(ok,'regular geometry')
  call static_valid(target+debt)
  call static_valid(donor-donor_debt)
  do mode=1,9
    a=-17;b=-19;packet=-23;calls=0;target_calls=0;donor_calls=0
    select case(mode)
    case(1)
      call supply_reactive_eb_gas_eos_donor_3d(species,geometry,[1,1,1],1000._dp,geometry,[1,1,1], &
        1000._dp,target,debt,1._dp,donor,donor_debt,2._dp,a,b,packet,complete,ok,reactor_energy_bias=bias)
    case(2)
      call invoke(.false.)
    case(6)
      call supply_reactive_eb_gas_eos_donor_3d(species,geometry,[1,1,1],1000._dp,geometry,[1,1,1], &
        1000._dp,target,debt,1._dp,donor,donor_debt,2._dp,a,b,packet,complete,ok, &
        reactor_energy_bias=bias,native_projection=admit,check_native_target=.true.)
    case(7)
      call supply_reactive_eb_gas_eos_donor_3d(species,geometry,[1,1,1],1000._dp,geometry,[1,1,1], &
        1000._dp,target,debt,1._dp,donor,donor_debt,2._dp,a,b,packet,complete,ok, &
        reactor_energy_bias=bias,check_native_target=.true.)
    case(9)
      call supply_reactive_eb_gas_eos_donor_3d(species,geometry,[1,1,1],1000._dp,geometry,[1,1,1], &
        1000._dp,target,debt,1._dp,donor,donor_debt,2._dp,a,b,packet,complete,ok, &
        reactor_energy_bias=bias,native_projection=admit,inventory_nodes=[11,22])
    case default
      call invoke(.true.)
    end select
    print *, 'NATIVE_DONOR',mode,ok,complete,calls,target_calls,donor_calls
    if(mode==4) then
      call require(ok.and.complete,'static-valid but native-denied target receives material')
      call require(a(1)>=1.25_dp.and.packet(1)>0.and.b(1)<donor(1),'nonzero real donor receipt')
      call require(target_calls>1.and.donor_calls>0,'both native inventories queried')
      scale_error=maxval(abs((a-target)+2*(b-donor))/max(1._dp,abs(target)))
      call require(scale_error<128*epsilon(1._dp),'paired gas/species/held conservation')
      call require(maxval(abs((a-target)-packet)/max(1._dp,abs(target)))<128*epsilon(1._dp), &
        'target receipt equals integrated packet')
      call require(maxval(abs(2*(donor-b)-packet)/max(1._dp,abs(target)))<128*epsilon(1._dp), &
        'donor debit equals integrated packet')
      call require(maxval(abs(packet-packet(1)*donor)/max(1._dp,abs(donor)))<128*epsilon(1._dp), &
        'actual full material vector transports independent held energy')
      call static_valid(a+debt)
      call static_valid(b-donor_debt)
    else
      call require(.not.ok.and..not.complete,'no packet is published for no-op or refusal')
      call require(all(a==-17).and.all(b==-19).and.all(packet==-23),'all refused outputs unchanged')
    endif
    if(mode==1.or.mode==2.or.mode==6.or.mode==7.or.mode==9) &
      call require(calls==0,'legacy/disabled/invalid config no callback')
    if(mode==3) call require(calls==1,'accepted untouched endpoint needs no donor search')
    if(mode==5) call require(calls>1.and.donor_calls==0,'target refusal never publishes or spends donor')
    if(mode==8) call require(donor_calls>0,'donor refusal is exercised')
    call require(all(target==saved_target).and.all(donor==saved_donor).and. &
      all(debt==saved_debt).and.all(donor_debt==saved_donor_debt),'input vectors remain exact')
  enddo
  print *, 'PASS: opt-in native target selection, exact pre-debt query, paired material receipt and rollback'
contains
  subroutine invoke(enabled)
    logical, intent(in) :: enabled
    call supply_reactive_eb_gas_eos_donor_3d(species,geometry,[1,1,1],1000._dp,geometry,[1,1,1], &
      1000._dp,target,debt,1._dp,donor,donor_debt,2._dp,a,b,packet,complete,ok,reactor_energy_bias=bias, &
      native_projection=admit,inventory_nodes=[11,22],check_native_target=enabled)
  end subroutine invoke
  subroutine admit(node,cell,pre,require_post,good)
    integer, intent(in) :: node,cell(3)
    real(dp), intent(in) :: pre(:)
    logical, intent(in) :: require_post
    logical, intent(out) :: good
    calls=calls+1
    call require(all(cell==1).and.require_post,'exact cell and post-admission request')
    call require(size(pre)==nv+1,'complete gas and held packet')
    if(node==11) then
      target_calls=target_calls+1
      if(target_calls==1) call require(all(pre==target),'initial query is PRE debt, not static endpoint')
      good=mode==3.or.(pre(1)>=1.25_dp.and.mode/=5)
    else
      call require(node==22,'donor node identity')
      donor_calls=donor_calls+1
      good=mode/=8
    endif
  end subroutine admit
  subroutine static_valid(value)
    real(dp), intent(in) :: value(:)
    logical :: valid
    call recover_reactive_eb_cell_3d(species,geometry,[1,1,1],value(:nv),1000._dp,pv,t,sound,valid)
    call require(valid,'gas endpoint statically admissible')
    call recover_reactive_eb_cell_3d(species,geometry,[1,1,1],value(:nv),700._dp,pv,t,sound,valid, &
      internal_energy_density=value(nv+1))
    call require(valid,'independent held endpoint statically admissible')
  end subroutine static_valid
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    endif
  end subroutine require
end program test_eb_native_donor_admission_3d
