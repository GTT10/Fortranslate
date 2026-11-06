program test_mol_riemann_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties,nasa7_transient_mass_properties
  use mixture_thermo_mod, only: mixture_density
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved,reconstructed_pelec_state
  use mol_reconstruction_3d_mod, only: mol_characteristic_slope_3d,limit_mol_positive_slope_3d,mol_reconstruct_face_3d
  use mol_riemann_3d_mod, only: mol_riemann_flux_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp), allocatable :: q(:,:),center(:),jump(:),face(:),state(:),flux(:),expected(:),other_state(:)
  real(dp) :: t,c,sound(2),nan,diffusion_speed,potential(2),plain_energy,aux,expected_aux,rhoe
  integer :: status,lax
  character(len=32) :: mode
  logical :: ok,connected(2,2)
  integer :: axis,order,i,checks
  checks=0
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  allocate(q(reactive_nprim(size(species)),-1:2))
  allocate(center(size(q,1)),jump(size(q,1)),face(size(q,1)))
  allocate(state(reactive_nvar(size(species))))
  allocate(flux(size(state)),expected(size(state)),other_state(size(state)))
  if(command_argument_count()>0) then
    call get_command_argument(1,mode)
    select case(trim(mode))
    case('--thermo-stream')
      do i=1,size(species)
        write(*,'(16(es26.17e3,1x))') species(i)%molecular_weight,species(i)%temperature_mid, &
          species(i)%low_coefficients,species(i)%high_coefficients
      end do
    case('--rotating-stream','--internal-energy-stream')
      do
        read(*,*,iostat=status) q(:,0),q(:,1),sound,potential,lax
        if(status<0) exit
        if(status/=0) error stop 'Invalid rotating MOL input'
        q(:,-1)=q(:,0);q(:,2)=q(:,1);connected=.true.
        call mol_riemann_flux_3d(species,q,sound,1,connected,1,flux,ok, &
          lax_friedrichs=lax==1,rotating_potential=potential,internal_energy_flux=aux)
        if(.not.ok) error stop 'Rejected rotating MOL input'
        if(trim(mode)=='--internal-energy-stream') then
          write(*,'(*(es26.17e3,1x))') flux,aux
        else
          write(*,'(*(es26.17e3,1x))') flux
        endif
      end do
    case default
      error stop 'Unknown MOL test mode'
    end select
    stop
  end if
  center=0;center(1:5)=[2._dp,10._dp,-20._dp,30._dp,200000._dp]
  do i=1,size(species)
    if(species(i)%name=='N2') center(5+i)=.7_dp
    if(species(i)%name=='O2') center(5+i)=.3_dp
  end do
  nan=ieee_value(0._dp,ieee_quiet_nan)
  connected=.true.
  call check_pressure_expansion()
  call check_reconstructed_lf_expansion()
  call check_transient_mixture_heat_capacity()
  do axis=1,3
    do order=1,2
      jump=0
      if(order==2) then
        jump(1)=.02_dp;jump(2:4)=[.2_dp,-.3_dp,.4_dp];jump(5)=100
      end if
      do i=-1,2
        q(:,i)=center+real(i,dp)*jump
      end do
      face=center+.5_dp*jump
      call reactive_primitive_to_conserved(species,face,state,t,c,ok)
      call require(ok,'face EOS')
      rhoe=state(5)-.5_dp*state(1)*sum(face(2:4)**2)
      expected_aux=rhoe*face(axis+1)
      sound=c
      expected=state*face(axis+1)
      expected(axis+1)=expected(axis+1)+face(5)
      expected(5)=(state(5)+face(5))*face(axis+1)
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok,internal_energy_flux=aux)
      call require(ok,'MOL flux')
      call require(abs(aux-expected_aux)<2e-13_dp*max(1._dp,abs(expected_aux)),'internal-energy advective flux')
      call require(maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'analytic physical face flux')
      call require(abs(sum(flux(6:))-flux(1))<1e-12_dp,'species flux closes mass')
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok,lax_friedrichs=.true.,internal_energy_flux=aux)
      call require(ok,'MOL LF flux')
      call require(abs(aux-expected_aux)<2e-13_dp*max(1._dp,abs(expected_aux)),'LF internal-energy advective flux')
      call require(maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'LF analytic physical face flux')
      potential=[400._dp,1600._dp];plain_energy=expected(5)
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok,rotating_potential=potential)
      expected(5)=plain_energy-expected(1)*merge(potential(1),potential(2),face(axis+1)>0)
      call require(ok.and.maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'rotating acoustic energy flux')
      expected(5)=plain_energy-.5_dp*expected(1)*sum(potential)+ &
        .5_dp*(abs(face(axis+1))+c)*face(1)*(potential(2)-potential(1))
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok, &
        lax_friedrichs=.true.,rotating_potential=potential)
      call require(ok.and.maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'rotating LF energy diffusion')
      flux=19
      call mol_riemann_flux_3d(species,q,sound,axis,connected,order,flux,ok,rotating_potential=[-1._dp,0._dp])
      call require(.not.ok.and.all(flux==19),'invalid potential rollback')

    end do
    ! Stationary equal-pressure contact: the reference contact speed is zero,
    ! and the LF diffusion speed is the arithmetic mean of EOS sound speeds.
    face=center;face(2:4)=0
    q(:,0)=face
    call reactive_primitive_to_conserved(species,face,state,t,sound(1),ok)
    call require(ok,'LF left stationary contact')
    face(1)=1.5_dp; q(:,1)=face
    call reactive_primitive_to_conserved(species,face,other_state,t,sound(2),ok)
    call require(ok,'LF right stationary contact')
    expected=-.25_dp*sum(sound)*(other_state-state)
    expected(axis+1)=expected(axis+1)+face(5)
    q(:,-1)=nan;q(:,2)=nan
    call mol_riemann_flux_3d(species,q,sound,axis,connected,1,flux,ok,lax_friedrichs=.true.)
    call require(ok.and.maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'LF stationary contact diffusion')
    call require(abs(sum(flux(6:))-flux(1))<1e-12_dp,'LF species mass closure')
    ! Equal density and zero velocities isolate pressure-driven contact speed.
    face=q(:,0);face(5)=300000; q(:,1)=face
    call reactive_primitive_to_conserved(species,face,other_state,t,sound(2),ok)
    call require(ok,'LF pressure jump EOS')
    diffusion_speed=100000._dp/(face(1)*sum(sound))+.5_dp*sum(sound)
    expected=-.5_dp*diffusion_speed*(other_state-state)
    expected(axis+1)=250000
    call mol_riemann_flux_3d(species,q,sound,axis,connected,1,flux,ok,lax_friedrichs=.true.)
    call require(ok.and.maxval(abs(flux-expected)/max(1._dp,abs(expected)))<2e-13_dp,'LF pressure-driven speed')
    q(:,-1)=nan;q(:,2)=nan;q(:,0)=center;q(:,1)=center
    connected(:,1)=[.false.,.true.];connected(:,2)=[.true.,.false.]
    call mol_riemann_flux_3d(species,q,sound,axis,connected,2,flux,ok)
    call require(ok,'inactive outer donors ignored')
    connected=.true.;flux=19
    call mol_riemann_flux_3d(species,q,sound,axis,connected,2,flux,ok)
    call require(.not.ok.and.all(flux==19),'active NaN rollback')
    call mol_riemann_flux_3d(species,q,sound,axis,connected,2,flux,ok,lax_friedrichs=.true.)
    call require(.not.ok.and.all(flux==19),'LF active NaN rollback')
  end do
  print *, 'PASS MOL Riemann checks:',checks
contains
  subroutine check_transient_mixture_heat_capacity()
    real(dp) :: cp,cv,h,e,s,r,transient(size(center)),u(size(flux)),cc,tt,observed(size(flux))
    integer :: m,a
    m=0
    do a=1,size(species)
      if(species(a)%name=='O2') m=a
    end do
    call require(m>0,'transient O2 species')
    call nasa7_transient_mass_properties(species(m),7000._dp,cp,cv,h,e,s,ok)
    call require(ok.and.cv<0,'transient signed species heat capacity')
    r=sum(center(6:)*[(8314.46261815324_dp/species(a)%molecular_weight,a=1,size(species))])
    transient=center;transient(5)=transient(1)*r*7000._dp
    call reconstructed_pelec_state(species,transient,u,tt,cc,ok)
    call require(ok.and.cc>0,'transient positive mixture heat capacity')
    do a=1,3
      q=spread(transient,2,4);sound=cc
      expected=u*transient(a+1);expected(a+1)=expected(a+1)+transient(5)
      expected(5)=(u(5)+transient(5))*transient(a+1)
      call mol_riemann_flux_3d(species,q,sound,a,connected,1,observed,ok,lax_friedrichs=.true.)
      call require(ok.and.maxval(abs(observed-expected)/max(1._dp,abs(expected)))<2e-13_dp, &
        'transient LF mixture flux with signed species Cv')
    end do
    transient(5)=transient(1)*r*8000._dp
    call reconstructed_pelec_state(species,transient,u,tt,cc,ok)
    call require(.not.ok,'transient nonpositive mixture heat capacity rejects')
    q=spread(transient,2,4);observed=19
    call mol_riemann_flux_3d(species,q,sound,1,connected,1,observed,ok,lax_friedrichs=.true.)
    call require(.not.ok.and.all(observed==19),'transient invalid mixture flux rollback')
  end subroutine check_transient_mixture_heat_capacity

  subroutine check_reconstructed_lf_expansion()
    real(dp) :: pressures(2),density(4),temperatures(4),r,slopes(size(species)+4,2),l(size(center)),rr(size(center))
    real(dp) :: ul(size(flux)),ur(size(flux)),cl,cr,tl,tr,contact,speed,observed(size(flux)),saved(size(q,1),4)
    integer :: a,m
    density=[.25_dp,.24_dp,1.4_dp,.68_dp];temperatures=[900._dp,301._dp,301._dp,1500._dp]
    r=sum(center(6:)*[(8314.46261815324_dp/species(m)%molecular_weight,m=1,size(species))])
    do a=1,3
      do m=-1,2
        q(:,m)=center;q(1,m)=density(m+2);q(2:4,m)=0
        q(5,m)=density(m+2)*r*temperatures(m+2)
        call reactive_primitive_to_conserved(species,q(:,m),state,t,c,ok)
        call require(ok,'reconstructed LF valid strict cell EOS')
        if(m==0.or.m==1) sound(m+1)=c
      end do
      saved=q
      do m=0,1
        call mol_characteristic_slope_3d(q(:,m-1),q(:,m),q(:,m+1),sound(m+1),a, &
          connected(:,m+1),slopes(:,m+1),ok)
        call require(ok,'reconstructed LF characteristic slope')
        call limit_mol_positive_slope_3d(q(:,m),sound(m+1),a,slopes(:,m+1),ok)
        call require(ok,'reconstructed LF positive slope')
      end do
      call mol_reconstruct_face_3d(q(:,0),sound(1),a,1,slopes(:,1),l,ok)
      call require(ok,'reconstructed LF left face')
      call mol_reconstruct_face_3d(q(:,1),sound(2),a,-1,slopes(:,2),rr,ok)
      call require(ok,'reconstructed LF right face')
      call reactive_primitive_to_conserved(species,l,state,t,c,ok)
      call require(.not.ok,'reconstructed LF left face outside strict fit range')
      call reconstructed_pelec_state(species,l,ul,tl,cl,ok)
      call require(ok.and.tl>0.and.tl<300,'reconstructed LF transient left EOS')
      call reconstructed_pelec_state(species,rr,ur,tr,cr,ok)
      call require(ok,'reconstructed LF transient right EOS')
      contact=(l(5)-rr(5))/(l(1)*cl+rr(1)*cr)
      speed=abs(contact)+.5_dp*(cl+cr)
      expected=-.5_dp*speed*(ur-ul)
      expected(a+1)=expected(a+1)+.5_dp*(l(5)+rr(5))
      call mol_riemann_flux_3d(species,q,sound,a,connected,2,observed,ok, &
        lax_friedrichs=.true.,reconstructed_pressure=pressures)
      call require(all(pressures==[l(5),rr(5)]),"diagnostic actual reconstructed pressures")
      call require(ok,'LF flux accepts positive transient reconstructed faces')
      call require(maxval(abs(observed-expected)/max(1._dp,abs(expected)))<2e-13_dp, &
        'reconstructed LF stationary diffusive flux')
      potential=[400._dp,1600._dp]
      expected(5)=expected(5)+.5_dp*speed*(ur(1)*potential(2)-ul(1)*potential(1))
      call mol_riemann_flux_3d(species,q,sound,a,connected,2,observed,ok, &
        lax_friedrichs=.true.,rotating_potential=potential)
      call require(ok.and.maxval(abs(observed-expected)/max(1._dp,abs(expected)))<2e-13_dp, &
        'reconstructed LF rotating transient energy diffusion')
      call require(all(q==saved),'reconstructed LF immutable original cells')
      pressures=-17
      call mol_riemann_flux_3d(species,q,sound,0,connected,2,observed,ok, &
        reconstructed_pressure=pressures)
      call require(.not.ok.and.all(pressures==-17),'diagnostic pressures rollback')
    end do
  end subroutine check_reconstructed_lf_expansion

  subroutine check_pressure_expansion()
    real(dp) :: left(size(center)),right(size(center)),mirror(size(flux)),observed(size(flux))
    real(dp) :: cp,cv,h,e,s,density,ratio
    integer :: a,p
    ! Strict accepted-cell properties keep the NASA fit range. The separate
    ! transient evaluator supplies the star state below that range.
    call nasa7_mass_properties(species(size(species)),20._dp,cp,cv,h,e,s,ok)
    call require(.not.ok,'strict NASA rejects below-range temperature')
    call nasa7_transient_mass_properties(species(size(species)),20._dp,cp,cv,h,e,s,ok)
    call require(ok.and.cv>0,'transient NASA evaluates positive below-range temperature')
    call nasa7_transient_mass_properties(species(size(species)),0._dp,cp,cv,h,e,s,ok)
    call require(.not.ok,'transient NASA rejects zero temperature')
    left=center;left(2:4)=0;left(5)=left(1)*20._dp* &
      sum(center(6:)*[(8314.46261815324_dp/species(p)%molecular_weight,p=1,size(species))])
    call reactive_primitive_to_conserved(species,left,state,t,c,ok)
    call require(.not.ok,'accepted-cell conversion retains strict NASA range')
    do a=1,3
    do p=1,3
      ratio=10._dp**p
      left=center;left(2:4)=0;left(5)=1e6_dp
      density=mixture_density(species,left(6:),left(5),300._dp,ok)
      call require(ok,'pressure expansion valid left cell')
      left(1)=density;right=left;right(1)=density/ratio;right(5)=left(5)/ratio
      call reactive_primitive_to_conserved(species,left,state,t,sound(1),ok)
      call require(ok,'pressure expansion strict left EOS')
      call reactive_primitive_to_conserved(species,right,other_state,t,sound(2),ok)
      call require(ok,'pressure expansion strict right EOS')
      q(:,0)=left;q(:,1)=right;q(:,-1)=left;q(:,2)=right
      call mol_riemann_flux_3d(species,q,sound,a,connected,1,observed,ok)
      call require(ok,'strong-pressure MOL acoustic expansion')
      call require(abs(sum(observed(6:))-observed(1))<1e-12_dp*max(1._dp,abs(observed(1))), &
        'pressure expansion species flux closure')
      q(:,0)=right;q(:,1)=left;q(:,-1)=right;q(:,2)=left
      call mol_riemann_flux_3d(species,q,sound(2:1:-1),a,connected,1,mirror,ok)
      call require(ok,'mirrored strong-pressure MOL expansion')
      expected=-observed;expected(a+1)=observed(a+1)
      call require(maxval(abs(mirror-expected)/max(1._dp,abs(expected)))<2e-13_dp, &
        'pressure expansion mirror symmetry')
    end do
    end do
  end subroutine check_pressure_expansion

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine
end program
