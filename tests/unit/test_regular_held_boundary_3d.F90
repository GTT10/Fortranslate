program test_regular_held_boundary_3d
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use transport_database_mod, only: load_h2o2_full_transport
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_mass_properties
  use dual_energy_mod, only: dual_energy_policy
  use reactive_mol_3d_mod, only: compute_reactive_mol_3d,advance_reactive_mol_3d
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d, &
    fill_mach_inlet_boundary_ghosts_3d,fill_pressure_outlet_ghosts_3d
  use les_mpi_boundary_mod, only: fill_mpi_les_boundary_ghosts_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  implicit none
  type :: fields
    real(dp), allocatable :: gas(:,:,:,:),temp(:,:,:),held(:,:,:),ghost(:,:,:,:),gt(:,:,:),gh(:,:,:)
    real(dp), allocatable :: rates(:,:,:,:),x(:,:,:,:),y(:,:,:,:),z(:,:,:,:),arate(:,:,:),boundary(:)
    type(eb_integrated_flux_3d) :: aflux
    real(dp) :: prefix
  end type
  type(fields) :: run(2),saved
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(dual_energy_policy) :: policy
  real(dp), allocatable :: composition(:)
  real(dp) :: h(3),dt,mach(6),tin(6),pressure(6),walls(3,6),seen
  integer :: rank,nrank,status,nv,n(3),width,case_id,a,c,i,j,k,calls
  character(len=16) :: faces(6)
  logical :: good,fail_new_halo
  call MPI_Init(status);call MPI_Comm_rank(MPI_COMM_WORLD,rank,status);call MPI_Comm_size(MPI_COMM_WORLD,nrank,status)
  call require(nrank==3,'three ranks')
  call load_h2o2_full_thermo(species,good);call require(good,'thermo')
  call load_h2o2_full_transport(transport,good);call require(good,'transport')
  nv=5+size(species);allocate(composition(size(species)));composition=0
  do c=1,size(species)
    if(species(c)%name=='N2') composition(c)=1
  enddo
  h=.001_dp;dt=2e-8_dp;width=3;n=[2,4,4]
  faces='periodic';mach=-1;mach(1)=.2_dp;tin=600;pressure=0;walls=0
  ! All MPI ranks own a physical x slab. The two bulk-energy choices have the
  ! same actual scalar and temperature. They must share every scalar donor.
  do a=1,2
    call initialize(run(a),n,a==1,rank*2)
    call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,run(a)%gas,run(a)%temp,[6,4,4],faces,width, &
      run(a)%ghost,run(a)%gt,good,internal_energy_density=run(a)%held,internal_energy_ghost=run(a)%gh)
    call require(good,'MPI periodic actualH halo')
  enddo
  call require(all(run(1)%gh==run(2)%gh),'MPI actualH halo independent of bulkenergy')
  call require(all(run(1)%gh(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))==run(1)%held), &
    'MPI actualH valid cells retained')
  do case_id=1,2
    faces(1)='inflow';faces(2)='outflow';pressure=0
    if(case_id==2) pressure(2)=80000
    do a=1,2
      call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,run(a)%gas,run(a)%temp,[6,4,4],faces,width, &
        run(a)%ghost,run(a)%gt,good,inflow_temperature=tin,inflow_mach=mach,outlet_pressure=pressure, &
        internal_energy_density=run(a)%held,internal_energy_ghost=run(a)%gh)
      call require(good,'MPI actualH Mach/pressure boundary')
    enddo
    call require(all(run(1)%gh==run(2)%gh).and.all(run(1)%gt==run(2)%gt),'MPI actualH prescribed boundary temperature')
    if(rank==0) call require(all(run(1)%ghost(:,:width,:,:)==run(2)%ghost(:,:width,:,:)), &
      'MPI Mach state uses actualH donor EOS')
    if(rank==2.and.case_id==2) call require(all(run(1)%ghost(:,width+n(1)+1:,:,:)== &
      run(2)%ghost(:,width+n(1)+1:,:,:)),'MPI pressure state uses actualH donor EOS')
  enddo
  saved=run(1)
  if(rank==2) run(1)%held(1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
  call fill_mpi_les_boundary_ghosts_3d(MPI_COMM_WORLD,species,run(1)%gas,run(1)%temp,[6,4,4],faces,width, &
    run(1)%ghost,run(1)%gt,good,inflow_temperature=tin,inflow_mach=mach,outlet_pressure=pressure, &
    internal_energy_density=run(1)%held,internal_energy_ghost=run(1)%gh)
  call require(.not.good,'MPI malformed actualH refuses collectively')
  call require(all(run(1)%ghost==saved%ghost).and.all(run(1)%gt==saved%gt).and.all(run(1)%gh==saved%gh), &
    'MPI actualH halo transactional rollback')
  run(1)=saved
  ! Serial physical constructor, callback constructor and the native MOL step
  ! run only on rank zero. They use identical nonuniform actualH fields.
  if(rank==0) then
    n=[4,4,4];policy%eta2=2;policy%allow_negative_energy=.false.;policy%update_total_from_internal=.false.
    do a=1,2
      call initialize(run(a),n,a==1,0)
      allocate(run(a)%ghost(nv,n(1)+6,n(2)+6,n(3)+6),run(a)%gt(n(1)+6,n(2)+6,n(3)+6), &
        run(a)%gh(n(1)+6,n(2)+6,n(3)+6))
    enddo
    do case_id=1,3
      faces='periodic';pressure=0
      if(case_id==2) then
        faces(1)='inflow';faces(2)='outflow';pressure(2)=80000
      else if(case_id==3) then
        faces(1:2)='no_slip_wall';walls(:,1)=[0._dp,2._dp,0._dp]
      endif
      do a=1,2
        if(case_id==2) then
          call fill_mach_inlet_boundary_ghosts_3d(species,run(a)%gas,run(a)%temp,faces,mach,run(a)%ghost, &
            run(a)%gt,good,inflow_temperature=tin,ghost_width=3, &
            internal_energy_density=run(a)%held,internal_energy_ghost=run(a)%gh)
          call require(good,'physical actualH Mach halo')
          call fill_pressure_outlet_ghosts_3d(species,n,faces,pressure,run(a)%ghost,run(a)%gt,good,ghost_width=3, &
            internal_energy_ghost=run(a)%gh)
        else
          call fill_reactive_boundary_ghosts_3d(run(a)%gas,run(a)%temp,faces,run(a)%ghost,run(a)%gt,good, &
            ghost_width=3,viscous_wall=.true.,wall_velocity=walls, &
            internal_energy_density=run(a)%held,internal_energy_ghost=run(a)%gh)
        endif
        call require(good,'physical actualH donor')
      enddo
      call require(all(run(1)%gh==run(2)%gh).and.all(run(1)%gt==run(2)%gt), &
        'physical actualH independent of bulkenergy')
    enddo
    faces='periodic';walls=0
    ! Spatial output UEINT rates and fluxes must use H, including an explicitly
    ! nonuniform H field, in both physical and callback boundary modes.
    do case_id=1,2
      do a=1,2
        calls=0;seen=0;fail_new_halo=.false.
        if(case_id==1) then
          call compute_reactive_mol_3d(species,run(a)%gas,run(a)%temp,h,dt,0._dp,run(a)%rates, &
            run(a)%x,run(a)%y,run(a)%z,run(a)%boundary,good,internal_energy_density=run(a)%held, &
            internal_energy_rate=run(a)%arate,internal_energy_integrated_flux=run(a)%aflux)
        else
          call compute_reactive_mol_3d(species,run(a)%gas,run(a)%temp,h,dt,0._dp,run(a)%rates, &
            run(a)%x,run(a)%y,run(a)%z,run(a)%boundary,good,internal_energy_density=run(a)%held, &
            internal_energy_rate=run(a)%arate,internal_energy_integrated_flux=run(a)%aflux,fill_halos=halo)
        endif
        call require(good,'standalone actualH MOL spatial')
      enddo
      call require(all(run(1)%arate==run(2)%arate),'standalone actualH rates independent of Eden')
      call require(all(run(1)%aflux%x==run(2)%aflux%x).and.all(run(1)%aflux%y==run(2)%aflux%y).and. &
        all(run(1)%aflux%z==run(2)%aflux%z),'standalone actualH fluxes independent of Eden')
      call require(maxval(abs(run(1)%arate))>1e-6_dp,'nonuniform actualH spatial operator reached')
    enddo
    do a=1,2
      calls=0;seen=0;fail_new_halo=.false.
      call advance_reactive_mol_3d(species,run(a)%gas,run(a)%temp,h,dt,0._dp,2,.false.,1e-6_dp,1e-14_dp,good, &
        transport=transport,transport_flags=[.true.,.true.,.true.,.false.],fill_halos=halo, &
        internal_energy_density=run(a)%held,energy_correction=run(a)%prefix,energy_policy=policy)
      call require(good.and.calls==4.and.seen>1e-8_dp,'standalone old/new H hydro/transport halo chronology')
    enddo
    call require(all(abs(run(1)%held-run(2)%held)<=128*epsilon(1._dp)*max(1._dp,abs(run(2)%held))), &
      'standalone actualH accepted endpoint')
    call require(all(abs(run(1)%temp-run(2)%temp)<=128*epsilon(1._dp)*max(1._dp,abs(run(2)%temp))), &
      'standalone actualH accepted EOS')
    call require(run(1)%prefix==3.5_dp.and.run(2)%prefix==3.5_dp,'standalone numerical prefix retained')
    saved=run(1);calls=0;seen=0;fail_new_halo=.true.
    call advance_reactive_mol_3d(species,run(1)%gas,run(1)%temp,h,dt,0._dp,2,.false.,1e-6_dp,1e-14_dp,good, &
      transport=transport,transport_flags=[.true.,.true.,.true.,.false.],fill_halos=halo, &
      internal_energy_density=run(1)%held,energy_correction=run(1)%prefix,energy_policy=policy)
    call require(.not.good.and.calls==3,'standalone new hydro halo rejection')
    call require(all(run(1)%gas==saved%gas).and.all(run(1)%temp==saved%temp).and.all(run(1)%held==saved%held).and. &
      run(1)%prefix==saved%prefix,'standalone late actualH rollback')
    print *, 'PASS: actualH periodic/wall/Mach/pressure halos, MPI rollback, standalone MOL stages and late rollback'
  endif
  call MPI_Finalize(status)
contains
  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message
    if(ok) return
    print *, 'FAIL rank',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
  subroutine initialize(f,n,bulk_hot,xoffset)
    type(fields), intent(inout) :: f
    integer, intent(in) :: n(3),xoffset
    logical, intent(in) :: bulk_hot
    real(dp) :: xx(3),v(3),tt,mw,rg,cp,cv,gamma,enthalpy,e,e1000,entropy
    if(allocated(f%gas)) deallocate(f%gas,f%temp,f%held,f%rates,f%x,f%y,f%z,f%arate,f%boundary)
    if(allocated(f%ghost)) deallocate(f%ghost,f%gt,f%gh)
    allocate(f%gas(nv,n(1),n(2),n(3)),f%temp(n(1),n(2),n(3)),f%held(n(1),n(2),n(3)))
    allocate(f%rates(nv,n(1),n(2),n(3)),f%x(nv,0:n(1),n(2),n(3)), &
      f%y(nv,n(1),0:n(2),n(3)),f%z(nv,n(1),n(2),0:n(3)),f%arate(n(1),n(2),n(3)),f%boundary(nv))
    call mixture_mass_properties(species,composition,1000._dp,mw,rg,cp,cv,gamma,enthalpy,e1000,entropy,good)
    call require(good,'bulk initial EOS')
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      xx=real([xoffset+i,j,k],dp)-.5_dp
      tt=700+25*sin(2*acos(-1._dp)*xx(1)/6)
      v=[8*sin(2*acos(-1._dp)*xx(1)/6),2*sin(2*acos(-1._dp)*xx(2)/4),0._dp]
      call mixture_mass_properties(species,composition,tt,mw,rg,cp,cv,gamma,enthalpy,e,entropy,good)
      call require(good,'actual initial EOS')
      f%gas(:,i,j,k)=0;f%gas(1,i,j,k)=.4_dp;f%gas(2:4,i,j,k)=.4_dp*v
      f%gas(5,i,j,k)=.4_dp*merge(e1000,e,bulk_hot)+.2_dp*sum(v*v);f%gas(6:,i,j,k)=.4_dp*composition
      f%temp(i,j,k)=tt;f%held(i,j,k)=.4_dp*e
    enddo
    enddo
    enddo
    f%rates=0;f%x=0;f%y=0;f%z=0;f%arate=0;f%boundary=0;f%prefix=3.5_dp
  end subroutine
  subroutine halo(u,t,stamp,width,viscous,ghost,gt,ok,internal_energy_density,internal_energy_ghost)
    real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),stamp
    integer, intent(in) :: width
    logical, intent(in) :: viscous
    real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    real(dp), optional, intent(out) :: internal_energy_ghost(:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable, save :: first_h(:,:,:)
    calls=calls+1
    call require(present(internal_energy_density).and.present(internal_energy_ghost),'native actualH halo mode')
    if(calls==1) first_h=internal_energy_density
    seen=max(seen,maxval(abs(internal_energy_density-first_h)))
    call fill_reactive_boundary_ghosts_3d(u,t,[character(len=8)::'periodic','periodic','periodic', &
      'periodic','periodic','periodic'],ghost,gt,ok,ghost_width=width,viscous_wall=viscous, &
      internal_energy_density=internal_energy_density,internal_energy_ghost=internal_energy_ghost)
    if(fail_new_halo.and.calls==3) ok=.false.
  end subroutine
end program
