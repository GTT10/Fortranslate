program test_mpi_eb_tree_transfer_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties,universal_gas_constant
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_elementary_transport
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use dual_energy_mod, only: initialize_nasa7_internal_energy_patch
  use reactive_1d_mod, only: reactive_nvar
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,fill_reactive_eb_tree_halo_3d, &
    build_eb_tree_extended_geometry_3d,reactive_eb_tree_timestep_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use mpi_amr_patch_tree_distribution_3d_mod, only: scatter_mpi_amr_tree_3d,gather_mpi_amr_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,integrate_mpi_eb_tree_3d
  use mpi_eb_tree_initialize_3d_mod, only: initialize_mpi_eb_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use mpi_eb_tree_timestep_3d_mod, only: mpi_eb_tree_timestep_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse,geometry(5),original(5)
  type(eb_reactive_tree_3d) :: tree,serial
  type(rotating_frame_3d) :: frame,mismatched_frame
  type(eb_geometry_3d), allocatable :: initialized_geometry(:)
  type(reactive_eb_3d_config) :: config
  type(amr_tree_state_3d) :: global,gt,local,lt,result,rt,saved,st
  type(amr_tree_state_3d) :: initialized,initialized_t,serial_held,initialized_held,saved_initialized_held
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:),expected(:),actual(:)
  real(dp) :: error,position(3),potential,cp,cv,h,e,entropy,expected_e
  integer :: nv,node,count,rank,nranks,status,owners(5),root,mode,i,j,k,sp,layout
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,status)
  root=nranks-1;count=5
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(expected(nv),actual(nv))
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok)
  call require(ok,'child left')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok)
  call require(ok,'child right')
  call initialize_amr_patch_3d(4,4,4,3,4,3,4,2,3,2,boxes(3),ok)
  call require(ok,'grandchild left')
  call initialize_amr_patch_3d(4,4,4,1,2,3,4,2,3,2,boxes(4),ok)
  call require(ok,'grandchild right')
  if (count==5) then
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2,3],boxes,topology,ok)
  else
    call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],boxes(:3),topology,ok)
  end if
  call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'z',.31_dp,finest,ok)
  call require(ok,'finest geometry')
  call coarsen_eb_geometry_3d(finest,2,middle,ok)
  call require(ok,'middle geometry')
  call coarsen_eb_geometry_3d(middle,2,coarse,ok)
  call require(ok,'coarse geometry')
  allocate(u(nv,4,4,4),t(4,4,4))
  config%nx=4;config%ny=4;config%nz=4
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,u,t,ok)
  call require(ok,'root state')
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,tree,ok)
  call require(ok,'tree initialization')

  do node=2,5
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
  end do
  serial=tree
  call average_down_reactive_eb_tree_3d(species,serial,ok)
  call require(ok,'serial average')
  call integrate_reactive_eb_tree_3d(serial,expected,ok)
  call require(ok,'serial inventory')
  do node=1,5
    geometry(node)=tree%nodes(node)%geometry
  end do
  original=geometry
  if (rank==root) then
    allocate(global%fields(5),gt%fields(5))
    do node=1,5
      global%fields(node)%values=tree%nodes(node)%state
      allocate(gt%fields(node)%values(1,4,4,4))
      gt%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
    end do
  end if
  do mode=1,2
    owners=[(mod(node-1,nranks),node=1,5)]
    if (mode==2) owners=0
    call scatter_mpi_amr_tree_3d(topology,owners,global,local,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter U')
    call scatter_mpi_amr_tree_3d(topology,owners,gt,lt,MPI_COMM_WORLD,root,ok)
    call require(ok,'scatter T')
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(ok,'MPI average')
    do node=1,5
      call require(allocated(local%fields(node)%values).eqv.(rank==owners(node)),'owner-local U')
      call require(allocated(lt%fields(node)%values).eqv.(rank==owners(node)),'owner-local T')
    end do
    call gather_mpi_amr_tree_3d(topology,owners,local,result,MPI_COMM_WORLD,root,ok)
    call require(ok,'gather U')
    call gather_mpi_amr_tree_3d(topology,owners,lt,rt,MPI_COMM_WORLD,root,ok)
    call require(ok,'gather T')
    if (rank==root) then
      do node=1,5
        call require(all(result%fields(node)%values==serial%nodes(node)%state),'exact serial U including solids')
        call require(all(rt%fields(node)%values(1,:,:,:)==serial%nodes(node)%temperature),'exact serial T')
      end do
    end if
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,actual,MPI_COMM_WORLD,ok)
    call require(ok,'MPI inventory')
    error=maxval(abs(actual-expected)/max(1._dp,abs(expected)))
    call require(error<3.e-14_dp,'composite inventory')
    if (rank==0) print *, 'ranks, mode, inventory error: ',nranks,mode,error
    saved=local;st=lt
    ! Reject a late child after another child has already restricted privately.
    if (rank==owners(4)) local%fields(4)%values(1,4,4,4)=-1._dp
    saved=local
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(.not.ok,'bad fine EOS rejected collectively')
    call unchanged()
    local=saved
    if (rank==owners(4)) local%fields(4)%values(1,4,4,4)=serial%nodes(4)%state(1,4,4,4)
    saved=local
    if (rank==0) then
      do node=1,5
        geometry(node)%x_lower=geometry(node)%x_lower+.01_dp
        geometry(node)%x_upper=geometry(node)%x_upper+.01_dp
        call require(geometry(node)%is_valid(),'locally valid translated geometry')
      end do
    end if
    call average_down_mpi_eb_tree_3d(species,topology,owners,geometry,local,lt,MPI_COMM_WORLD,ok)
    call require(.not.ok,'geometry mismatch rejected collectively')
    call unchanged()
    geometry=original
    if (rank==owners(5)) local%fields(5)%values(1,4,4,4)=ieee_value(0._dp,ieee_quiet_nan)
    actual=17._dp
    call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,actual,MPI_COMM_WORLD,ok)
    call require(.not.ok.and.all(actual==17._dp),'invalid exposed payload preserves inventory')
  end do
  frame%axis=1;frame%angular_speed=100._dp;frame%axis_point=[10._dp,10._dp,10._dp]
  do k=1,4
  do j=1,4
  do i=1,4
    if(coarse%volume_fraction(i,j,k)<=0) cycle
    position=(real([i,j,k],dp)-.5_dp)/4
    call rotating_potential_3d(frame,position,potential,ok)
    call require(ok,'rotating root potential')
    u(5,i,j,k)=u(5,i,j,k)-u(1,i,j,k)*potential
  end do
  end do
  end do
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,serial,ok,frame)
  call require(ok,'serial rotating initialization')
  owners=[(mod(node-1,nranks),node=1,5)]
  if(nranks>1) then
    mismatched_frame=frame
    if(rank==1) mismatched_frame%angular_speed=mismatched_frame%angular_speed+1._dp
    call initialize_mpi_eb_tree_3d(species,topology,owners,finest,u,t, &
      initialized_geometry,initialized,initialized_t,MPI_COMM_WORLD,ok,mismatched_frame)
    call require(.not.ok,'rank-mismatched initialization frame rejected')
    call require(.not.allocated(initialized_geometry),'mismatched initialization rollback')
  end if
  call initialize_mpi_eb_tree_3d(species,topology,owners,finest,u,t, &
    initialized_geometry,initialized,initialized_t,MPI_COMM_WORLD,ok,frame)
  call require(ok,'MPI rotating initialization')
  do node=1,5
    if(owners(node)/=rank) cycle
    call require(all(initialized%fields(node)%values==serial%nodes(node)%state),'rotating initial U serial parity')
    call require(all(initialized_t%fields(node)%values(1,:,:,:)==serial%nodes(node)%temperature), &
      'rotating initial T serial parity')
  end do
  call initialize_reactive_eb_tree_3d(species,topology,finest,u,t,serial,ok,frame,serial_held)
  call require(ok,'serial rotating held initialization')
  do layout=1,2
    owners=[(mod(node-1,nranks),node=1,5)]
    if(layout==2) owners=0
    call initialize_mpi_eb_tree_3d(species,topology,owners,finest,u,t, &
      initialized_geometry,initialized,initialized_t,MPI_COMM_WORLD,ok,frame,initialized_held)
    call require(ok,'native held initialization on distributed and owner-zero layouts')
    do node=1,5
      if(owners(node)/=rank) then
        call require(.not.allocated(initialized_held%fields(node)%values),'idle ranks retain no held payload')
        cycle
      endif
      call require(all(initialized_held%fields(node)%values==serial_held%fields(node)%values), &
        'native serial/MPI held initialization agrees exactly')
      do k=1,4
      do j=1,4
      do i=1,4
        if(initialized_geometry(node)%volume_fraction(i,j,k)<=0) then
          call require(initialized_held%fields(node)%values(1,i,j,k)==0,'initial covered scalar is zero')
          cycle
        endif
        expected_e=0
        do sp=1,size(species)
          call nasa7_mass_properties(species(sp),initialized_t%fields(node)%values(1,i,j,k), &
            cp,cv,h,e,entropy,ok)
          call require(ok,'independent species EOS initialization reference')
          expected_e=expected_e+initialized%fields(node)%values(5+sp,i,j,k)*e
        enddo
        call require(abs(initialized_held%fields(node)%values(1,i,j,k)-expected_e)< &
          3.e-14_dp*max(1._dp,abs(expected_e)),'held initialization matches independent species energy sum')
      enddo
      enddo
      enddo
    enddo
    if(nranks>1) then
      saved_initialized_held=initialized_held
      if(rank==0) then
        call initialize_mpi_eb_tree_3d(species,topology,owners,finest,u,t, &
          initialized_geometry,initialized,initialized_t,MPI_COMM_WORLD,ok,frame,initialized_held)
      else
        call initialize_mpi_eb_tree_3d(species,topology,owners,finest,u,t, &
          initialized_geometry,initialized,initialized_t,MPI_COMM_WORLD,ok,frame)
      endif
      call require(.not.ok,'rank-disagreeing held initialization presence rejects collectively')
      do node=1,5
        if(owners(node)/=rank) cycle
        call require(all(initialized_held%fields(node)%values==saved_initialized_held%fields(node)%values).and. &
          all(initialized%fields(node)%values==serial%nodes(node)%state).and. &
          all(initialized_t%fields(node)%values(1,:,:,:)==serial%nodes(node)%temperature), &
          'held presence rejection preserves prior gas/T/scalar outputs')
      enddo
    endif
  enddo
  call check_held_tree_average()
  call check_held_root_seed()
  if (rank==0) print *, 'PASS: MPI EB restriction and composite inventory'
  call MPI_Finalize(status)
contains
  subroutine check_held_root_seed()
    type(eb_reactive_tree_3d) :: reference,trial,saved_tree
    type(amr_tree_state_3d) :: reference_h,trial_h,saved_h,gas,temp,held,g0,t0,h0
    type(eb_geometry_3d), allocatable :: geom(:),geom0(:)
    type(rotating_frame_3d), allocatable :: rotation
    real(dp), allocatable :: seed_u(:,:,:,:),seed_t(:,:,:),empty_u(:,:,:,:),empty_t(:,:,:),bad_t(:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: owner(count),p,x,y,z,condition,layout,donor(3),factor
    real(dp) :: rho,phi,position(3),expected_e,cp,cv,h,e,s,thermal_error
    logical :: good
    allocate(empty_u(0,0,0,0),empty_t(0,0,0))
    call amr_tree_integer_layout_3d(topology,origin,domain,ratios,good)
    call require(good,'independent root seed donor coordinates')
    do condition=1,2
      if(condition==2) then
        allocate(rotation);rotation%axis=2;rotation%angular_speed=1.e10_dp
        rotation%axis_point=[10._dp,10._dp,10._dp]
      endif
      seed_u=tree%nodes(1)%state;seed_t=tree%nodes(1)%temperature
      do z=1,4
      do y=1,4
      do x=1,4
        if(coarse%volume_fraction(x,y,z)<=0) cycle
        rho=seed_u(1,x,y,z);seed_t(x,y,z)=700._dp+5*x+3*y+z
        seed_u(2:4,x,y,z)=[rho*1.e10_dp,0._dp,0._dp]
        seed_u(5,x,y,z)=.5_dp*rho*1.e20_dp
        if(allocated(rotation)) then
          position=(real([x,y,z],dp)-.5_dp)/4
          call rotating_potential_3d(rotation,position,phi,good)
          call require(good,'root seed rotating potential')
          seed_u(5,x,y,z)=seed_u(5,x,y,z)-rho*phi
        endif
      enddo
      enddo
      enddo
      call initialize_reactive_eb_tree_3d(species,topology,finest,seed_u,seed_t,reference,good,rotation,reference_h)
      call require(good,'held native root seed survives unresolved total-minus-kinetic recovery')
      do p=1,count
        factor=product(ratios(2:topology%nodes(p)%level+1))
        do z=1,4
        do y=1,4
        do x=1,4
          if(reference%nodes(p)%geometry%volume_fraction(x,y,z)<=0) then
            call require(reference_h%fields(p)%values(1,x,y,z)==0,'root seeded covered UEINT remains zero')
            cycle
          endif
          donor=1+int((origin(:,p)+int([x,y,z]-1,int64))/int(factor,int64))
          call require(abs(reference%nodes(p)%temperature(x,y,z)-seed_t(donor(1),donor(2),donor(3)))<1.e-10_dp, &
            'independent root temperature retained throughout all descendants')
          expected_e=0
          do sp=1,size(species)
            call nasa7_mass_properties(species(sp),seed_t(donor(1),donor(2),donor(3)),cp,cv,h,e,s,good)
            call require(good,'independent root temperature EOS')
            expected_e=expected_e+reference%nodes(p)%state(5+sp,x,y,z)*e
          enddo
          call require(abs(reference_h%fields(p)%values(1,x,y,z)-expected_e)< &
            3.e-14_dp*max(1._dp,abs(expected_e)),'root UEINT equals independent native species sum')
        enddo
        enddo
        enddo
      enddo
      ! The same gas-only constructor cannot infer the independently supplied thermal seed.
      trial=reference;saved_tree=trial
      call initialize_reactive_eb_tree_3d(species,topology,finest,seed_u,seed_t,trial,good,rotation)
      if(good) then
        ! Native formation energies can map the canceled gas energy to another
        ! admissible temperature. Acceptance alone cannot retain the seed.
        thermal_error=0
        do p=1,count
          thermal_error=max(thermal_error,maxval(abs(trial%nodes(p)%temperature-reference%nodes(p)%temperature), &
            mask=reference%nodes(p)%geometry%volume_fraction>0))
        enddo
        call require(thermal_error>1._dp,'gas-only canceled energy cannot preserve the supplied thermal seed')
      else
        do p=1,count
          call require(all(trial%nodes(p)%state==saved_tree%nodes(p)%state).and. &
            all(trial%nodes(p)%temperature==saved_tree%nodes(p)%temperature),'gas-only seed refusal is atomic')
        enddo
      endif
      trial=reference
      bad_t=seed_t;bad_t(4,4,4)=1.e30_dp;trial_h=reference_h;saved_h=trial_h
      call initialize_reactive_eb_tree_3d(species,topology,finest,seed_u,bad_t,trial,good,rotation,trial_h)
      call require(.not.good,'late finite thermal seed outside EOS bounds refuses initialization')
      do p=1,count
        call require(all(trial%nodes(p)%state==saved_tree%nodes(p)%state).and. &
          all(trial%nodes(p)%temperature==saved_tree%nodes(p)%temperature).and. &
          all(trial_h%fields(p)%values==saved_h%fields(p)%values),'bad root thermal seed preserves gas/T/held')
      enddo
      do layout=1,2
        owner=[(mod(p+nranks-2,nranks),p=1,count)];if(layout==2) owner=0
        if(rank==owner(1)) then
          call initialize_mpi_eb_tree_3d(species,topology,owner,finest,seed_u,seed_t, &
            geom,gas,temp,MPI_COMM_WORLD,good,rotation,held)
        else
          call initialize_mpi_eb_tree_3d(species,topology,owner,finest,empty_u,empty_t, &
            geom,gas,temp,MPI_COMM_WORLD,good,rotation,held)
        endif
        call require(good,'owner-local held root initialization including empty seed on other ranks')
        do p=1,count
          call require(matching_eb_metrics_3d(geom(p),reference%nodes(p)%geometry),'replicated seed geometry')
          if(owner(p)/=rank) then
            call require(.not.allocated(held%fields(p)%values),'idle root seed rank has no scalar payload')
            cycle
          endif
          call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
            all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
            all(held%fields(p)%values==reference_h%fields(p)%values),'high kinetic root seed exact serial/MPI parity')
        enddo
        call check_held_halo(reference,reference_h,geom,gas,temp,held,owner,rotation)
        call check_held_timestep(reference,reference_h,geom,gas,temp,held,owner,rotation)
        g0=gas;t0=temp;h0=held;geom0=geom
        if(rank==owner(1)) then
          call initialize_mpi_eb_tree_3d(species,topology,owner,finest,seed_u,bad_t, &
            geom,gas,temp,MPI_COMM_WORLD,good,rotation,held)
        else
          call initialize_mpi_eb_tree_3d(species,topology,owner,finest,empty_u,empty_t, &
            geom,gas,temp,MPI_COMM_WORLD,good,rotation,held)
        endif
        call require(.not.good,'owner-only late thermal seed EOS refusal is collective')
        do p=1,count
          call require(matching_eb_metrics_3d(geom(p),geom0(p)),'refused root seed preserves replicated geometry')
          if(owner(p)/=rank) cycle
          call require(all(gas%fields(p)%values==g0%fields(p)%values).and. &
            all(temp%fields(p)%values==t0%fields(p)%values).and.all(held%fields(p)%values==h0%fields(p)%values), &
            'failed high kinetic seed restores prior gas/T/held outputs')
        enddo
      enddo
    enddo
  end subroutine check_held_root_seed
  subroutine check_held_halo(seed,seed_h,geom,gas,temp,held,owner,rotation)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_h,gas,temp,held
    type(eb_geometry_3d), intent(in) :: geom(:)
    integer, intent(in) :: owner(:)
    type(rotating_frame_3d), optional, intent(in) :: rotation
    type(amr_tree_state_3d) :: gh,gt,hh,gv,metrics,gh0,gt0,hh0,late_h,late_serial
    type(amr_tree_state_3d), allocatable :: as,ae,ats,ate,ahs,ahe,serial_h0,serial_h1
    type(eb_reactive_tree_3d), allocatable :: serial0,serial1
    type(eb_geometry_3d), allocatable :: extended(:)
    type(eb_geometry_3d) :: padded
    real(dp), allocatable :: u(:,:,:,:),t(:,:,:),h(:,:,:),vf(:,:,:),u0(:,:,:,:),t0(:,:,:),h0(:,:,:)
    real(dp), allocatable :: mach(:),tin(:),inlet(:,:),fractions(:),serial_fraction(:)
    character(len=16), allocatable :: boundary(:)
    character(len=16) :: faces(6)
    integer :: mode,p,width,species_index
    real(dp) :: cp,cv,enthalpy,energy,entropy,expected_h,pressure,expected_pressure
    logical :: good
    do mode=1,4
      width=1;if(mode==4) width=2
      faces='outflow';if(mode==2) faces='no_slip_wall'
      if(allocated(u)) deallocate(u,t,h,vf)
      allocate(u(nv,1-width:4+width,1-width:4+width,1-width:4+width), &
        t(1-width:4+width,1-width:4+width,1-width:4+width), &
        h(1-width:4+width,1-width:4+width,1-width:4+width), &
        vf(1-width:4+width,1-width:4+width,1-width:4+width))
      if(mode<3) then
        call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
          boundary_faces=faces,frame=rotation,internal_energy=held,ghost_internal_energy=hh)
      else if(mode==3) then
        call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
          boundary_faces=faces,frame=rotation,internal_energy=held,ghost_internal_energy=hh,ghost_volume_fraction=gv)
      else
        call build_axis_plane_eb_geometry_3d(32,32,32,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp,'z',.31_dp,padded,good)
        call require(good,'padded finest geometry for width-two held halos')
        call build_eb_tree_extended_geometry_3d(topology,geom,padded,extended,good,halo_width=width)
        call require(good,'held wide halo finest-resolution metric construction')
        allocate(metrics%fields(count))
        do p=1,count
          if(owner(p)/=rank) cycle
          allocate(metrics%fields(p)%values(1,1-width:4+width,1-width:4+width,1-width:4+width))
          metrics%fields(p)%values(1,:,:,:)=extended(p)%volume_fraction
        enddo
        call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
          boundary_faces=faces,frame=rotation,internal_energy=held,ghost_internal_energy=hh, &
          ghost_volume_fraction=gv,halo_width=width,target_volume_fraction=metrics)
      endif
      call require(good,'native high kinetic held halo accepts owner/idle, physical and peer donors')
      do p=1,count
        if(owner(p)/=rank) then
          call require(.not.allocated(hh%fields(p)%values),'idle rank has no held halo payload')
          cycle
        endif
        u=17;t=19;h=23;vf=29
        if(mode<3) then
          call fill_reactive_eb_tree_halo_3d(species,seed,p,u,t,good,boundary_faces=faces, &
            frame=rotation,internal_energy=seed_h,ghost_internal_energy=h)
        else if(mode==3) then
          call fill_reactive_eb_tree_halo_3d(species,seed,p,u,t,good,boundary_faces=faces, &
            frame=rotation,internal_energy=seed_h,ghost_internal_energy=h,ghost_volume_fraction=vf)
        else
          call fill_reactive_eb_tree_halo_3d(species,seed,p,u,t,good,boundary_faces=faces, &
            frame=rotation,internal_energy=seed_h,ghost_internal_energy=h,ghost_volume_fraction=vf, &
            halo_width=width,target_volume_fraction=metrics%fields(p)%values(1,:,:,:))
        endif
        call require(good,'serial paired native held halo recovers supplied thermal seed')
        call require(all(u==gh%fields(p)%values).and.all(t==gt%fields(p)%values(1,:,:,:)).and. &
          all(h==hh%fields(p)%values(1,:,:,:)),'exact serial/MPI gas/T/UEINT halo parity including peer overlay')
        if(mode>=3) call require(all(vf==gv%fields(p)%values(1,:,:,:)),'exact serial/MPI held halo volume fractions')
        call require(all(t>500._dp.or.u(1,:,:,:)==0),'held halo EOS retains the high kinetic thermal seed')
      enddo
      if(allocated(metrics%fields)) deallocate(metrics%fields)
      gh0=gh;gt0=gt;hh0=hh
      ! A finite bad scalar in a late uncovered donor rejects the private
      ! restriction after earlier descendants, preserving accepted halos.
      late_h=held;late_serial=seed_h
      if(owner(2)==rank) late_h%fields(2)%values(1,4,4,4)=1.e30_dp
      late_serial%fields(2)%values(1,4,4,4)=1.e30_dp
      call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
        boundary_faces=faces,frame=rotation,internal_energy=late_h,ghost_internal_energy=hh)
      call require(.not.good,'finite late held donor EOS refusal is collective')
      do p=1,count
        if(owner(p)/=rank) cycle
        call require(all(gh%fields(p)%values==gh0%fields(p)%values).and. &
          all(gt%fields(p)%values==gt0%fields(p)%values).and.all(hh%fields(p)%values==hh0%fields(p)%values), &
          'failed held halo leaves every accepted owner payload unchanged')
      enddo
      ! Width-one refusal must not replace previously accepted wide buffers.
      u0=u;t0=t;h0=h
      if(width==1) then
        call fill_reactive_eb_tree_halo_3d(species,seed,1,u,t,good,boundary_faces=faces, &
          frame=rotation,internal_energy=late_serial,ghost_internal_energy=h)
        call require(.not.good.and.all(u==u0).and.all(t==t0).and.all(h==h0),'serial late held halo refusal is atomic')
      endif
      if(nranks>1) then
        if(rank==nranks-1) then
          call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
            boundary_faces=faces,frame=rotation)
        else
          call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
            boundary_faces=faces,frame=rotation,internal_energy=held,ghost_internal_energy=hh)
        endif
        call require(.not.good,'rank-disagreeing held halo presence refuses collectively')
      endif
    enddo
    deallocate(u,t,h,vf)
    allocate(u(nv,0:5,0:5,0:5),t(0:5,0:5,0:5),h(0:5,0:5,0:5))
    do mode=5,7
      if(mode==5) then
        allocate(boundary(6),tin(6),inlet(nv,6));boundary='outflow';boundary(1)='inflow';tin=750._dp
        do p=1,6
          inlet(:,p)=seed%nodes(1)%state(:,4,4,4)
          inlet(5,p)=.5_dp*inlet(1,p)*1.e20_dp
        enddo
      else if(mode==6) then
        deallocate(boundary,inlet)
        allocate(mach(6));mach=-1;mach(1)=.3_dp
      else
        deallocate(mach,tin)
        allocate(fractions(2),as,ae,ats,ate,ahs,ahe,serial0,serial1,serial_h0,serial_h1)
        fractions=.25_dp;as=gas;ae=gas;ats=temp;ate=temp;ahs=held;ahe=held
        serial0=seed;serial1=seed;serial_h0=seed_h;serial_h1=seed_h
        do p=1,count
          serial_h1%fields(p)%values=1.1_dp*serial_h1%fields(p)%values
          if(owner(p)==rank) ahe%fields(p)%values=1.1_dp*ahe%fields(p)%values
        enddo
      endif
      call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
        boundary_faces=boundary,inflow_state=inlet,inflow_mach=mach,inflow_temperature=tin,frame=rotation, &
        internal_energy=held,ghost_internal_energy=hh,restrict_snapshot=.false., &
        ancestor_start=as,ancestor_end=ae,temperature_start=ats,temperature_end=ate,ancestor_fraction=fractions, &
        ancestor_internal_energy_start=ahs,ancestor_internal_energy_end=ahe)
      call require(good,'held halo accepts high kinetic prescribed/Mach inlet and timed ancestor donors')
      do p=1,count
        if(owner(p)/=rank) cycle
        if(mode==7) serial_fraction=fractions(:seed%topology%nodes(p)%level)
        call fill_reactive_eb_tree_halo_3d(species,seed,p,u,t,good,boundary_faces=boundary, &
          inflow_state=inlet,inflow_mach=mach,inflow_temperature=tin,frame=rotation, &
          internal_energy=seed_h,ghost_internal_energy=h,restrict_snapshot=.false., &
          ancestor_start=serial0,ancestor_end=serial1,ancestor_fraction=serial_fraction, &
          ancestor_internal_energy_start=serial_h0,ancestor_internal_energy_end=serial_h1)
        if(.not.good.or.any(u/=gh%fields(p)%values).or.any(t/=gt%fields(p)%values(1,:,:,:)).or. &
          any(h/=hh%fields(p)%values(1,:,:,:))) then
          print *, 'Held inlet/temporal halo mismatch rank/mode/node/rotation/accepted: ',rank,mode,p,present(rotation),good
          print *, 'Held inlet/temporal halo max gas/T/H difference: ',maxval(abs(u-gh%fields(p)%values)), &
            maxval(abs(t-gt%fields(p)%values(1,:,:,:))),maxval(abs(h-hh%fields(p)%values(1,:,:,:)))
          print *, 'Held halo largest T/H array indices: ',maxloc(abs(t-gt%fields(p)%values(1,:,:,:))), &
            maxloc(abs(h-hh%fields(p)%values(1,:,:,:)))
        endif
        call require(good.and.all(u==gh%fields(p)%values).and.all(t==gt%fields(p)%values(1,:,:,:)).and. &
          all(h==hh%fields(p)%values(1,:,:,:)),'inlet/temporal held halo exactly matches native serial and MPI')
        if(p==1.and.mode<7) then
          call require(abs(t(0,2,4)-750._dp)<1.e-10_dp,'held inlet retains its prescribed thermal seed')
          expected_h=0;pressure=0;expected_pressure=0
          do species_index=1,size(species)
            call nasa7_mass_properties(species(species_index),750._dp,cp,cv,enthalpy,energy,entropy,good)
            call require(good,'independent held inlet species EOS')
            expected_h=expected_h+u(5+species_index,0,2,4)*energy
            pressure=pressure+u(5+species_index,0,2,4)*750._dp*universal_gas_constant/species(species_index)%molecular_weight
            expected_pressure=expected_pressure+seed%nodes(1)%state(5+species_index,1,2,4)* &
              seed%nodes(1)%temperature(1,2,4)*universal_gas_constant/species(species_index)%molecular_weight
          enddo
          call require(abs(h(0,2,4)-expected_h)<3.e-14_dp*max(1._dp,abs(expected_h)), &
            'inlet held UEINT matches independent EOS after rotating energy conversion')
          if(mode==6) call require(abs(pressure-expected_pressure)<3.e-14_dp*max(1._dp,abs(expected_pressure)), &
            'Mach inlet extrapolates pressure from held donor energy without kinetic cancellation')
        endif
        if(p==2.and.mode==7) then
          expected_h=1.025_dp*seed_h%fields(1)%values(1,1,3,3)
          call require(abs(h(1,5,4)-expected_h)<3.e-14_dp*max(1._dp,abs(expected_h)), &
            'coarse donor UEINT uses independent actual ancestor time interpolation')
          call require(abs(t(1,5,4)-seed%nodes(1)%temperature(1,3,3))>.1_dp, &
            'timed ancestor scalar changes the recovered ghost temperature')
        endif
      enddo
      if(mode==5) then
        gh0=gh;gt0=gt;hh0=hh;tin(1)=1.e30_dp
        call fill_mpi_eb_tree_halo_3d(species,topology,owner,geom,gas,temp,gh,gt,MPI_COMM_WORLD,good, &
          boundary_faces=boundary,inflow_state=inlet,inflow_temperature=tin,frame=rotation, &
          internal_energy=held,ghost_internal_energy=hh,restrict_snapshot=.false.)
        call require(.not.good,'EOS-invalid prescribed thermal seed rejects halo collectively')
        do p=1,count
          if(owner(p)/=rank) cycle
          call require(all(gh%fields(p)%values==gh0%fields(p)%values).and. &
            all(gt%fields(p)%values==gt0%fields(p)%values).and.all(hh%fields(p)%values==hh0%fields(p)%values), &
            'failed inlet thermal recovery preserves accepted gas/T/held halos')
        enddo
        tin=750._dp
      endif
    enddo
  end subroutine check_held_halo
  subroutine check_held_timestep(seed,seed_h,geom,gas,temp,held,owner,rotation)
    type(eb_reactive_tree_3d), intent(in) :: seed
    type(amr_tree_state_3d), intent(in) :: seed_h,gas,temp,held
    type(eb_geometry_3d), intent(in) :: geom(:)
    integer, intent(in) :: owner(:)
    type(rotating_frame_3d), optional, intent(in) :: rotation
    type(gas_transport_species), allocatable :: transport(:)
    type(amr_tree_state_3d) :: bad_h,serial_bad
    character(len=16), allocatable :: faces(:)
    real(dp), allocatable :: wall_t(:)
    real(dp) :: coefficients(4),donors(4),maximum_d,cp,cv,enthalpy,energy,entropy,cp_mix,cv_mix,r_mix
    real(dp) :: rho,y,tt,sound,velocity(3),dmax,metric,maximum_metric,rate,maximum_rate,wall_rate,max_wall
    real(dp) :: reference_dt,serial_dt,mpi_dt,factor,transport_factor,spacing(3),volume
    integer :: mode,clock,p,x,j,z,s
    logical :: good,viscosity,conduction,diffusion,recursive,coupled
    call load_h2o2_elementary_transport(transport,good)
    call require(good,'native held timestep transport database')
    coefficients=[1.e8_dp,1.e18_dp,1.e7_dp,1.e8_dp]
    do mode=0,4
      viscosity=mode==1.or.mode==4;conduction=mode==2.or.mode==4;diffusion=mode==3.or.mode==4
      if(mode==1) then
        allocate(faces(6));faces='outflow'
      endif
      if(mode==4) then
        faces(1)='isothermal_wall';allocate(wall_t(6));wall_t=800
      endif
      ! Independent NASA species sums give sound and each diffusivity. Hidden
      ! parent cells and all donor levels participate in the root interval.
      donors=0;maximum_d=0
      do p=1,count
      do z=1,4
      do j=1,4
      do x=1,4
        if(geom(p)%volume_fraction(x,j,z)<=0) cycle
        rho=seed%nodes(p)%state(1,x,j,z);tt=seed%nodes(p)%temperature(x,j,z)
        velocity=seed%nodes(p)%state(2:4,x,j,z)/rho;cp_mix=0;cv_mix=0;r_mix=0
        do s=1,size(species)
          call nasa7_mass_properties(species(s),tt,cp,cv,enthalpy,energy,entropy,good)
          call require(good,'independent timestep native species heat capacities')
          y=seed%nodes(p)%state(5+s,x,j,z)/sum(seed%nodes(p)%state(6:,x,j,z))
          cp_mix=cp_mix+y*cp;cv_mix=cv_mix+y*cv;r_mix=r_mix+y*universal_gas_constant/species(s)%molecular_weight
        enddo
        sound=sqrt(cp_mix/cv_mix*r_mix*tt)
        donors(:3)=max(donors(:3),abs(velocity)+sound);donors(4)=max(donors(4),sum(abs(velocity))+sound)
        dmax=0
        if(viscosity) dmax=max(dmax,((4._dp/3._dp)*coefficients(1)+coefficients(3))/rho)
        if(conduction) dmax=max(dmax,coefficients(2)/(rho*cv_mix))
        if(diffusion) dmax=max(dmax,coefficients(4)/rho)
        maximum_d=max(maximum_d,dmax)
      enddo
      enddo
      enddo
      enddo
      do clock=1,3
        recursive=clock>1;coupled=clock/=3;reference_dt=huge(1._dp)
        do p=1,count
          associate(g=>geom(p))
          spacing=[g%dx,g%dy,g%dz];volume=product(spacing);maximum_rate=0;maximum_metric=0;max_wall=0
          do z=1,4
          do j=1,4
          do x=1,4
            if(g%volume_fraction(x,j,z)<=0) cycle
            rate=(donors(1)*(g%x_face_fraction(x-1,j,z)+g%x_face_fraction(x,j,z))/g%dx+ &
              donors(2)*(g%y_face_fraction(x,j-1,z)+g%y_face_fraction(x,j,z))/g%dy+ &
              donors(3)*(g%z_face_fraction(x,j,z-1)+g%z_face_fraction(x,j,z))/g%dz+ &
              donors(4)*g%boundary_area(x,j,z)/volume)/g%volume_fraction(x,j,z)
            metric=.5_dp*((g%x_face_fraction(x-1,j,z)+g%x_face_fraction(x,j,z))/g%dx**2+ &
              (g%y_face_fraction(x,j-1,z)+g%y_face_fraction(x,j,z))/g%dy**2+ &
              (g%z_face_fraction(x,j,z-1)+g%z_face_fraction(x,j,z))/g%dz**2)/g%volume_fraction(x,j,z)
            maximum_rate=max(maximum_rate,rate);maximum_metric=max(maximum_metric,metric)
            if(mode==4.and.x==1.and.g%x_lower==geom(1)%x_lower) then
              rho=seed%nodes(p)%state(1,x,j,z);tt=seed%nodes(p)%temperature(x,j,z);cv_mix=0
              do s=1,size(species)
                call nasa7_mass_properties(species(s),tt,cp,cv,enthalpy,energy,entropy,good)
                call require(good,'independent held wall timestep heat capacity')
                y=seed%nodes(p)%state(5+s,x,j,z)/sum(seed%nodes(p)%state(6:,x,j,z));cv_mix=cv_mix+y*cv
              enddo
              wall_rate=2*coefficients(2)/(rho*cv_mix)*g%x_face_fraction(0,j,z)/g%dx**2/g%volume_fraction(x,j,z)
              max_wall=max(max_wall,wall_rate)
            endif
          enddo
          enddo
          enddo
          factor=1;if(recursive) factor=real(2**topology%nodes(p)%level,dp)
          transport_factor=factor;if(.not.coupled) transport_factor=factor*factor
          reference_dt=min(reference_dt,.4_dp/maximum_rate*factor)
          if(mode>0) reference_dt=min(reference_dt,.2_dp/(maximum_d*maximum_metric)*transport_factor)
          if(max_wall>0) reference_dt=min(reference_dt,.2_dp/max_wall*transport_factor)
          end associate
        enddo
        serial_dt=17;mpi_dt=19
        call reactive_eb_tree_timestep_3d(species,transport,seed,.4_dp,.2_dp,mode>0,viscosity,conduction,diffusion, &
          serial_dt,good,subcycling=recursive,coupled_subcycling=coupled,frame=rotation,boundary_faces=faces, &
          wall_temperature=wall_t,constant_coefficients=coefficients,internal_energy=seed_h)
        call require(good.and.abs(serial_dt-reference_dt)<3.e-14_dp*reference_dt, &
          'serial held hydro/transport/wall root timestep matches independent geometry/NASA bound')
        call mpi_eb_tree_timestep_3d(species,transport,topology,owner,geom,gas,temp,.4_dp,.2_dp,mode>0, &
          viscosity,conduction,diffusion,mpi_dt,MPI_COMM_WORLD,good,subcycling=recursive,coupled_subcycling=coupled, &
          frame=rotation,boundary_faces=faces,wall_temperature=wall_t,constant_coefficients=coefficients,internal_energy=held)
        call require(good.and.mpi_dt==serial_dt,'held root timestep exactly matches serial/MPI including idle ranks')
      enddo
    enddo
    bad_h=held;serial_bad=seed_h
    if(owner(2)==rank) bad_h%fields(2)%values(1,4,4,4)=1.e30_dp
    serial_bad%fields(2)%values(1,4,4,4)=1.e30_dp;serial_dt=17;mpi_dt=19
    call reactive_eb_tree_timestep_3d(species,transport,seed,.4_dp,.2_dp,.false.,.false.,.false.,.false., &
      serial_dt,good,frame=rotation,internal_energy=serial_bad)
    call require(.not.good.and.serial_dt==17,'late held EOS refusal preserves serial caller timestep')
    call mpi_eb_tree_timestep_3d(species,transport,topology,owner,geom,gas,temp,.4_dp,.2_dp, &
      .false.,.false.,.false.,.false.,mpi_dt,MPI_COMM_WORLD,good,frame=rotation,internal_energy=bad_h)
    call require(.not.good.and.mpi_dt==19,'late owner held EOS refusal preserves every MPI caller timestep')
    if(nranks>1) then
      if(rank==nranks-1) then
        call mpi_eb_tree_timestep_3d(species,transport,topology,owner,geom,gas,temp,.4_dp,.2_dp, &
          .false.,.false.,.false.,.false.,mpi_dt,MPI_COMM_WORLD,good,frame=rotation)
      else
        call mpi_eb_tree_timestep_3d(species,transport,topology,owner,geom,gas,temp,.4_dp,.2_dp, &
          .false.,.false.,.false.,.false.,mpi_dt,MPI_COMM_WORLD,good,frame=rotation,internal_energy=held)
      endif
      call require(.not.good.and.mpi_dt==19,'rank-disagreeing held timestep presence refuses collectively and atomically')
    endif
  end subroutine check_held_timestep
  subroutine check_held_tree_average()
    type(eb_reactive_tree_3d) :: seed,reference,serial_trial,serial_saved
    type(amr_tree_state_3d) :: seed_held,reference_held,sh,ss,gas,temp,held,g0,t0,h0
    type(rotating_frame_3d), allocatable :: rotation
    integer :: p,x,y,z,condition,layout,owner(count)
    real(dp) :: rho,phi,position(3)
    logical :: good
    do condition=1,2
      if(condition==2) then
        allocate(rotation);rotation%axis=2;rotation%angular_speed=1.e10_dp
        rotation%axis_point=[10._dp,10._dp,10._dp]
      endif
      seed=tree
      if(allocated(seed_held%fields)) deallocate(seed_held%fields)
      allocate(seed_held%fields(count))
      do p=1,count
        seed%nodes(p)%temperature=seed%nodes(p)%temperature+12._dp
        allocate(seed_held%fields(p)%values(1,4,4,4),source=0._dp)
        call initialize_nasa7_internal_energy_patch(species,seed%nodes(p)%state,seed%nodes(p)%temperature, &
          seed%nodes(p)%geometry%volume_fraction,seed_held%fields(p)%values(1,:,:,:),good)
        call require(good,'tree held fixture uses independent native species EOS')
        do z=1,4
        do y=1,4
        do x=1,4
          if(seed%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          rho=seed%nodes(p)%state(1,x,y,z)
          seed%nodes(p)%state(2:4,x,y,z)=[rho*1.e10_dp,0._dp,0._dp]
          seed%nodes(p)%state(5,x,y,z)=.5_dp*rho*1.e20_dp
          if(allocated(rotation)) then
            associate(g=>seed%nodes(p)%geometry)
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_potential_3d(rotation,position,phi,good)
            call require(good,'tree held fixture rotating energy')
            seed%nodes(p)%state(5,x,y,z)=seed%nodes(p)%state(5,x,y,z)-rho*phi
          endif
        enddo
        enddo
        enddo
      enddo
      reference=seed;reference_held=seed_held
      call average_down_reactive_eb_tree_3d(species,reference,good,rotation,reference_held)
      call require(good,'native serial tree gas/T/held average survives unresolved thermal subtraction')
      serial_trial=seed;sh=seed_held;sh%fields(2)%values(1,4,4,4)=1.e30_dp
      serial_saved=serial_trial;ss=sh
      call average_down_reactive_eb_tree_3d(species,serial_trial,good,rotation,sh)
      call require(.not.good,'late native held EOS refuses entire serial tree average')
      do p=1,count
        call require(all(serial_trial%nodes(p)%state==serial_saved%nodes(p)%state).and. &
          all(serial_trial%nodes(p)%temperature==serial_saved%nodes(p)%temperature).and. &
          all(sh%fields(p)%values==ss%fields(p)%values),'serial late held average restores all completed descendants')
      enddo
      do layout=1,2
        owner=[(mod(p-1,nranks),p=1,count)];if(layout==2) owner=0
        allocate(gas%fields(count),temp%fields(count),held%fields(count))
        do p=1,count
          if(owner(p)/=rank) cycle
          gas%fields(p)%values=seed%nodes(p)%state
          allocate(temp%fields(p)%values(1,4,4,4))
          temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          held%fields(p)=seed_held%fields(p)
        enddo
        call average_down_mpi_eb_tree_3d(species,topology,owner,geometry,gas,temp,MPI_COMM_WORLD,good,rotation,held)
        call require(good,'native MPI tree gas/T/held restriction on distributed/owner-zero layouts')
        do p=1,count
          if(owner(p)/=rank) cycle
          call require(all(gas%fields(p)%values==reference%nodes(p)%state).and. &
            all(temp%fields(p)%values(1,:,:,:)==reference%nodes(p)%temperature).and. &
            all(held%fields(p)%values==reference_held%fields(p)%values), &
            'native serial/MPI tree gas/T/held average agrees exactly')
        enddo
        do p=1,count
          if(owner(p)/=rank) cycle
          gas%fields(p)%values=seed%nodes(p)%state;temp%fields(p)%values(1,:,:,:)=seed%nodes(p)%temperature
          held%fields(p)=seed_held%fields(p)
          if(p==2) held%fields(p)%values(1,4,4,4)=1.e30_dp
        enddo
        g0=gas;t0=temp;h0=held
        call average_down_mpi_eb_tree_3d(species,topology,owner,geometry,gas,temp,MPI_COMM_WORLD,good,rotation,held)
        call require(.not.good,'remote late held EOS refuses MPI tree average collectively')
        do p=1,count
          if(owner(p)/=rank) cycle
          call require(all(gas%fields(p)%values==g0%fields(p)%values).and. &
            all(temp%fields(p)%values==t0%fields(p)%values).and.all(held%fields(p)%values==h0%fields(p)%values), &
            'MPI failed held average preserves all prior gas/T/scalar fields')
        enddo
        if(nranks>1) then
          if(rank==nranks-1) then
            call average_down_mpi_eb_tree_3d(species,topology,owner,geometry,gas,temp,MPI_COMM_WORLD,good,rotation)
          else
            call average_down_mpi_eb_tree_3d(species,topology,owner,geometry,gas,temp,MPI_COMM_WORLD,good,rotation,held)
          endif
          call require(.not.good,'rank-disagreeing held tree average presence refuses collectively')
        do p=1,count
          if(owner(p)/=rank) cycle
          call require(all(gas%fields(p)%values==g0%fields(p)%values).and. &
            all(temp%fields(p)%values==t0%fields(p)%values).and.all(held%fields(p)%values==h0%fields(p)%values), &
            'MPI failed held average preserves all prior gas/T/scalar fields')
        enddo
        endif
        deallocate(gas%fields,temp%fields,held%fields)
      enddo
    enddo
  end subroutine check_held_tree_average
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: ierr
    if (condition) return
    print *, 'FAIL rank ',rank,': ',label
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  end subroutine
  subroutine unchanged()
    integer :: n
    do n=1,5
      if (owners(n)/=rank) cycle
      call require(all(local%fields(n)%values==saved%fields(n)%values),'state rollback')
      call require(all(lt%fields(n)%values==st%fields(n)%values),'temperature rollback')
    end do
  end subroutine
end program
