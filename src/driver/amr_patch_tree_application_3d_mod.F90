module amr_patch_tree_application_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use amr_patch_tree_remap_3d_mod, only: remap_amr_patch_tree_reactive_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use mixture_thermo_mod, only: mixture_density,mass_fractions_from_mole_fractions
  use selected_composition_mod, only: resolve_selected_composition
  use reactive_1d_mod, only: reactive_nvar,reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_boundaries_3d
  use reactive_entropy_wave_3d_problem_mod, only: initialize_reactive_problem_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d,average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d
  use amr_patch_tree_transport_3d_mod, only: compute_amr_patch_tree_transport_timestep_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d
  use reactive_csv_io_3d_mod, only: reactive_3d_csv_header,format_reactive_3d_csv_row
  implicit none
  private
  public :: run_amr_patch_tree_application_3d
contains
  subroutine run_amr_patch_tree_application_3d(topology_path,config,species,reactions,transport,mole_fractions,schedule_path)
    character(len=*), intent(in) :: topology_path
    character(len=*), optional, intent(in) :: schedule_path
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(reactive_3d_config) :: local
    type(amr_patch_tree_3d) :: tree
    type(amr_patch_tree_3d), allocatable :: events(:)
    type(amr_tree_state_3d) :: state,temperature
    integer, allocatable :: ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    real(dp), allocatable :: x(:),y(:),z(:),root_h(:),h(:,:),recovered(:,:,:),mass_fractions(:),inlet_x(:)
    real(dp), allocatable :: initial(:),final(:),budget(:),change(:),inlet(:,:),q(:),inlet_cube(:,:,:,:)
    real(dp), allocatable :: event_times(:)
    real(dp) :: stop_time,max_remap_error
    real(dp) :: time,dt,transport_dt,diffusivity,base_density,correction,theta,max_correction,min_theta
    real(dp) :: inlet_t(6),inlet_temperature(2,2,2),inlet_dt,face_dt,sound,min_h(3)
    integer :: count,unit,status,node,n(3),lo(3),hi(3),i,j,k,step,attempt,face,c,child,nevents,event,slash,regrids
    character(len=32) :: faces(6)
    character(len=1024) :: message,event_path
    character(len=:), allocatable :: row
    logical :: ok,covered,viscosity,conduction,diffusion,barodiffusion
    call read_topology(topology_path,tree)
    count=size(tree%nodes)-1
    nevents=0
    if (present(schedule_path)) then
      open(newunit=unit,file=schedule_path,status='old',action='read',iostat=status)
      if (status/=0) error stop 'Cannot open tree regrid schedule'
      read(unit,*,iostat=status) nevents
      if (status/=0) error stop 'Invalid tree regrid event count'
      if (nevents<0) error stop 'Invalid tree regrid event count'
    end if
    allocate(events(nevents),event_times(nevents),stat=status)
    if (status/=0) error stop 'Cannot allocate tree regrid schedule'
    do event=1,nevents
      read(unit,*,iostat=status) event_times(event),event_path
      if (status/=0) error stop 'Invalid tree regrid event'
      if (.not.ieee_is_finite(event_times(event))) error stop 'Nonfinite tree regrid time'
      if (event_times(event)<=0.or.event_times(event)>config%final_time) error stop 'Tree regrid time outside run'
      if (event>1) then
        if (event_times(event)<=event_times(event-1)) error stop 'Tree regrid times must increase'
      end if
      if (len_trim(event_path)==0) error stop 'Empty regrid topology path'
      if (event_path(1:1)/='/') then
        slash=scan(schedule_path,'/',back=.true.)
        event_path=schedule_path(:slash)//trim(event_path)
      end if
      call read_topology(trim(event_path),events(event))
    end do
    if (present(schedule_path)) close(unit)
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) error stop 'Tree coordinate overflow'
    root_h=[(config%x_upper-config%x_lower)/config%nx,(config%y_upper-config%y_lower)/config%ny, &
      (config%z_upper-config%z_lower)/config%nz]
    allocate(h(3,count+1),state%fields(count+1),temperature%fields(count+1))
    allocate(mass_fractions(size(species)),inlet_x(size(species)),q(reactive_nvar(size(species))))
    allocate(initial(size(q)),final(size(q)),budget(size(q)),change(size(q)),inlet(size(q),6),inlet_cube(size(q),2,2,2))
    do node=1,count+1
      h(:,node)=root_h*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
      n=tree%nodes(node)%shape
      call coordinates(node)
      allocate(state%fields(node)%values(size(q),n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)))
      local=config;local%nx=n(1);local%ny=n(2);local%nz=n(3)
      call initialize_reactive_problem_3d(species,local,x,y,z,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),base_density,mass_fractions,ok,mole_fractions)
      if (.not.ok) error stop 'Cannot initialize tree node'
    end do
    call average_down_amr_patch_tree_3d(tree,state,ok)
    if (.not.ok) error stop 'Cannot synchronize initial tree'
    do node=1,count+1
      n=tree%nodes(node)%shape;allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,ok)
      if (.not.ok) error stop 'Invalid initial tree temperature'
      temperature%fields(node)%values(1,:,:,:)=recovered;deallocate(recovered)
    end do
    faces=resolved_reactive_boundaries_3d(config)
    viscosity=config%transport_enabled.and.config%viscosity_enabled
    conduction=config%transport_enabled.and.config%thermal_conduction_enabled
    diffusion=config%transport_enabled.and.config%species_diffusion_enabled
    barodiffusion=config%transport_enabled.and.config%barodiffusion_enabled
    min_h=minval(h,dim=2)
    ! Reservoir bounds remain conservative for every scheduled refinement.
    do event=1,nevents
      call include_event_spacing(events(event))
    end do
    inlet=0;inlet_t=0;inlet_dt=huge(1.0_dp)
    if (any(faces=='inflow').and.config%inflow_composition_count/=0) then
      call resolve_selected_composition('Inflow',config%inflow_composition_count,config%inflow_composition_species, &
        config%inflow_composition_mole_fractions,species,inlet_x,ok,message)
      if (.not.ok) error stop 'Invalid inflow composition'
      call mass_fractions_from_mole_fractions(species,inlet_x,mass_fractions,ok)
      if (.not.ok) error stop 'Invalid inflow mass fractions'
    end if
    do face=1,6
      if (faces(face)/='inflow') cycle
      q=0;q(1)=mixture_density(species,mass_fractions,config%inflow_pressure(face),config%inflow_temperature(face),ok)
      if (.not.ok) error stop 'Invalid inflow density'
      q(2:4)=config%inflow_velocity(:,face);q(5)=config%inflow_pressure(face);q(6:)=mass_fractions
      call reactive_primitive_to_conserved(species,q,inlet(:,face),inlet_t(face),sound,ok)
      if (.not.ok) error stop 'Invalid inflow state'
      inlet_dt=min(inlet_dt,config%cfl/sum((abs(q(2:4))+sound)/min_h))
      if (config%transport_enabled) then
        do k=1,2
          do j=1,2
            do i=1,2
              inlet_cube(:,i,j,k)=inlet(:,face)
            end do
          end do
        end do
        inlet_temperature=inlet_t(face)
        call reactive_transport_timestep_3d(species,transport,inlet_cube,inlet_temperature,2,2,2, &
          min_h(1),min_h(2),min_h(3),config%transport_cfl,viscosity,conduction,diffusion,face_dt,diffusivity,ok)
        if (.not.ok) error stop 'Invalid inflow transport bound'
        inlet_dt=min(inlet_dt,face_dt)
      end if
    end do
    call composite_integrals_amr_patch_tree_3d(tree,state,root_h,initial,ok)
    if (.not.ok) error stop 'Invalid initial composite inventory'
    time=0;step=0;budget=0;max_correction=0;min_theta=1;event=1;regrids=0;max_remap_error=0
    write(*,'(a)') 'PeleF '//pelef_version//' experimental serial regular 3D patch tree'
    write(*,'(a,i0,a,i0)') 'Nodes: ',count+1,' Levels: ',tree%level_count()
    do
      if (event<=nevents) then
        if (time==event_times(event)) then
          call apply_event(events(event))
          regrids=regrids+1;event=event+1
        end if
      end if
      if (time>=config%final_time) exit
      stop_time=config%final_time
      if (event<=nevents) stop_time=event_times(event)
      if (step>=config%maximum_steps) error stop 'Tree maximum_steps reached'
      call compute_amr_patch_tree_cfl_3d(tree,species,state,temperature,root_h,config%cfl,dt,ok)
      if (.not.ok) error stop 'Invalid tree hydro timestep'
      dt=min(dt,inlet_dt,stop_time-time)
      if (config%transport_enabled) then
        call compute_amr_patch_tree_transport_timestep_3d(tree,species,transport,state,temperature,root_h, &
          config%transport_cfl,viscosity,conduction,diffusion,transport_dt,diffusivity,ok,faces,config%wall_temperature)
        if (.not.ok) error stop 'Invalid tree transport timestep'
        dt=min(dt,transport_dt)
      end if
      do attempt=1,20
        if (time+dt<=time) error stop 'Tree timestep cannot advance time'
        call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,state,temperature,root_h,dt,faces, &
          config%reconstruction,config%limiter,config%riemann_solver,config%chemistry_enabled, &
          viscosity,conduction,diffusion,barodiffusion,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          change,correction,theta,ok,inflow_state=inlet,inflow_temperature=inlet_t, &
          wall_temperature=config%wall_temperature,wall_velocity=config%wall_velocity)
        if (ok) exit
        dt=.5_dp*dt
      end do
      if (.not.ok) error stop 'Tree step failed after 20 attempts'
      if (dt==stop_time-time) then
        time=stop_time
      else
        time=time+dt
      end if
      step=step+1;budget=budget+change
      max_correction=max(max_correction,correction);min_theta=min(min_theta,theta)
    end do
    call composite_integrals_amr_patch_tree_3d(tree,state,root_h,final,ok)
    if (.not.ok) error stop 'Invalid final composite inventory'
    open(newunit=unit,file=trim(config%output_file),status='replace',action='write',iostat=status)
    if (status/=0) error stop 'Cannot open tree CSV'
    write(unit,'(a)',iostat=status) 'node,level,dx,dy,dz,'//reactive_3d_csv_header(species)
    if (status/=0) error stop 'Cannot write tree CSV header'
    allocate(character(len=26*(14+2*size(species))-1)::row)
    do node=1,count+1
      n=tree%nodes(node)%shape;call coordinates(node)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            covered=.false.
            do child=node+1,count+1
              if (tree%nodes(child)%parent/=node) cycle
              lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
                tree%nodes(child)%patch%coarse_k_lower]
              hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
                tree%nodes(child)%patch%coarse_k_upper]
              if (all([i,j,k]>=lo).and.all([i,j,k]<=hi)) covered=.true.
            end do
            if (covered) cycle
            call format_reactive_3d_csv_row(species,[x(i),y(j),z(k)],time,state%fields(node)%values(:,i,j,k), &
              temperature%fields(node)%values(1,i,j,k),row,ok)
            if (.not.ok) error stop 'Invalid tree output state'
            write(unit,'(i0,",",i0,3(",",es25.16e3),",",a)',iostat=status) node,tree%nodes(node)%level,h(:,node),row
            if (status/=0) error stop 'Cannot write tree CSV row'
          end do
        end do
      end do
    end do
    close(unit,iostat=status)
    if (status/=0) error stop 'Cannot close tree CSV'
    open(newunit=unit,file=trim(config%output_file)//'.budget.csv',status='replace',action='write',iostat=status)
    if (status/=0) error stop 'Cannot open tree budget CSV'
    write(unit,'(a)',iostat=status) 'component,initial,final,boundary_change'
    if (status/=0) error stop 'Cannot write tree budget header'
    do c=1,size(q)
      write(unit,'(i0,3(",",es25.16e3))',iostat=status) c,initial(c),final(c),budget(c)
      if (status/=0) error stop 'Cannot write tree budget row'
    end do
    close(unit,iostat=status)
    if (status/=0) error stop 'Cannot close tree budget CSV'
    write(*,'(a,i0,a,es12.4)') 'Regrids: ',regrids,' Maximum remap error: ',max_remap_error
    write(*,'(a,i0,a,es24.16)') 'Steps: ',step,' Final time: ',time
    write(*,'(a,es12.4,a,es12.4)') 'Maximum reflux: ',max_correction,' Minimum transport theta: ',min_theta
  contains
    subroutine read_topology(path,result)
      character(len=*), intent(in) :: path
      type(amr_patch_tree_3d), intent(out) :: result
      type(amr_patch_3d), allocatable :: patches(:)
      integer, allocatable :: parents(:),shapes(:,:)
      integer :: input,ios,nchild,peer,parent,dims(3),lower(3),upper(3),ratio
      logical :: valid
      open(newunit=input,file=path,status='old',action='read',iostat=ios)
      if (ios/=0) error stop 'Cannot open tree topology'
      read(input,*,iostat=ios) nchild
      if (ios/=0) error stop 'Invalid tree child count'
      if (nchild<0.or.nchild==huge(nchild)) error stop 'Invalid tree child count'
      allocate(patches(nchild),parents(nchild),shapes(3,nchild+1),stat=ios)
      if (ios/=0) error stop 'Cannot allocate tree topology'
      shapes(:,1)=[config%nx,config%ny,config%nz]
      do peer=1,nchild
        read(input,*,iostat=ios) parent,lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),ratio
        if (ios/=0) error stop 'Invalid tree patch record'
        if (parent<1.or.parent>peer) error stop 'Tree parents must precede children'
        parents(peer)=parent;dims=shapes(:,parent)
        call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lower(1),upper(1),lower(2),upper(2), &
          lower(3),upper(3),ratio,patches(peer),valid)
        if (.not.valid) error stop 'Invalid tree patch geometry'
        if (any(int(upper-lower+1,int64)*int(ratio,int64)>int(huge(1),int64))) error stop 'Tree patch extent overflow'
        shapes(:,peer+1)=[patches(peer)%fine_nx(),patches(peer)%fine_ny(),patches(peer)%fine_nz()]
      end do
      close(input)
      call initialize_amr_patch_tree_3d(shapes(:,1),parents,patches,result,valid)
      if (.not.valid) error stop 'Invalid tree topology'
    end subroutine
    subroutine include_event_spacing(next_tree)
      type(amr_patch_tree_3d), intent(in) :: next_tree
      integer(int64), allocatable :: o(:,:),d(:,:)
      integer, allocatable :: r(:)
      logical :: valid
      call amr_tree_integer_layout_3d(next_tree,o,d,r,valid)
      if (.not.valid) error stop 'Scheduled tree coordinate overflow'
      min_h=min(min_h,root_h*real(d(:,1),dp)/real(d(:,size(r)),dp))
    end subroutine
    subroutine apply_event(next_tree)
      type(amr_patch_tree_3d), intent(in) :: next_tree
      type(amr_tree_state_3d) :: candidate,temps,absolute_state
      integer(int64), allocatable :: o(:,:),d(:,:)
      integer, allocatable :: r(:)
      real(dp), allocatable :: spacing(:,:)
      real(dp) :: before(size(q)),after(size(q)),scale(size(q)),error
      integer :: peer
      logical :: valid
      call remap_amr_patch_tree_reactive_3d(tree,state,temperature,next_tree,species,candidate,temps,valid)
      if (.not.valid) error stop 'Scheduled tree remap failed'
      call composite_integrals_amr_patch_tree_3d(tree,state,root_h,before,valid)
      if (.not.valid) error stop 'Invalid pre-remap inventory'
      call composite_integrals_amr_patch_tree_3d(next_tree,candidate,root_h,after,valid)
      if (.not.valid) error stop 'Invalid post-remap inventory'
      absolute_state=state
      do peer=1,size(tree%nodes)
        absolute_state%fields(peer)%values=abs(state%fields(peer)%values)
      end do
      call composite_integrals_amr_patch_tree_3d(tree,absolute_state,root_h,scale,valid)
      if (.not.valid) error stop 'Invalid remap inventory scale'
      error=maxval(abs(after-before)/max(tiny(1.0_dp),scale))
      if (error>5e-12_dp) error stop 'Scheduled tree remap lost conservation'
      call amr_tree_integer_layout_3d(next_tree,o,d,r,valid)
      if (.not.valid) error stop 'Scheduled tree coordinate overflow'
      allocate(spacing(3,size(next_tree%nodes)))
      do peer=1,size(next_tree%nodes)
        spacing(:,peer)=root_h*real(d(:,1),dp)/real(d(:,next_tree%nodes(peer)%level+1),dp)
      end do
      tree=next_tree;state=candidate;temperature=temps
      origin=o;domain=d;ratios=r;h=spacing;count=size(tree%nodes)-1
      max_remap_error=max(max_remap_error,error)
      write(*,'(a,es24.16,a,i0,a,i0)') 'Regrid time: ',time,' Nodes: ',count+1,' Levels: ',tree%level_count()
    end subroutine
    subroutine coordinates(peer)
      integer, intent(in) :: peer
      integer :: a,dims(3)
      dims=tree%nodes(peer)%shape
      if (allocated(x)) deallocate(x,y,z)
      allocate(x(dims(1)),y(dims(2)),z(dims(3)))
      x=[(config%x_lower+(real(origin(1,peer),dp)+real(a,dp)-.5_dp)*h(1,peer),a=1,dims(1))]
      y=[(config%y_lower+(real(origin(2,peer),dp)+real(a,dp)-.5_dp)*h(2,peer),a=1,dims(2))]
      z=[(config%z_lower+(real(origin(3,peer),dp)+real(a,dp)-.5_dp)*h(3,peer),a=1,dims(3))]
    end subroutine
  end subroutine
end module amr_patch_tree_application_3d_mod
