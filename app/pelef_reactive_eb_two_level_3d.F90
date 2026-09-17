program pelef_reactive_eb_two_level_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  use eb_reactive_two_level_3d_mod, only: advance_reactive_eb_two_level_full_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(reactive_eb_3d_config) :: config
  type(eb_geometry_3d) :: full,coarse,fine
  integer :: nx=8,ny=8,nz=8,refinement_ratio=2,fine_lower(3)=[2,2,2],fine_upper(3)=[3,3,3],maximum_steps=100
  character(len=16) :: geometry_type='plane',plane_axis='z',extrusion_axis='z',riemann_solver='rusanov'
  character(len=16) :: chemistry_integrator='explicit'
  character(len=512) :: level_set_file='',section_file='',output_file='eb_two_level.csv'
  real(dp) :: domain_lower(3)=0._dp,domain_upper(3)=1._dp,plane_position=.31_dp,time_step=1.e-8_dp,final_time=3.e-8_dp
  real(dp) :: rtol=2.e-7_dp,atol=1.e-12_dp
  logical :: chemistry_enabled=.true.,transport_enabled=.true.,viscosity_enabled=.true.
  logical :: thermal_conduction_enabled=.true.,species_diffusion_enabled=.true.,barodiffusion_enabled=.true.
  namelist /eb_two_level/ nx,ny,nz,refinement_ratio,fine_lower,fine_upper,maximum_steps,geometry_type,plane_axis, &
    extrusion_axis,level_set_file,section_file,output_file,domain_lower,domain_upper,plane_position,time_step,final_time, &
    rtol,atol,chemistry_enabled,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
    species_diffusion_enabled,barodiffusion_enabled,chemistry_integrator,riemann_solver
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),uf(:,:,:,:),tf(:,:,:),initial(:),current(:),scale(:), &
    boundary(:),wall(:),chemical(:),total_boundary(:),total_wall(:),total_chemical(:)
  logical, allocatable :: refined(:,:,:)
  real(dp) :: time,dt,theta,error,minimum_theta
  integer :: unit,status,close_status,nvar,step,dims(3),nf(3),c
  character(len=1024) :: input_path,message,budget_path
  logical :: ok
  if (command_argument_count()/=1) error stop 'Usage: pelef_reactive_eb_two_level_3d <input.nml>'
  call get_command_argument(1,input_path)
  open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=status)
  if (status/=0) error stop 'Cannot open two-level EB input'
  message=""
  read(unit,nml=eb_two_level,iostat=status,iomsg=message)
  close(unit,iostat=close_status)
  if (status/=0.or.close_status/=0) then
    write(*,'(a)') trim(message)
    error stop 'Invalid two-level EB input'
  end if
  if (min(nx,ny,nz)<2.or.refinement_ratio<2.or.maximum_steps<1) error stop 'Invalid grid or step count'
  if (any(mod([nx,ny,nz],refinement_ratio)/=0)) error stop 'Finest dimensions must divide by refinement ratio'
  dims=[nx,ny,nz]/refinement_ratio
  if (any(fine_lower<1).or.any(fine_upper<fine_lower).or.any(fine_upper>dims)) error stop 'Invalid fine patch box'
  if (.not.all(ieee_is_finite([domain_lower,domain_upper,time_step,final_time,rtol,atol]))) &
    error stop 'Nonfinite two-level controls'
  if (any(domain_upper<=domain_lower).or.min(time_step,final_time,rtol,atol)<=0._dp) &
    error stop 'Invalid two-level physical controls'
  if (barodiffusion_enabled.and..not.species_diffusion_enabled.and.transport_enabled) &
    error stop 'Barodiffusion requires species diffusion'
  select case(trim(riemann_solver))
  case('rusanov','hllc','pelec')
  case default;error stop 'Invalid Riemann solver'
  end select
  if (chemistry_integrator/='explicit'.and.chemistry_integrator/='implicit') error stop 'Invalid chemistry integrator'
  if (len_trim(output_file)==0) error stop 'Empty output path'
  budget_path=trim(output_file)//'.budget.csv'
  if (output_file==input_path.or.budget_path==input_path.or.output_file==level_set_file.or. &
    budget_path==level_set_file.or.output_file==section_file.or.budget_path==section_file) &
    error stop 'Output aliases an input file'
  config%nx=nx;config%ny=ny;config%nz=nz
  config%x_lower=domain_lower(1);config%y_lower=domain_lower(2);config%z_lower=domain_lower(3)
  config%x_upper=domain_upper(1);config%y_upper=domain_upper(2);config%z_upper=domain_upper(3)
  config%geometry_type=geometry_type;config%plane_axis=plane_axis;config%plane_position=plane_position
  config%extrusion_axis=extrusion_axis;config%section_file=section_file;config%level_set_file=level_set_file
  call build_configured_eb_geometry_3d(config,full,ok)
  if (.not.ok) error stop 'Cannot construct finest EB geometry'
  call coarsen_eb_geometry_3d(full,refinement_ratio,coarse,ok)
  if (.not.ok) error stop 'Cannot represent coarse EB geometry'
  call extract_eb_geometry_3d(full,(fine_lower-1)*refinement_ratio+1,fine_upper*refinement_ratio,fine,ok)
  if (.not.ok) error stop 'Cannot extract fine EB patch'
  call load_h2o2_elementary_thermo(species,ok)
  if (.not.ok) error stop 'Thermodynamic load failed'
  call load_h2o2_elementary_mechanism(reactions,ok)
  if (.not.ok) error stop 'Mechanism load failed'
  call load_h2o2_elementary_transport(transport,ok)
  if (.not.ok) error stop 'Transport load failed'
  nvar=reactive_nvar(size(species));nf=[fine%nx,fine%ny,fine%nz]
  allocate(uc(nvar,dims(1),dims(2),dims(3)),tc(dims(1),dims(2),dims(3)), &
    uf(nvar,nf(1),nf(2),nf(3)),tf(nf(1),nf(2),nf(3)),refined(dims(1),dims(2),dims(3)))
  config%nx=dims(1);config%ny=dims(2);config%nz=dims(3)
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,ok)
  if (.not.ok) error stop 'Coarse initialization failed'
  config%nx=nf(1);config%ny=nf(2);config%nz=nf(3)
  call initialize_reactive_eb_density_sheet_3d(species,config,fine,uf,tf,ok)
  if (.not.ok) error stop 'Fine initialization failed'
  refined=.false.;refined(fine_lower(1):fine_upper(1),fine_lower(2):fine_upper(2),fine_lower(3):fine_upper(3))=.true.
  allocate(initial(nvar),current(nvar),scale(nvar),boundary(nvar),wall(nvar),chemical(nvar), &
    total_boundary(nvar),total_wall(nvar),total_chemical(nvar))
  call inventory(initial)
  do c=1,nvar
    scale(c)=max(1.e-30_dp,coarse%dx*coarse%dy*coarse%dz*sum(abs(uc(c,:,:,:))*coarse%volume_fraction, &
      mask=.not.refined)+fine%dx*fine%dy*fine%dz*sum(abs(uf(c,:,:,:))*fine%volume_fraction))
  end do
  total_boundary=0._dp;total_wall=0._dp;total_chemical=0._dp
  time=0._dp;step=0;minimum_theta=1._dp
  do while(time<final_time)
    if (step>=maximum_steps) error stop 'Maximum steps reached before final time'
    dt=min(time_step,final_time-time)
    if (time+dt<=time) error stop 'Time step cannot advance time'
    call advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,fine_lower,refinement_ratio, &
      riemann_solver,dt,chemistry_enabled,rtol,atol,transport_enabled,viscosity_enabled,thermal_conduction_enabled, &
      species_diffusion_enabled,barodiffusion_enabled,uc,tc,uf,tf,boundary,wall,chemical,theta,ok,chemistry_integrator)
    if (.not.ok) error stop 'Two-level EB step rejected; reduce time_step or inspect input'
    total_boundary=total_boundary+boundary;total_wall=total_wall+wall;total_chemical=total_chemical+chemical
    time=time+dt;step=step+1;minimum_theta=min(minimum_theta,theta)
  end do
  call inventory(current)
  ! Use initial L1, final and exchange scales with an absolute floor of one.
  error=maxval(abs(current-initial-total_boundary-total_wall-total_chemical)/ &
    max(1._dp,scale,abs(current),abs(total_boundary)+abs(total_wall)+abs(total_chemical)))
  if (.not.ieee_is_finite(error)) error stop 'Nonfinite composite inventory residual'
  if (error>5.e-10_dp) error stop 'Composite inventory balance failed'
  call write_output()
  write(*,'(a,i0)') 'Completed steps: ',step
  write(*,'(a,es24.16)') 'Final time: ',time
  write(*,'(a,es24.16)') 'Composite inventory error: ',error
  write(*,'(a,es24.16)') 'Minimum transport theta: ',minimum_theta
contains
  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: component
    do component=1,size(values)
      values(component)=coarse%dx*coarse%dy*coarse%dz*sum(uc(component,:,:,:)*coarse%volume_fraction,mask=.not.refined) &
        +fine%dx*fine%dy*fine%dz*sum(uf(component,:,:,:)*fine%volume_fraction)
    end do
  end subroutine inventory

  subroutine write_output()
    integer :: stream,ios,component,level,i,j,k
    open(newunit=stream,file=trim(output_file),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open composite CSV'
    write(stream,'(a)',advance='no',iostat=ios) 'level,i,j,k,x,y,z,kappa,temperature,rho,rhou,rhov,rhow,rhoE'
    if (ios/=0) error stop 'CSV header failure'
    do component=1,size(species)
      write(stream,'(a)',advance='no',iostat=ios) ',rhoY_'//trim(species(component)%name)
      if (ios/=0) error stop 'CSV species header failure'
    end do
    write(stream,*,iostat=ios)
    if (ios/=0) error stop 'CSV header terminator failure'
    do level=0,1
      if (level==0) then
        do k=1,coarse%nz
          do j=1,coarse%ny
            do i=1,coarse%nx
              if (refined(i,j,k).or.coarse%volume_fraction(i,j,k)<=0._dp) cycle
              call row(stream,level,i,j,k,coarse,uc(:,i,j,k),tc(i,j,k))
            end do
          end do
        end do
      else
        do k=1,fine%nz
          do j=1,fine%ny
            do i=1,fine%nx
              if (fine%volume_fraction(i,j,k)<=0._dp) cycle
              call row(stream,level,i,j,k,fine,uf(:,i,j,k),tf(i,j,k))
            end do
          end do
        end do
      end if
    end do
    close(stream,iostat=ios)
    if (ios/=0) error stop 'CSV close failure'
    open(newunit=stream,file=trim(budget_path),status='replace',action='write',iostat=ios)
    if (ios/=0) error stop 'Cannot open budget CSV'
    write(stream,'(a)',iostat=ios) 'component,initial,final,boundary,wall,chemistry,residual'
    if (ios/=0) error stop 'Budget header failure'
    do component=1,nvar
      write(stream,'(i0,6(",",es26.17e3))',iostat=ios) component,initial(component),current(component), &
        total_boundary(component),total_wall(component),total_chemical(component), &
        current(component)-initial(component)-total_boundary(component)-total_wall(component)-total_chemical(component)
      if (ios/=0) error stop 'Budget row failure'
    end do
    close(stream,iostat=ios)
    if (ios/=0) error stop 'Budget close failure'
  end subroutine write_output

    subroutine row(stream,level,i,j,k,geometry,state,temperature)
      integer, intent(in) :: stream,level,i,j,k
      integer :: ios
      type(eb_geometry_3d), intent(in) :: geometry
      real(dp), intent(in) :: state(:),temperature
      write(stream,'(i0,3(",",i0),*(:,",",es26.17e3))',iostat=ios) level,i,j,k, &
        geometry%x_lower+(real(i,dp)-.5_dp+geometry%cell_centroid_x(i,j,k))*geometry%dx, &
        geometry%y_lower+(real(j,dp)-.5_dp+geometry%cell_centroid_y(i,j,k))*geometry%dy, &
        geometry%z_lower+(real(k,dp)-.5_dp+geometry%cell_centroid_z(i,j,k))*geometry%dz, &
        geometry%volume_fraction(i,j,k),temperature,state
      if (ios/=0) error stop 'Composite CSV row failure'
    end subroutine row
end program pelef_reactive_eb_two_level_3d
