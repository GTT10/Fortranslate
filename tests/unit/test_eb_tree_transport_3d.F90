program test_eb_tree_transport_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    reactive_eb_tree_timestep_3d, &
    average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d,advance_reactive_eb_tree_transport_3d, &
    compute_sdc_eb_tree_transport_3d,synchronize_reactive_eb_tree_flux_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(amr_patch_3d) :: boxes(4)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,middle,coarse
  type(eb_reactive_tree_3d) :: tree,saved
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: u(:,:,:,:),t(:,:,:)
  real(dp), allocatable :: before(:),after(:),boundary(:),primitive(:)
  real(dp) :: residual,theta,speed
  character(len=1), parameter :: axes(3)=['x','y','z']
  integer :: nv,node,axis,count,i,j,k
  character(len=16) :: periodic_faces(6)
  real(dp) :: wall_t(6),base_dt,thermal_dt
  integer :: bound_mode
  logical :: tightened=.false.
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(before(nv),after(nv),boundary(nv),primitive(reactive_nprim(size(species))))
  do axis=1,3
  do count=4,5
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
  call build_axis_plane_eb_geometry_3d(16,16,16,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.31_dp,finest,ok)
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
  do node=2,count
    tree%nodes(node)%state=real(node,dp)*tree%nodes(node)%state
    tree%nodes(node)%state(5,:,:,:)=1.1_dp*tree%nodes(node)%state(5,:,:,:)
    do k=1,4
    do j=1,4
    do i=1,4
      if (tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0._dp) cycle
      call reactive_conserved_to_primitive(species,tree%nodes(node)%state(:,i,j,k),t(i,j,k),primitive, &
        tree%nodes(node)%temperature(i,j,k),speed,ok)
      call require(ok,'temperature perturbation EOS')
    end do
    end do
    end do
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok)
  call require(ok,'initial average down')
  call integrate_reactive_eb_tree_3d(tree,before,ok)
  call require(ok,'initial inventory')
  saved=tree
  boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(ok,'tree transport step')
  call require(all(boundary==0._dp),'closed physical diffusion faces')
  call require(theta>=0._dp.and.theta<=1._dp,'limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'final inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'axis, nodes, inventory residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'transport conservation')
  call require(any(tree%nodes(count)%state/=saved%nodes(count)%state),'nontrivial deep transport update')
  if(axis==3.and.count==5) call check_sdc_transport()
  periodic_faces='no_slip_wall'
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(ok,'wall transport update')
  call require(boundary(1)==0.and.boundary(5)==0.and.all(boundary(6:)==0),'impermeable adiabatic wall flux')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'wall inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'wall axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'wall composite transport conservation')
  periodic_faces='isothermal_wall'
  wall_t=0
  do node=1,count
    wall_t=max(wall_t,maxval(saved%nodes(node)%temperature)+100)
  end do
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces,wall_temperature=wall_t)
  call require(ok,'isothermal transport update')
  call require(boundary(1)==0.and.boundary(5)>0.and.all(boundary(6:)==0),'impermeable heated wall flux')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'isothermal inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'isothermal axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'isothermal composite transport conservation')
  do bound_mode=0,1
    call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
      .false.,.true.,.false.,base_dt,ok,subcycling=bound_mode==1)
    call require(ok,'baseline tree timestep')
    call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
      .false.,.true.,.false.,thermal_dt,ok,subcycling=bound_mode==1, &
      boundary_faces=periodic_faces,wall_temperature=wall_t)
    call require(ok.and.thermal_dt>0.and.thermal_dt<=base_dt,'thermal timestep restriction')
    tightened=tightened.or.thermal_dt<base_dt
    print *, 'thermal timestep axis, nodes, mode, baseline, wall: ',axes(axis),count,bound_mode,base_dt,thermal_dt
    thermal_dt=17
    call reactive_eb_tree_timestep_3d(species,transport,saved,.5_dp,1.e-12_dp,.true., &
      .false.,.true.,.false.,thermal_dt,ok,subcycling=bound_mode==1,boundary_faces=periodic_faces)
    call require(.not.ok.and.thermal_dt==17,'missing wall temperature preserves timestep')
  end do

  periodic_faces='periodic';periodic_faces(2*axis-1:2*axis)='outflow'
  tree=saved
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(ok,'periodic transport update')
  call require(theta>=0.and.theta<=1,'periodic limiter range')
  call integrate_reactive_eb_tree_3d(tree,after,ok)
  call require(ok,'periodic inventory')
  residual=maxval(abs(after-before-boundary)/max(1._dp,abs(before)))
  print *, 'periodic axis, nodes, residual: ',axes(axis),count,residual
  call require(residual<3.e-14_dp,'periodic composite transport conservation')
  tree=saved;boundary=17;theta=19
  periodic_faces(2*mod(axis,3)+1)='outflow'
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok,boundary_faces=periodic_faces)
  call require(.not.ok.and.all(boundary==17).and.theta==19,'unpaired periodic transport rollback')
  call unchanged_tree()
  tree=saved;boundary=17._dp;theta=19._dp
  call advance_reactive_eb_tree_transport_3d(species,transport,-1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(.not.ok,'invalid interval rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.theta==19._dp,'diagnostic rollback')
  tree%nodes(count)%state(1,4,4,4)=-1._dp;saved=tree
  call advance_reactive_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
    tree,boundary,theta,ok)
  call require(.not.ok,'invalid deep state rejected')
  call unchanged_tree()
  call require(all(boundary==17._dp).and.theta==19._dp,'deep failure diagnostic rollback')
  deallocate(u,t)
  end do
  end do
  call require(tightened,'wall heat restricts an accepted tree timestep')
  print *, 'test_eb_tree_transport_3d: PASS'
contains
  subroutine check_sdc_transport()
    type(amr_tree_state_3d) :: rates,old_rates
    type(eb_integrated_flux_3d), allocatable :: flux(:),old_flux(:)
    type(eb_reactive_tree_3d) :: result,bad
    real(dp) :: limiter
    type(amr_tree_state_3d) :: level_rates,reference_rates
    type(eb_integrated_flux_3d), allocatable :: level_flux(:),reference_flux(:)
    type(eb_reactive_tree_3d) :: endpoint,interpolated
    integer :: p
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      saved,rates,flux,limiter,ok)
    call require(ok,'SDC provisional transport')
    call require(limiter==theta,'SDC transport limiter matches Euler')
    result=saved
    do p=1,count
      result%nodes(p)%state=result%nodes(p)%state+1.e-8_dp*rates%fields(p)%values
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,result,flux,ok)
    call require(ok,'SDC transport reconstructed synchronization')
    do p=1,count
      call require(all(abs(result%nodes(p)%state-tree%nodes(p)%state)<= &
        2.e-12_dp*max(1._dp,abs(tree%nodes(p)%state))),'SDC transport reconstructs Euler state')
    end do
    endpoint=saved;interpolated=saved
    do p=1,count
      if(saved%topology%nodes(p)%level>=2) cycle
      endpoint%nodes(p)%state=1.01_dp*saved%nodes(p)%state
      interpolated%nodes(p)%state=.75_dp*saved%nodes(p)%state+.25_dp*endpoint%nodes(p)%state
      interpolated%nodes(p)%temperature=.75_dp*saved%nodes(p)%temperature+.25_dp*endpoint%nodes(p)%temperature
    end do
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      interpolated,reference_rates,reference_flux,limiter,ok,target_depth=2)
    call require(ok,'explicit ancestor transport reference')
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      saved,level_rates,level_flux,limiter,ok,target_depth=2, &
      ancestor_start=saved,ancestor_end=endpoint,ancestor_fraction=[.25_dp,.25_dp])
    call require(ok,'selected level transport with ancestor time interpolation')
    do p=1,count
      if(saved%topology%nodes(p)%level==2) then
        call require(all(level_rates%fields(p)%values==reference_rates%fields(p)%values).and. &
          all(level_flux(p)%x==reference_flux(p)%x).and.all(level_flux(p)%y==reference_flux(p)%y).and. &
          all(level_flux(p)%z==reference_flux(p)%z),'exact interpolated level transport')
      else
        call require(all(level_rates%fields(p)%values==0).and..not.allocated(level_flux(p)%x), &
          'inactive level has zero rates and absent flux')
      end if
    end do
    old_rates=rates;old_flux=flux;bad=saved
    bad%nodes(count)%state(1,4,4,4)=-1
    limiter=17
    call compute_sdc_eb_tree_transport_3d(species,transport,1.e-8_dp,.true.,.true.,.true.,.true., &
      bad,rates,flux,limiter,ok)
    call require(.not.ok.and.limiter==17,'SDC transport invalid fluid density')
    do p=1,count
      call require(all(rates%fields(p)%values==old_rates%fields(p)%values).and. &
        all(flux(p)%x==old_flux(p)%x).and.all(flux(p)%y==old_flux(p)%y).and. &
        all(flux(p)%z==old_flux(p)%z),'SDC transport complete output rollback')
    end do
  end subroutine check_sdc_transport

  subroutine unchanged_tree()
    integer :: n
    do n=1,size(tree%nodes)
      call require(all(tree%nodes(n)%state==saved%nodes(n)%state).and. &
        all(tree%nodes(n)%temperature==saved%nodes(n)%temperature),'all state and temperature rollback')
    end do
  end subroutine unchanged_tree
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
  end subroutine require
end program test_eb_tree_transport_3d
