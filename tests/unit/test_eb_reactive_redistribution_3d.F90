program test_eb_reactive_redistribution_3d
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_value, ieee_quiet_nan,ieee_is_nan
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use eb_geometry_2d_mod, only: eb_geometry_2d,build_eb_geometry_2d
  use state_indices_mod, only: irho
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_mass_fraction_component, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: &
    eb_geometry_3d, eb_covered_cell_3d, eb_cut_cell_3d, &
    build_axis_plane_eb_geometry_3d,extrude_eb_geometry_3d,build_level_set_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_redistribution_3d_mod, only: &
    compute_state_redist_level_stage_3d,reactive_eb_flux_redistribute_3d,redistribute_sdc_eb_rate_3d, &
    advance_reactive_eb_redistributed_3d, &
    reactive_eb_weighted_state_redistribute_3d, &
    advance_reactive_eb_state_redistributed_3d,state_redistribution_graph_slopes_3d, &
    state_redist_box_3d,state_redist_register_3d,initialize_state_redist_register_3d,accumulate_state_redist_transfer_3d, &
    srd_coarse_regular,srd_coarse_boundary,srd_coarse_fine_covered, &
    restrict_state_redist_fine_delta_3d,build_state_redist_coarse_flags_3d,initialize_state_redist_tree_register_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d,initialize_amr_patch_tree_3d
  use eb_reactive_redistribution_3d_mod, only: accumulate_state_redist_register_3d, &
    accumulate_state_redist_coarse_register_3d,assemble_state_redist_interface_source_3d, &
    assemble_state_redist_interface_rereflux_3d
  use eb_geometry_3d_mod, only: eb_cell_connectivity_3d
  use eb_reactive_redistribution_3d_mod, only: post_eb_spatial_density_rate_3d
  implicit none
  character(len=32) :: test_mode

  call get_command_argument(1,test_mode)
  call check_post_density_rate()
  if(test_mode=='density_rate') then
    write(*,'(a)') 'EB post-density covered/source/rollback checks: PASS'
    stop
  endif
  call check_connectivity_bits()
  call check_interface_source()
  call check_register_accumulation()
  call check_register_restriction()
  call check_multilevel_registers()
  call check_register_valid_box()
  call check_tree_register_masks()
  call check_update_scales()
  call check_level_stages()
  call check_patch_boxes()
  call check_periodic_state_halos()
  call check_boundary_state_halos()
  call check_boundary_graph_stencil()
  call check_grown_graph_stencil()
  call check_high_order_state()
  call check_sdc_rates()
  call check_regular_and_covered()
  call check_axis("x", 8, 1, 1, [3, 1, 1], [4, 1, 1])
  call check_axis("y", 1, 8, 1, [1, 3, 1], [1, 4, 1])
  call check_axis("z", 1, 1, 8, [1, 1, 3], [1, 1, 4])
  call check_state_axis("x", 8, 1, 1, [3, 1, 1], [4, 1, 1])
  call check_state_axis("y", 1, 8, 1, [1, 3, 1], [1, 4, 1])
  call check_state_axis("z", 1, 1, 8, [1, 1, 3], [1, 1, 4])
  call check_transverse_neighbors()
  call check_reactive_advance()
  call check_curved_state()
  call check_connected_clusters()
  call check_weighted_boundary_pair()

  write(*, '(a)') "test_eb_reactive_redistribution_3d: PASS"

contains
  subroutine check_post_density_rate()
    type(eb_geometry_3d) :: geometry
    real(dp) :: rate(8,4,1,1),kept(8,4,1,1)
    logical :: good
    call build_axis_plane_eb_geometry_3d(4,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',.3_dp,geometry,good)
    call require(good,'post-density covered/cut/regular geometry')
    rate=7
    rate(1,:,:,:)=97
    rate(6,:,:,:)=.25_dp;rate(7,:,:,:)=-.125_dp;rate(8,:,:,:)=.5_dp
    kept=rate
    ! Covered garbage is discarded before the fluid-state finite check.
    rate(:,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call post_eb_spatial_density_rate_3d(geometry,rate,good)
    call require(good,'post-density native ordered species sum')
    call require(all(rate(:,1,1,1)==0),'post-density all covered sources zero')
    call require(all(rate(1,2:4,1,1)==.625_dp),'post-density density equals signed species sum')
    call require(all(rate(2:8,2:4,1,1)==kept(2:8,2:4,1,1)), &
      'post-density species/momentum/energy untouched')
    kept=rate;rate(:,1,1,1)=13;rate(1,2:4,1,1)=97
    rate(5,4,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call post_eb_spatial_density_rate_3d(geometry,rate,good)
    call require(.not.good.and.all(rate(:,1,1,1)==13).and.all(rate(1,2:4,1,1)==97), &
      'post-density late nonfinite source retains every prior cell')
    call require(ieee_is_nan(rate(5,4,1,1)),'post-density failure retains invalid caller data')
    rate=kept;rate(6:8,4,1,1)=huge(1._dp)
    call post_eb_spatial_density_rate_3d(geometry,rate,good)
    call require(.not.good.and.all(rate(1,:,:,:)==kept(1,:,:,:)), &
      'post-density overflowing sum rejected before publication')
  end subroutine check_post_density_rate
  subroutine check_connectivity_bits()
    type(eb_geometry_3d) :: geometry
    logical :: connected(3,3,3),kept(3,3,3),good
    call build_axis_plane_eb_geometry_3d(5,5,5,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'x',-1._dp,geometry,good)
    call require(good,'connectivity regular geometry')
    connected=.false.
    call eb_cell_connectivity_3d(geometry,[3,3,3],connected,good)
    call require(good.and.all(connected),'open apertures set all 27 connectivity bits')
    geometry%x_face_fraction=0;geometry%y_face_fraction=0;geometry%z_face_fraction=0
    geometry%y_face_fraction(3,3,3)=1
    geometry%x_face_fraction(3,4,3)=1
    geometry%y_face_fraction(4,3,3)=1
    geometry%z_face_fraction(4,4,3)=1
    call eb_cell_connectivity_3d(geometry,[3,3,3],connected,good)
    call require(good.and.connected(3,3,2).and.connected(3,3,3), &
      'two-face edge and three-face corner paths set connection bits')
    call require(.not.connected(3,2,2),'three-face detour does not set a closed axial bit')
    call require(.not.connected(1,1,1),'isolated opposite corner remains disconnected')
    kept=connected
    geometry%y_face_fraction(3,3,3)=ieee_value(0._dp,ieee_quiet_nan)
    call eb_cell_connectivity_3d(geometry,[3,3,3],connected,good)
    call require(.not.good.and.all(connected.eqv.kept),'invalid aperture preserves connectivity output')
    geometry%y_face_fraction(3,3,3)=1
    geometry%cell_type(3,3,3)=0
    call eb_cell_connectivity_3d(geometry,[3,3,3],connected,good)
    call require(good.and..not.any(connected),'covered root has no connection bits')
  end subroutine check_connectivity_bits
  subroutine check_interface_source()
    type(amr_patch_tree_3d) :: tree
    type(amr_patch_3d) :: patches(1)
    type(eb_geometry_3d) :: geometry(2)
    type(state_redist_register_3d) :: registers(2)
    type(amr_tree_state_3d) :: source,kept,increment,kept_increment
    character(len=16) :: faces(6)
    real(dp) :: denominator,expected,mass
    integer :: mode,p
    logical :: good
    call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patches(1),good)
    call require(good,'interface source child bounds')
    call initialize_amr_patch_tree_3d([4,4,4],[1],patches,tree,good)
    call require(good,'interface source topology')
    call build_axis_plane_eb_geometry_3d(16,16,16,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
      'y',.1_dp,geometry(1),good)
    call require(good,'interface source coarse prepared geometry')
    call build_axis_plane_eb_geometry_3d(16,16,16,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp, &
      'y',.1_dp,geometry(2),good)
    call require(good,'interface source fine prepared geometry')
    do mode=0,1
      faces='outflow'
      if(mode==1) faces(1:2)='periodic'
      do p=1,2
        call initialize_state_redist_tree_register_3d(tree,p,1,6,registers(p),good,faces,.5_dp)
        call require(good,'interface source prepared registers')
      enddo
      registers(1)%coarse_delta(1,8,7,8)=.25_dp
      registers(1)%coarse_delta(1,7,8,8)=.5_dp
      registers(2)%fine_delta(1,7:8,6,7:8)=2
      ! Outside the fixed dm valid+1 box: contributes neither numerator nor
      ! clipping to the full-block geometric denominator.
      registers(2)%fine_delta(1,7:8,5,7:8)=100
      denominator=sum(geometry(2)%volume_fraction(7:8,5:6,7:8))
      expected=.25_dp+8._dp/denominator
      call assemble_state_redist_interface_source_3d(tree,geometry,registers,1,0,source,good,faces)
      call require(good,'interface coarse and fine source assembly')
      call require(abs(source%fields(1)%values(1,8,7,8)-expected)<1.e-14_dp, &
        'interface source full-block kappa normalization and clipped fine dm')
      call require(source%fields(1)%values(1,7,8,8)==.5_dp,'interface source preserves weighted coarse delta')
      call require(.not.allocated(source%fields(2)%values),'interface source only coarse-level halo published')
      if(mode==1) then
        call require(source%fields(1)%values(1,11,8,8)==.5_dp,'interface periodic image takes accepted valid owner')
      else
        call require(source%fields(1)%values(1,11,8,8)==0,'interface physical exterior source absent')
      endif
      call assemble_state_redist_interface_rereflux_3d(tree,geometry,registers,1,0,increment,good,faces)
      call require(good,'interface source feeds EB re-reflux')
      mass=sum(geometry(1)%volume_fraction(7:10,7:10,7:10)*increment%fields(1)%values(1,:,:,:))
      call require(abs(mass-expected*geometry(1)%volume_fraction(8,7,8)-.5_dp)<1.e-13_dp, &
        'EB re-reflux conserves volume-weighted interface ledger')
      call require(increment%fields(1)%values(1,2,2,2)>0,'EB re-reflux retains fine-covered coarse recipient')
      call require(all(increment%fields(2)%values==0),'fine injection remains owned by tree flux synchronization')
      kept_increment=increment
      kept=source
      registers(2)%fine_delta(1,7,6,7)=ieee_value(0._dp,ieee_quiet_nan)
      call assemble_state_redist_interface_source_3d(tree,geometry,registers,1,0,source,good,faces)
      call require(.not.good.and.all(source%fields(1)%values==kept%fields(1)%values), &
        'interface invalid late fine ledger preserves output')
      call assemble_state_redist_interface_rereflux_3d(tree,geometry,registers,1,0,increment,good,faces)
      call require(.not.good.and.all(increment%fields(1)%values==kept_increment%fields(1)%values).and. &
        all(increment%fields(2)%values==kept_increment%fields(2)%values),'EB re-reflux joint output rollback')
      registers(2)%fine_delta=0
      registers(2)%fine_delta(1,7,6,7)=huge(1._dp)*.75_dp
      registers(1)%coarse_delta(1,8,7,8)=huge(1._dp)*.95_dp
      call assemble_state_redist_interface_source_3d(tree,geometry,registers,1,0,source,good,faces)
      call require(.not.good.and.all(source%fields(1)%values==kept%fields(1)%values), &
        'interface coarse-plus-fine overflow preserves output')
    enddo
  end subroutine check_interface_source
  subroutine check_register_accumulation()
    type(state_redist_register_3d) :: total,part,kept
    integer :: flags(2,2,2),mode
    logical :: ghost(2,2,2),cells(2,2,2),good
    flags=0;flags(1,1,1)=1;ghost=.false.;ghost(2,2,2)=.true.;cells=.not.ghost
    call initialize_state_redist_register_3d(2,[2,2,2],total,good,flags,ghost,.5_dp,cells)
    call require(good,'interval register constructor')
    part=total;part%coarse_delta=2;part%fine_delta=-3
    call accumulate_state_redist_register_3d(total,part,2,[2,2,2],good)
    call require(good.and.all(total%coarse_delta==2).and.all(total%fine_delta==-3), &
      'interval register already weighted seed')
    part%coarse_delta=-.5_dp;part%fine_delta=1
    call accumulate_state_redist_register_3d(total,part,2,[2,2,2],good)
    call require(good.and.all(total%coarse_delta==1.5_dp).and.all(total%fine_delta==-2).and.total%scale==.5_dp, &
      'interval signed sum without reapplying stage weight')
    ! Compiled fixed AMReX control gives coarse .125, .25, .375 after
    ! three identical calls without reset. Fine data remains endpoint-local.
    kept=total;part=total;part%coarse_delta=.125_dp;part%fine_delta=9
    total%coarse_delta=.125_dp;total%fine_delta=17
    call accumulate_state_redist_coarse_register_3d(total,part,2,[2,2,2],good)
    call require(good.and.all(total%coarse_delta==.25_dp).and.all(total%fine_delta==17), &
      'second coarse evaluation accumulates and preserves current fine')
    call accumulate_state_redist_coarse_register_3d(total,part,2,[2,2,2],good)
    call require(good.and.all(total%coarse_delta==.375_dp).and.all(total%fine_delta==17), &
      'third coarse evaluation follows compiled fixed lifecycle')
    total=kept
    do mode=1,6
      kept=total;part=total
      select case(mode)
      case(1)
        part%scale=1
      case(2)
        part%coarse_flags(2,1,1)=2
      case(3)
        part%fine_ghost(1,1,1)=.true.;part%valid_cells(1,1,1)=.false.
      case(4)
        part%valid_cells(2,1,1)=.false.
      case(5)
        part%coarse_delta(1,2,2,2)=ieee_value(0._dp,ieee_quiet_nan)
      case(6)
        total%coarse_delta=huge(1._dp)*.75_dp;kept=total
        part%coarse_delta=1
        part%coarse_delta(2,2,2,2)=huge(1._dp)*.75_dp
      end select
      call accumulate_state_redist_register_3d(total,part,2,[2,2,2],good)
      call require(.not.good,'interval register invalid or overflowing sum rejected')
      call require(all(total%coarse_delta==kept%coarse_delta).and.all(total%fine_delta==kept%fine_delta).and. &
        all(total%coarse_flags==kept%coarse_flags).and.all(total%fine_ghost.eqv.kept%fine_ghost).and. &
        all(total%valid_cells.eqv.kept%valid_cells).and.total%scale==kept%scale,'interval register joint rollback')
    enddo
  end subroutine check_register_accumulation
  subroutine check_update_scales()
    type(eb_geometry_3d) :: geometry
    type(state_redist_box_3d) :: box
    type(state_redist_register_3d) :: reg,direct,half,kept
    real(dp), allocatable :: u(:,:,:,:),rhs(:,:,:,:),rate(:,:,:,:),reference(:,:,:,:),provisional(:,:,:,:), &
      redistributed(:,:,:,:),scale(:,:,:),unit_scale(:,:,:),saved(:,:,:,:)
    integer, allocatable :: flags(:,:,:)
    logical, allocatable :: affected(:,:,:),valid_cells(:,:,:)
    integer :: mode,order,n(3),dims(3),i,j,k,bad
    real(dp) :: dt
    logical :: good
    character(len=16) :: plane_axis
    dt=.01_dp
    do mode=1,2
      n=[4,4,4];plane_axis='x'
      if(mode==2) then
        n=[6,1,1];plane_axis='y'
      endif
      dims=n+12
      call build_axis_plane_eb_geometry_3d(dims(1),dims(2),dims(3), &
        -6._dp/n(1),1._dp+6._dp/n(1),-6._dp/n(2),1._dp+6._dp/n(2), &
        -6._dp/n(3),1._dp+6._dp/n(3),trim(plane_axis),merge(.49_dp,.66_dp,mode==1),geometry,good)
      call require(good,'update scale geometry')
      box%valid_lower=7;box%valid_upper=n+6;box%domain_lower=7;box%domain_upper=n+6;box%periodic=.false.
      allocate(u(2,dims(1),dims(2),dims(3)),rhs(2,dims(1),dims(2),dims(3)),scale(dims(1),dims(2),dims(3)), &
        unit_scale(dims(1),dims(2),dims(3)),rate(2,dims(1),dims(2),dims(3)),reference(2,dims(1),dims(2),dims(3)), &
        provisional(2,dims(1),dims(2),dims(3)),redistributed(2,dims(1),dims(2),dims(3)), &
        affected(dims(1),dims(2),dims(3)),flags(dims(1),dims(2),dims(3)),valid_cells(dims(1),dims(2),dims(3)))
      flags=srd_coarse_fine_covered;flags(:8,:,:)=srd_coarse_boundary
      valid_cells=.false.;valid_cells(7:n(1)+6,7:n(2)+6,7:n(3)+6)=.true.
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        u(:,i,j,k)=[1._dp,2._dp]+.02_dp*(i+2*j+3*k)
        rhs(:,i,j,k)=[.2_dp,-.1_dp]+.01_dp*(2*i-j+k)
        scale(i,j,k)=.5_dp+.25_dp*modulo(i+j+k,5)
      enddo
      enddo
      enddo
      unit_scale=1
      do order=0,4,2
        call initialize_state_redist_register_3d(2,dims,reg,good,coarse_flags=flags,valid_cells=valid_cells)
        call require(good,'update scale register')
        direct=reg
        call redistribute_sdc_eb_rate_3d(geometry,u,rhs,dt,'state_redist',reference,good,max_order=order,registers=reg,box=box)
        call require(good,'update scale absent control')
        reg=direct
        call redistribute_sdc_eb_rate_3d(geometry,u,rhs,dt,'state_redist',rate,good,max_order=order, &
          registers=reg,box=box,update_scale=unit_scale)
        call require(good.and.all(rate==reference),'unit update scale exact default identity')
        provisional=u
        do k=1,dims(3)
        do j=1,dims(2)
        do i=1,dims(1)
          if(geometry%volume_fraction(i,j,k)<=0) then
            provisional(:,i,j,k)=0;cycle
          endif
          if(.not.box%contains_valid([i,j,k],3)) cycle
          if(box%contains_valid([i,j,k],1).and..not.box%contains_domain([i,j,k],1)) cycle
          provisional(:,i,j,k)=u(:,i,j,k)+dt*rhs(:,i,j,k)/scale(i,j,k)
        enddo
        enddo
        enddo
        call reactive_eb_weighted_state_redistribute_3d(geometry,provisional,redistributed,good,max_order=order, &
          registers=direct,box=box,affected=affected)
        call require(good,'nonuniform scale independent weighted control')
        do k=7,n(3)+6
        do j=7,n(2)+6
        do i=7,n(1)+6
          reference(:,i,j,k)=rhs(:,i,j,k)
          if(geometry%volume_fraction(i,j,k)<=0) then
            reference(:,i,j,k)=0
          else if(affected(i,j,k)) then
            reference(:,i,j,k)=scale(i,j,k)*(redistributed(:,i,j,k)-u(:,i,j,k))/dt
          endif
        enddo
        enddo
        enddo
        call initialize_state_redist_register_3d(2,dims,reg,good,coarse_flags=flags,valid_cells=valid_cells)
        call require(good,'nonuniform scale fresh register')
        call redistribute_sdc_eb_rate_3d(geometry,u,rhs,dt,'state_redist',rate,good,max_order=order, &
          registers=reg,box=box,update_scale=scale)
        call require(good,'nonuniform update scale evaluation')
        call require(maxval(abs(rate(:,7:n(1)+6,7:n(2)+6,7:n(3)+6)- &
          reference(:,7:n(1)+6,7:n(2)+6,7:n(3)+6)))<1.e-12_dp,'nonuniform scaled rate matches independent control')
        call require(maxval(abs(reg%coarse_delta-direct%coarse_delta))<1.e-12_dp,'update scale provisional ledger')
        call initialize_state_redist_register_3d(2,dims,half,good,coarse_flags=flags,valid_cells=valid_cells,scale=.5_dp)
        call require(good,'half register scale constructor')
        call redistribute_sdc_eb_rate_3d(geometry,u,rhs,dt,'state_redist',reference,good,max_order=order, &
          registers=half,box=box,update_scale=scale)
        call require(good.and.all(reference==rate),'half register scale leaves rate exact')
        call require(maxval(abs(half%coarse_delta-.5_dp*reg%coarse_delta))<1.e-12_dp,'half register scale ledger only')
      enddo
      saved=rate;kept=reg
      ! Unused/covered scale values are not used in scratch or rate arithmetic.
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        if(.not.box%contains_valid([i,j,k],3).or.geometry%volume_fraction(i,j,k)<=0) &
          scale(i,j,k)=ieee_value(0._dp,ieee_quiet_nan)
      enddo
      enddo
      enddo
      call initialize_state_redist_register_3d(2,dims,reg,good,coarse_flags=flags,valid_cells=valid_cells)
      call require(good,'poisoned update scale constructor')
      call redistribute_sdc_eb_rate_3d(geometry,u,rhs,dt,'state_redist',rate,good,max_order=4, &
        registers=reg,box=box,update_scale=scale)
      call require(good.and.all(rate==saved),'unused update scale poison ignored')
      kept=reg
      do bad=1,4
        scale(n(1)+6,7,7)=0
        if(bad==2) scale(n(1)+6,7,7)=-1
        if(bad==3) scale(n(1)+6,7,7)=ieee_value(0._dp,ieee_quiet_nan)
        if(bad==4) scale(n(1)+6,7,7)=tiny(1._dp)*epsilon(1._dp)
        call redistribute_sdc_eb_rate_3d(geometry,u,rhs,dt,'state_redist',rate,good,max_order=4, &
          registers=reg,box=box,update_scale=scale)
        call require(.not.good,'invalid or overflowing update scale rejected')
        call require(all(rate==saved).and.all(reg%coarse_delta==kept%coarse_delta),'update scale joint rollback')
      enddo
      deallocate(u,rhs,scale,unit_scale,rate,reference,provisional,redistributed,affected,flags,valid_cells,saved)
    enddo
  end subroutine check_update_scales

  subroutine check_level_stages()
    type(amr_patch_tree_3d) :: tree
    type(amr_patch_3d) :: patches(2)
    type(eb_geometry_3d) :: g(3)
    type(amr_tree_state_3d) :: old,rhs,rate,saved
    type(state_redist_register_3d), allocatable :: regs(:),kept(:)
    type(state_redist_register_3d) :: direct
    type(state_redist_box_3d) :: box
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: reference(:,:,:,:)
    real(dp) :: dx(3),lo(3),hi(3),dt
    integer :: p,n(3),dims(3),i,j,k,order
    logical :: good,owned(3)
    call initialize_amr_patch_3d(6,4,4,1,2,2,3,2,3,2,patches(1),good)
    call require(good,'level stage first child')
    call initialize_amr_patch_3d(6,4,4,3,4,2,3,2,3,2,patches(2),good)
    call require(good,'level stage peer child')
    call initialize_amr_patch_tree_3d([6,4,4],[1,1],patches,tree,good)
    call require(good,'level stage topology')
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,good)
    call require(good,'level stage layout')
    allocate(old%fields(3),rhs%fields(3),rate%fields(3))
    dt=.01_dp
    do p=1,3
      n=tree%nodes(p)%shape;dims=n+12
      dx=1._dp/real(domain(:,tree%nodes(p)%level+1),dp)
      lo=(real(origin(:,p),dp)-6._dp)*dx;hi=lo+real(dims,dp)*dx
      call build_axis_plane_eb_geometry_3d(dims(1),dims(2),dims(3),lo(1),hi(1),lo(2),hi(2),lo(3),hi(3), &
        'y',.34_dp,g(p),good)
      call require(good,'level stage geometry')
      allocate(old%fields(p)%values(2,dims(1),dims(2),dims(3)),rhs%fields(p)%values(2,dims(1),dims(2),dims(3)), &
        rate%fields(p)%values(2,n(1),n(2),n(3)))
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        old%fields(p)%values(:,i,j,k)=[1._dp,2._dp]+.03_dp*(lo(1)+(i-.5_dp)*dx(1))
        rhs%fields(p)%values(:,i,j,k)=[.2_dp,-.1_dp]+.02_dp*(lo(3)+(k-.5_dp)*dx(3))
      enddo
      enddo
      enddo
      rate%fields(p)%values=-7._dp
    enddo
    do order=0,4,2
      call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,0,dt,rate,regs,good,max_order=order)
      call require(good,'level stage root update')
      call require(all(rate%fields(2)%values==-7._dp).and.all(rate%fields(3)%values==-7._dp), &
        'level stage inactive rates preserved')
      call require(.not.allocated(regs(2)%fine_delta),'level stage inactive ledger preserved')
      call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,1,dt,rate,regs,good,max_order=order,register_scale=.25_dp)
      call require(good,'level stage fine update')
      do p=2,3
        n=tree%nodes(p)%shape;dims=n+12
        allocate(reference(2,dims(1),dims(2),dims(3)),source=0._dp)
        call initialize_state_redist_tree_register_3d(tree,p,2,6,direct,good,scale=.25_dp,box=box)
        call require(good,'level stage direct masks')
        call redistribute_sdc_eb_rate_3d(g(p),old%fields(p)%values,rhs%fields(p)%values,dt,'state_redist',reference, &
          good,max_order=order,registers=direct,box=box)
        call require(good,'level stage direct kernel')
        call require(maxval(abs(rate%fields(p)%values-reference(:,7:6+n(1),7:6+n(2),7:6+n(3))))<1.e-14_dp, &
          'level stage matches patch rate')
        call require(all(regs(p)%fine_delta==direct%fine_delta),'level stage matches patch ledger')
        deallocate(reference)
      enddo
      if(order<4) then
        rate%fields(2)%values=-7;rate%fields(3)%values=-7
        deallocate(regs)
      endif
    enddo
    saved=rate;kept=regs
    rhs%fields(3)%values(:,7,7,7)=ieee_value(1._dp,ieee_quiet_nan)
    call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,1,dt,rate,regs,good,max_order=4)
    call require(.not.good,'level stage later patch failure rejected')
    do p=1,3
      call require(all(rate%fields(p)%values==saved%fields(p)%values),'level stage all rates rollback')
      if(p==1) then
        call require(all(regs(p)%coarse_delta==kept(p)%coarse_delta),'level stage coarse ledger rollback')
      else
        call require(all(regs(p)%fine_delta==kept(p)%fine_delta),'level stage fine ledger rollback')
      endif
    enddo
    owned=[.false.,.true.,.false.]
    deallocate(old%fields(3)%values,rhs%fields(3)%values,rate%fields(3)%values)
    call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,1,dt,rate,regs,good,owned=owned)
    call require(good,'level stage owner-local absent peer')
    call require(.not.allocated(rate%fields(3)%values),'level stage absent peer stays absent')
    call require(all(regs(3)%fine_delta==kept(3)%fine_delta),'level stage unowned ledger retained')
    saved=rate;kept=regs
    g(2)%x_lower=g(2)%x_lower+.01_dp;g(2)%x_upper=g(2)%x_upper+.01_dp
    call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,1,dt,rate,regs,good,owned=owned)
    call require(.not.good,'level stage displaced halo geometry rejected')
    call require(all(rate%fields(2)%values==saved%fields(2)%values),'level stage geometry failure rollback')
    g(2)%x_lower=g(2)%x_lower-.01_dp;g(2)%x_upper=g(2)%x_upper-.01_dp
    owned=.false.
    call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,1,dt,rate,regs,good,owned=owned)
    call require(good,'level stage idle participant')
    call compute_state_redist_level_stage_3d(tree,g,old,rhs,2,1,dt,rate,regs,good,max_order=3,owned=owned)
    call require(.not.good,'level stage idle participant validates controls')
  end subroutine check_level_stages

  subroutine check_patch_boxes()
    type(eb_geometry_3d) :: g,hg,wide_g
    type(state_redist_box_3d) :: box,invalid_box,internal_box
    type(state_redist_register_3d) :: reg,native,kept
    type(amr_patch_tree_3d) :: tree
    type(amr_patch_3d) :: patch(1)
    real(dp), allocatable :: u(:,:,:,:),out(:,:,:,:),saved(:,:,:,:),poison(:,:,:,:),rhs(:,:,:,:),rate(:,:,:,:), &
      reference(:,:,:,:),halo(:,:,:,:),provisional(:,:,:,:)
    integer, allocatable :: flags(:,:,:),native_flags(:,:,:)
    logical, allocatable :: owned(:,:,:)
    integer :: mode,bc,order,n(3),dims(3),lo(3),hi(3),i,j,k,ci,cj,ck
    real(dp) :: error
    logical :: good
    character(len=16) :: faces(6),plane_axis
    do mode=1,2
      n=[4,4,4]
      if(mode==2) n=[6,1,1]
      dims=n+12;lo=7;hi=lo+n-1
      plane_axis='x'
      if(mode==2) plane_axis='y'
      call build_axis_plane_eb_geometry_3d(dims(1),dims(2),dims(3), &
        -6._dp/n(1),1._dp+6._dp/n(1),-6._dp/n(2),1._dp+6._dp/n(2), &
        -6._dp/n(3),1._dp+6._dp/n(3),trim(plane_axis),merge(.49_dp,.66_dp,mode==1),wide_g,good)
      call require(good,'patch-box wide geometry')
      call extract_eb_geometry_3d(wide_g,lo,hi,g,good)
      call require(good,'patch-box valid geometry')
      call extract_eb_geometry_3d(wide_g,lo-2,hi+2,hg,good)
      call require(good,'patch-box old boundary geometry')
      allocate(u(1,dims(1),dims(2),dims(3)),out(1,dims(1),dims(2),dims(3)), &
        flags(dims(1),dims(2),dims(3)),owned(dims(1),dims(2),dims(3)), &
        reference(1,n(1),n(2),n(3)),halo(1,n(1)+4,n(2)+4,n(3)+4), &
        native_flags(n(1),n(2),n(3)))
      owned=.false.;owned(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      flags=0;native_flags=0
      if(mode==1) then
        call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,patch(1),good)
        call require(good,'patch-box refined coarse bounds')
        call initialize_amr_patch_tree_3d(n,[1],patch,tree,good)
        call require(good,'patch-box topology')
      else
        flags(lo(1):lo(1)+1,lo(2),lo(3))=srd_coarse_boundary
        flags(lo(1)+2:hi(1),lo(2),lo(3))=srd_coarse_fine_covered
        native_flags=flags(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
      endif
      do bc=1,3
        faces='outflow'
        if(bc==2) faces(3)='inflow'
        if(bc==3) faces(3:4)='periodic'
        if(mode==2.and.bc==3) cycle
        do k=1,dims(3)
        do j=1,dims(2)
        do i=1,dims(1)
          ci=i-6;cj=j-6;ck=k-6
          if(bc==3) cj=1+modulo(cj-1,n(2))
          u(1,i,j,k)=3+.02_dp*ci*ci+.03_dp*cj+.04_dp*ck
          if(wide_g%volume_fraction(i,j,k)<=0) u(:,i,j,k)=0
        enddo
        enddo
        enddo
        halo=u(:,lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2)
        do order=0,4,2
          if(mode==1) then
            call initialize_state_redist_tree_register_3d(tree,1,1,6,reg,good,faces,box=box)
            call require(good.and.box%is_valid(dims),'patch register and box assembled jointly')
            native_flags=reg%coarse_flags(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
          else
            box%valid_lower=lo;box%valid_upper=hi;box%domain_lower=lo;box%domain_upper=hi;box%periodic=.false.
            call initialize_state_redist_register_3d(1,dims,reg,good,flags,valid_cells=owned)
            call require(good,'BFS patch-box coarse register')
          endif
          call initialize_state_redist_register_3d(1,n,native,good,coarse_flags=native_flags)
          call require(good,'native patch reference register')
          call reactive_eb_weighted_state_redistribute_3d(g,u(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
            reference,good,max_order=order,boundary_state=halo,boundary_geometry=hg,boundary_faces=faces,registers=native)
          call require(good,'native physical-domain StateRedist control')
          call reactive_eb_weighted_state_redistribute_3d(wide_g,u,out,good,max_order=order, &
            boundary_faces=faces,registers=reg,box=box)
          call require(good,'prepared patch-box StateRedist')
          error=maxval(abs(out(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))-reference))
          call require(error<2.e-12_dp,'prepared patch matches native physical/periodic boundary reconstruction')
          do k=1,dims(3)
          do j=1,dims(2)
          do i=1,dims(1)
            if(box%contains_valid([i,j,k],0)) cycle
            call require(all(out(:,i,j,k)==u(:,i,j,k)),'patch scatter preserves every halo cell')
          enddo
          enddo
          enddo
          call require(maxval(abs(reg%coarse_delta(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))- &
            native%coarse_delta))<2.e-12_dp,'prepared patch matches native coarse exchange ledger')
          saved=out;kept=reg;poison=u
          do k=1,dims(3)
          do j=1,dims(2)
          do i=1,dims(1)
            if(box%contains_valid([i,j,k],3)) cycle
            if(wide_g%volume_fraction(i,j,k)>0) poison(:,i,j,k)=.9_dp*huge(1._dp)
          enddo
          enddo
          enddo
          reg%coarse_delta=0
          call reactive_eb_weighted_state_redistribute_3d(wide_g,poison,out,good,max_order=order, &
            boundary_faces=faces,registers=reg,box=box)
          call require(good.and.maxval(abs(out(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))- &
            saved(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))))<2.e-12_dp,'unused outer state halo cannot alter valid update')
          call require(maxval(abs(reg%coarse_delta-kept%coarse_delta))<2.e-12_dp, &
            'unused outer state halo cannot alter register corrections')
          if(mode==2) then
            ! An interior patch needs cent_hat in grow-3 but no U_in from
            ! grow-4 to build it. Huge unused state values must never enter
            ! a discarded Qhat average or trigger an overflow trap.
            internal_box=box;internal_box%domain_lower(1)=1;internal_box%domain_upper(1)=dims(1)
            reg%coarse_delta=0
            call reactive_eb_weighted_state_redistribute_3d(wide_g,u,out,good,max_order=order, &
              boundary_faces=faces,registers=reg,box=internal_box)
            call require(good,'interior BFS patch-box control')
            saved=out;kept=reg;reg%coarse_delta=0
            call reactive_eb_weighted_state_redistribute_3d(wide_g,poison,out,good,max_order=order, &
              boundary_faces=faces,registers=reg,box=internal_box)
            call require(good.and.maxval(abs(out(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))- &
              saved(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))))<2.e-12_dp,'interior BFS ignores unused outer state halo')
            call require(maxval(abs(reg%coarse_delta-kept%coarse_delta))<2.e-12_dp, &
              'interior BFS ledger ignores unused outer state halo')
          endif
          rhs=0*u
          do k=1,dims(3)
          do j=1,dims(2)
          do i=1,dims(1)
            if(box%contains_domain([i,j,k],3).and.box%contains_valid([i,j,k],3)) rhs(:,i,j,k)=.01_dp*u(:,i,j,k)
          enddo
          enddo
          enddo
          if(mode==1) then
            rhs(:,hi(1),lo(2):hi(2),lo(3):hi(3))=1.e-30_dp
            if(bc<3) then
              rhs(:,hi(1),lo(2)-2,lo(3))=.123_dp
              rhs(:,hi(1),lo(2)-3,lo(3))=.456_dp
            endif
          endif
          provisional=u+.25_dp*rhs;reg%coarse_delta=0
          call reactive_eb_weighted_state_redistribute_3d(wide_g,provisional,out,good,max_order=order, &
            boundary_faces=faces,registers=reg,box=box)
          call require(good,'patch-box provisional rate control')
          kept=reg;reg%coarse_delta=0;rate=0*u
          call redistribute_sdc_eb_rate_3d(wide_g,u,rhs,.25_dp,'state_redist',rate,good,max_order=order, &
            boundary_faces=faces,registers=reg,box=box)
          call require(good.and.maxval(abs(u+.25_dp*rate-out))<2.e-12_dp,'SDC patch-box rate forwarding')
          if(mode==1) then
            call require(all(rate(:,hi(1),lo(2):hi(2),lo(3):hi(3))==1.e-30_dp), &
              'unaffected regular cells return tiny original forcing exactly')
            if(bc<3) then
              call require(all(rate(:,hi(1),lo(2)-2,lo(3))==.123_dp).and. &
                all(rate(:,hi(1),lo(2)-3,lo(3))==.456_dp),'second/third physical source ghosts retain supplied forcing')
            endif
          endif
          call require(maxval(abs(reg%coarse_delta-kept%coarse_delta))<2.e-12_dp,'SDC patch-box register forwarding')
          saved=rate;reg%coarse_delta=0
          do k=1,dims(3)
          do j=1,dims(2)
          do i=1,dims(1)
            if(box%contains_valid([i,j,k],3)) then
              if(.not.box%contains_valid([i,j,k],1).or.box%contains_domain([i,j,k],1)) cycle
            endif
            rhs(:,i,j,k)=.9_dp*huge(1._dp)
          enddo
          enddo
          enddo
          call redistribute_sdc_eb_rate_3d(wide_g,u,rhs,.25_dp,'state_redist',rate,good,max_order=order, &
            boundary_faces=faces,registers=reg,box=box)
          call require(good.and.maxval(abs(rate-saved))<2.e-12_dp,'SDC resets grow-1 physical and unused source ghosts')
          call require(maxval(abs(reg%coarse_delta-kept%coarse_delta))<2.e-12_dp, &
            'SDC reset ghost forcing cannot enter register')
          invalid_box=box;invalid_box%valid_lower(1)=5;kept=reg;saved=rate
          call redistribute_sdc_eb_rate_3d(wide_g,u,rhs,.25_dp,'state_redist',rate,good,max_order=order, &
            boundary_faces=faces,registers=reg,box=invalid_box)
          call require(.not.good.and.all(rate==saved).and.all(reg%coarse_delta==kept%coarse_delta), &
            'insufficient patch halo rolls back rate and register together')
        enddo
      enddo
      deallocate(u,out,flags,owned,reference,halo,native_flags)
      if(allocated(rhs)) deallocate(rhs)
      if(allocated(rate)) deallocate(rate)
    enddo
  end subroutine check_patch_boxes

  subroutine check_register_valid_box()
    type(state_redist_register_3d) :: reg,kept
    logical :: valid(5,1,1),ghost(5,1,1),good
    integer :: flags(5,1,1)
    ! Cell 1 is a same-level peer halo, 2/3 are this patch's valid cells,
    ! 4 is another peer halo, and 5 is a coarse-interpolated fine ghost.
    valid=.false.;valid(2:3,1,1)=.true.;ghost=.false.;ghost(5,1,1)=.true.
    flags(:,1,1)=[srd_coarse_boundary,srd_coarse_fine_covered,srd_coarse_boundary, &
      srd_coarse_fine_covered,srd_coarse_regular]
    call initialize_state_redist_register_3d(1,[5,1,1],reg,good,flags,ghost,valid_cells=valid)
    call require(good,'valid-box register construction')
    call accumulate_state_redist_transfer_3d(reg,[2,1,1],[1,1,1],[2._dp],.5_dp,.25_dp,good)
    call require(good.and.all(reg%coarse_delta==0),'covered valid donor does not debit peer-halo recipient')
    call accumulate_state_redist_transfer_3d(reg,[1,1,1],[2,1,1],[2._dp],.25_dp,.5_dp,good)
    call require(good.and.all(reg%coarse_delta==0),'peer-halo donor does not credit coarse register')
    call accumulate_state_redist_transfer_3d(reg,[4,1,1],[3,1,1],[2._dp],.5_dp,.25_dp,good)
    call require(good.and.reg%coarse_delta(1,3,1,1)==-2,'covered peer halo debits valid boundary recipient')
    call accumulate_state_redist_transfer_3d(reg,[3,1,1],[4,1,1],[2._dp],.25_dp,.5_dp,good)
    call require(good.and.reg%coarse_delta(1,3,1,1)==2,'valid boundary donor credits covered peer halo exchange')
    call accumulate_state_redist_transfer_3d(reg,[5,1,1],[4,1,1],[2._dp],.5_dp,.25_dp,good)
    call require(good.and.all(reg%fine_delta==0),'fine ghost to peer halo is outside valid recipient box')
    call accumulate_state_redist_transfer_3d(reg,[4,1,1],[5,1,1],[2._dp],.25_dp,.5_dp,good)
    call require(good.and.all(reg%fine_delta==0),'peer halo to fine ghost is outside valid donor box')
    call accumulate_state_redist_transfer_3d(reg,[5,1,1],[3,1,1],[2._dp],.5_dp,.25_dp,good)
    call require(good.and.reg%fine_delta(1,5,1,1)==-.5_dp,'fine ghost to valid uses recipient fluid fraction')
    call accumulate_state_redist_transfer_3d(reg,[3,1,1],[5,1,1],[2._dp],.25_dp,.5_dp,good)
    call require(good.and.reg%fine_delta(1,5,1,1)==.5_dp,'valid to fine ghost uses ghost fluid fraction')
    kept=reg;valid(5,1,1)=.true.
    call initialize_state_redist_register_3d(1,[5,1,1],reg,good,flags,ghost,valid_cells=valid)
    call require(.not.good.and.all(reg%fine_delta==kept%fine_delta).and.all(reg%valid_cells.eqv.kept%valid_cells), &
      'invalid overlapping ghost/valid masks roll back register construction')
    deallocate(reg%fine_ghost);allocate(reg%fine_ghost(2,1,1),source=.false.)
    call require(.not.reg%is_valid(1,[5,1,1]),'malformed ghost mask rejected before mask intersection')
  end subroutine check_register_valid_box

  subroutine check_tree_register_masks()
    type(amr_patch_tree_3d) :: tree
    type(amr_patch_3d) :: patches(3)
    type(state_redist_register_3d) :: root,left,right,leaf,kept
    logical :: refined(6,4,4),good
    integer, allocatable :: flags(:,:,:)
    character(len=16) :: faces(6)
    call initialize_amr_patch_3d(6,4,4,1,2,2,3,2,3,2,patches(1),good)
    call require(good,'register left peer bounds')
    call initialize_amr_patch_3d(6,4,4,3,4,2,3,2,3,2,patches(2),good)
    call require(good,'register right peer bounds')
    call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(3),good)
    call require(good,'register nested bounds')
    call initialize_amr_patch_tree_3d([6,4,4],[1,1,2],patches,tree,good)
    call require(good,'register nested topology')
    call initialize_state_redist_tree_register_3d(tree,1,1,3,root,good)
    call require(good.and.root%is_valid(1,[12,10,10]),'coarse patch halo register')
    call require(.not.allocated(root%fine_ghost).and.count(root%valid_cells)==6*4*4,'root has no fine-ghost ledger')
    refined=.false.;refined(1:4,2:3,2:3)=.true.
    call build_state_redist_coarse_flags_3d(refined,flags,good)
    call require(good.and.all(root%coarse_flags(4:9,4:7,4:7)==flags),'root refined union equals coarse flags')
    call require(all(root%coarse_delta==0),'fresh root register has zero correction')
    call initialize_state_redist_tree_register_3d(tree,2,1,3,left,good)
    call require(good.and.left%is_valid(1,[10,10,10]),'intermediate left register')
    call require(count(left%valid_cells)==64.and..not.left%valid_cells(8,4,4),'peer halo is outside left valid box')
    call require(.not.left%fine_ghost(8,4,4),'left high-x halo is a same-level peer cell')
    call require(left%fine_ghost(4,3,4),'interior level gap is a coarse-interpolated fine ghost')
    call require(.not.left%fine_ghost(3,4,4),'nonperiodic physical exterior is not a fine ghost')
    call require(left%coarse_flags(6,4,4)==srd_coarse_fine_covered,'nested coverage marks owned coarse cell')
    call initialize_state_redist_tree_register_3d(tree,3,1,3,right,good)
    call require(good.and.right%is_valid(1,[10,10,10]),'intermediate right register')
    call require(.not.right%fine_ghost(3,4,4),'right low-x halo belongs to left peer')
    call require(right%fine_ghost(8,4,4),'right high-x halo lies in an unrefined level gap')
    call require(right%coarse_flags(3,4,4)==srd_coarse_fine_covered,'fine coverage from peer children marks halo')
    call require(right%coarse_flags(4,4,4)==srd_coarse_boundary,'peer child makes owned cell a coarse boundary')
    call initialize_state_redist_tree_register_3d(tree,4,1,3,leaf,good)
    call require(good.and..not.allocated(leaf%coarse_flags).and.allocated(leaf%fine_ghost),'finest patch has fine ledger only')
    faces='outflow';faces(1:2)='periodic'
    call initialize_state_redist_tree_register_3d(tree,1,1,3,root,good,faces)
    call require(good.and.root%coarse_flags(9,5,5)==srd_coarse_boundary,'periodic seam sees opposite fine coverage')
    call require(root%coarse_flags(3,5,5)==srd_coarse_boundary,'negative halo indices wrap to high-domain coarse cell')
    call initialize_state_redist_tree_register_3d(tree,2,1,3,left,good,faces)
    call require(good.and.left%fine_ghost(3,4,4),'periodic low-x halo wraps to an unrefined level gap')
    kept=left;faces(2)='outflow'
    call initialize_state_redist_tree_register_3d(tree,2,1,3,left,good,faces)
    call require(.not.good.and.all(left%coarse_flags==kept%coarse_flags).and. &
      all(left%fine_ghost.eqv.kept%fine_ghost),'unpaired periodic face preserves prior register')
    call initialize_state_redist_tree_register_3d(tree,0,1,3,left,good)
    call require(.not.good.and.all(left%coarse_flags==kept%coarse_flags),'invalid node preserves prior register')
  end subroutine check_tree_register_masks

  subroutine check_register_restriction()
    type(state_redist_register_3d) :: reg
    real(dp) :: vf(4,3,2),delta(1,2,1,1),saved(1,2,1,1)
    logical :: ghost(4,3,2),active(2,1,1),refined(5,5,5),good
    integer, allocatable :: flags(:,:,:)
    integer :: i,j,k
    refined=.false.;refined(3,3,3)=.true.
    call build_state_redist_coarse_flags_3d(refined,flags,good)
    call require(good.and.flags(3,3,3)==srd_coarse_fine_covered.and. &
      flags(4,4,4)==srd_coarse_boundary.and.flags(1,1,1)==srd_coarse_regular,'coarse register coverage mask')
    ghost=.true.
    call initialize_state_redist_register_3d(1,[4,3,2],reg,good,fine_ghost=ghost)
    call require(good,'fine register restriction fixture')
    do k=1,2
    do j=1,3
    do i=1,4
      vf(i,j,k)=.05_dp*i+.02_dp*j+.01_dp*k
      reg%fine_delta(1,i,j,k)=3*vf(i,j,k)
    enddo
    enddo
    enddo
    delta=7;active=.false.;active(1,1,1)=.true.
    call restrict_state_redist_fine_delta_3d(reg,vf,[2,3,2],[-1,0,0],[-2,0,0],delta,active,good)
    call require(good,'anisotropic fine register projection with negative global origin')
    call assert_close(delta(1,1,1,1),10._dp,2.e-14_dp,'fine register uses full fluid-volume denominator')
    call assert_close(delta(1,2,1,1),7._dp,1.e-14_dp,'same-level peer mask excludes inactive register block')
    saved=delta
    call restrict_state_redist_fine_delta_3d(reg,vf,[2,3,2],[-1,0,0],[-1,0,0],delta,active,good)
    call require(.not.good.and.all(delta==saved),'incomplete fine volume block restriction rollback')
    vf=.25_dp;reg%fine_delta=1
    call restrict_state_redist_fine_delta_3d(reg,vf,[2,3,2],[-1,0,0],[-2,0,0],delta,active,good,threshold=.25_dp)
    call require(good.and.all(delta==saved),'fine register strict redistribution threshold')
    vf=.5_dp
    call restrict_state_redist_fine_delta_3d(reg,vf,[2,3,2],[-1,0,0],[-2,0,0],delta,active,good,threshold=.25_dp)
    call require(good,'fine register above redistribution threshold')
    call assert_close(delta(1,1,1,1),12._dp,2.e-14_dp,'fine register contribution above threshold')
  end subroutine check_register_restriction
  subroutine check_multilevel_registers()
    type(state_redist_register_3d) :: reg,saved,part_a,part_b
    type(eb_geometry_3d) :: g
    integer :: flags(4,1,1),mode,order,i
    logical :: ghost(4,1,1),good,valid_a(4,1,1),valid_b(4,1,1)
    real(dp) :: u(1,4,1,1),v(1,4,1,1),reference(1,4,1,1),coarse_mass,fine_mass,gas_change,rate(1,4,1,1),rhs(1,4,1,1)
    ! Direct bidirectional formula checks, including scale and volume ratio.
    flags=srd_coarse_regular;flags(1,1,1)=srd_coarse_boundary;flags(2,1,1)=srd_coarse_fine_covered
    ghost=.false.;ghost(2,1,1)=.true.
    call initialize_state_redist_register_3d(1,[4,1,1],reg,good,flags,ghost,2._dp)
    call require(good,'ML StateRedist register initialization')
    call accumulate_state_redist_transfer_3d(reg,[2,1,1],[1,1,1],[3._dp],.8_dp,.2_dp,good)
    call require(good,'covered to uncovered StateRedist contribution')
    call assert_close(reg%coarse_delta(1,1,1,1),-6._dp,1.e-14_dp,'coarse covered-to-boundary sign')
    call assert_close(reg%fine_delta(1,2,1,1),-1.2_dp,1.e-14_dp,'fine ghost-to-valid volume factor')
    call accumulate_state_redist_transfer_3d(reg,[1,1,1],[2,1,1],[1._dp],.2_dp,.8_dp,good)
    call require(good,'uncovered to covered StateRedist contribution')
    call assert_close(reg%coarse_delta(1,1,1,1),2._dp,1.e-14_dp,'coarse boundary-to-covered volume ratio')
    call assert_close(reg%fine_delta(1,2,1,1),.4_dp,1.e-14_dp,'fine valid-to-ghost volume factor')
    saved=reg
    call accumulate_state_redist_transfer_3d(reg,[1,1,1],[2,1,1],[1._dp],0._dp,.8_dp,good)
    call require(.not.good.and.all(reg%coarse_delta==saved%coarse_delta).and. &
      all(reg%fine_delta==saved%fine_delta),'invalid ML contribution rollback')
    reg%scale=1;reg%coarse_delta(1,1,1,1)=-.75_dp*huge(1._dp);saved=reg
    call accumulate_state_redist_transfer_3d(reg,[2,1,1],[1,1,1],[.5_dp*huge(1._dp)],.8_dp,.2_dp,good)
    call require(.not.good.and.all(reg%coarse_delta==saved%coarse_delta).and. &
      all(reg%fine_delta==saved%fine_delta),'overflowing ML ledger addition rollback')

    do mode=1,2
      if(mode==1) then
        call build_axis_plane_eb_geometry_3d(4,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.49_dp,g,good)
        flags=srd_coarse_regular;flags(2,1,1)=srd_coarse_boundary;flags(3,1,1)=srd_coarse_fine_covered
        ghost=.false.;ghost(3,1,1)=.true.
      else
        ! A thin connected region has no full-volume receiver, exercising BFS.
        call build_axis_plane_eb_geometry_3d(4,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'y',.84_dp,g,good)
        flags=srd_coarse_boundary;flags(3:4,1,1)=srd_coarse_fine_covered
        ghost=.false.;ghost(3:4,1,1)=.true.
      endif
      call require(good,'ML register physical geometry')
      do i=1,4
        u(1,i,1,1)=1+real(i*i,dp)
        if(g%volume_fraction(i,1,1)<=0) u(:,i,1,1)=0
      enddo
      do order=0,4,2
        call reactive_eb_weighted_state_redistribute_3d(g,u,reference,good,max_order=order)
        call require(good,'ML register gas reference')
        call initialize_state_redist_register_3d(1,[4,1,1],reg,good,flags,ghost,.5_dp)
        call require(good,'ML register neighborhood initialization')
        call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order,registers=reg)
        call require(good.and.all(v==reference),'ML register preserves reconstructed gas update')
        gas_change=sum(g%volume_fraction*(v(1,:,:,:)-u(1,:,:,:)),mask=flags==srd_coarse_boundary)
        coarse_mass=sum(g%volume_fraction*reg%coarse_delta(1,:,:,:))
        fine_mass=sum(reg%fine_delta(1,:,:,:))
        call assert_close(coarse_mass,-.5_dp*gas_change,2.e-12_dp,'ML coarse ledger cancels interface gas transfer')
        call assert_close(fine_mass,coarse_mass,2.e-12_dp,'ML fine mass and coarse density ledgers agree')
        call require(abs(coarse_mass)>1.e-10_dp,'ML ledger is active')
        saved=reg
        ! The two valid boxes partition level-valid cells. Peer halos remain
        ! non-ghost but are outside each patch's ownership. Their assembled
        ! register sum must equal a single box, for primary and BFS operators.
        valid_a=.false.;valid_a(1,1,1)=.not.ghost(1,1,1)
        valid_b=.not.ghost;valid_b(1,1,1)=.false.
        call initialize_state_redist_register_3d(1,[4,1,1],part_a,good,flags,ghost,.5_dp,valid_a)
        call require(good,'first partition valid-box register')
        call initialize_state_redist_register_3d(1,[4,1,1],part_b,good,flags,ghost,.5_dp,valid_b)
        call require(good,'second partition valid-box register')
        call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order,registers=part_a)
        call require(good.and.all(v==reference),'first valid-box register preserves reconstructed state')
        call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order,registers=part_b)
        call require(good.and.all(v==reference),'second valid-box register preserves reconstructed state')
        call require(maxval(abs(part_a%coarse_delta+part_b%coarse_delta-saved%coarse_delta))< &
          64*epsilon(1._dp)*max(1._dp,maxval(abs(saved%coarse_delta))), &
          'partitioned valid boxes preserve full coarse interface ledger')
        call require(maxval(abs(part_a%fine_delta+part_b%fine_delta-saved%fine_delta))< &
          64*epsilon(1._dp)*max(1._dp,maxval(abs(saved%fine_delta))), &
          'partitioned valid boxes preserve full fine ghost ledger')
        call initialize_state_redist_register_3d(1,[4,1,1],reg,good,flags,ghost,.5_dp)
        rhs=0
        call redistribute_sdc_eb_rate_3d(g,u,rhs,.25_dp,'state_redist',rate,good,max_order=order,registers=reg)
        call require(good.and.maxval(abs(u+.25_dp*rate-reference))<1.e-13_dp,'SDC ML reconstructed rate')
        call require(maxval(abs(reg%coarse_delta-saved%coarse_delta))<1.e-13_dp.and. &
          maxval(abs(reg%fine_delta-saved%fine_delta))<1.e-13_dp,'SDC ML register forwarding')
        saved=reg;reg%scale=-1
        call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order,registers=reg)
        call require(.not.good.and.all(reg%coarse_delta==saved%coarse_delta).and. &
          all(reg%fine_delta==saved%fine_delta),'invalid ML register preserves accumulated ledger')
      enddo
    enddo
  end subroutine check_multilevel_registers
  subroutine check_periodic_state_halos()
    type(eb_geometry_3d) :: g,hg
    real(dp) :: phi(0:8,0:8,0:8),hp(0:12,0:12,0:12),wave(0:7)
    real(dp) :: u(1,8,8,8),v(1,8,8,8),halo(1,12,12,12),reference(1,8,8,8), &
      saved(1,8,8,8),before,after
    character(len=16) :: faces(6)
    integer :: i,j,k,ci,cj,ck,order,shift,mode
    logical :: good
    wave=[0._dp,sqrt(.5_dp),1._dp,sqrt(.5_dp),0._dp,-sqrt(.5_dp),-1._dp,-sqrt(.5_dp)]
    faces='outflow';faces(1:2)='periodic'
    do mode=1,2
    do order=0,4,2
    do shift=0,2,2
      do i=0,8
        phi(i,:,:)=wave(modulo(i-shift,8))+.02_dp
      enddo
      do i=0,12
        hp(i,:,:)=wave(modulo(i-2-shift,8))+.02_dp
      enddo
      if(mode==1) then
      call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],g,good)
      call require(good,'periodic seam physical geometry')
      call build_level_set_eb_geometry_3d(hp,[-.25_dp,-.25_dp,-.25_dp],[1.25_dp,1.25_dp,1.25_dp],hg,good)
      call require(good,'periodic seam halo geometry')
      else
        ! No cell reaches the target alone: BFS neighborhoods must cross the seam.
        call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'y',.98_dp,g,good)
        call require(good,'periodic connected thin geometry')
        call build_axis_plane_eb_geometry_3d(12,12,12,-.25_dp,1.25_dp,-.25_dp,1.25_dp, &
          -.25_dp,1.25_dp,'y',.98_dp,hg,good)
        call require(good,'periodic connected thin halo')
      endif
      do k=1,8
      do j=1,8
      do i=1,8
        u(1,i,j,k)=2+.1_dp*wave(modulo(i-1-shift,8))+.02_dp*j*j+.03_dp*k
        if(g%volume_fraction(i,j,k)<=0) u(:,i,j,k)=0
      enddo
      enddo
      enddo
      do k=1,12
      do j=1,12
      do i=1,12
        ci=1+modulo(i-3,8);cj=max(1,min(8,j-2));ck=max(1,min(8,k-2))
        halo(:,i,j,k)=u(:,ci,cj,ck)
      enddo
      enddo
      enddo
      before=sum(g%volume_fraction*u(1,:,:,:))
      call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order, &
        boundary_state=halo,boundary_geometry=hg,boundary_faces=faces)
      call require(good,'periodic seam receiver and Qhat reconstruction')
      after=sum(g%volume_fraction*v(1,:,:,:))
      call assert_close(after,before,2.e-12_dp,'periodic seam conservative scatter')
      if(shift==0) then
        reference=v
      else
        do i=1,8
          call require(maxval(abs(v(:,i,:,:)-reference(:,1+modulo(i-3,8),:,:)))<2.e-12_dp, &
            'periodic StateRedist translation across seam')
        enddo
      endif
      saved=v
      ! Periodic interior donors must use updated Qhat, never stale exterior U.
      halo(:,1:2,3:10,3:10)=99;halo(:,11:12,3:10,3:10)=-99
      call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order, &
        boundary_state=halo,boundary_geometry=hg,boundary_faces=faces)
      call require(good.and.maxval(abs(v-saved))<1.e-14_dp,'periodic ghost state independent Qhat')
    enddo
    enddo
    enddo
    faces(2)='outflow'
    call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=2, &
      boundary_state=halo,boundary_geometry=hg,boundary_faces=faces)
    call require(.not.good.and.all(v==0),'unpaired periodic StateRedist rejection')
  end subroutine check_periodic_state_halos
  subroutine check_boundary_state_halos()
    type(eb_geometry_3d) :: g,hg
    real(dp) :: u(1,8,8,8),v(1,8,8,8),truncated(1,8,8,8),halo(1,12,12,12),rhs(1,8,8,8),rate(1,8,8,8), &
      x,y,z,before,after
    character(len=16) :: faces(6)
    integer :: i,j,k,ci,cj,ck,order
    logical :: good
    call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,g,good)
    call require(good,'physical boundary state geometry')
    call build_axis_plane_eb_geometry_3d(12,12,12,-.25_dp,1.25_dp,-.25_dp,1.25_dp, &
      -.25_dp,1.25_dp,'x',.37_dp,hg,good)
    call require(good,'physical boundary state halo geometry')
    u=0
    do k=1,8
    do j=1,8
    do i=1,8
      if(g%volume_fraction(i,j,k)<=0) cycle
      x=(i-.5_dp+g%cell_centroid_x(i,j,k))/8
      y=(j-.5_dp+g%cell_centroid_y(i,j,k))/8
      z=(k-.5_dp+g%cell_centroid_z(i,j,k))/8
      u(1,i,j,k)=2+.1_dp*x+.2_dp*y*y+.03_dp*z
    enddo
    enddo
    enddo
    do k=1,12
    do j=1,12
    do i=1,12
      ci=max(1,min(8,i-2));cj=max(1,min(8,j-2));ck=max(1,min(8,k-2))
      halo(:,i,j,k)=u(:,ci,cj,ck)
      if(j<=2.and.g%volume_fraction(ci,cj,ck)>0) then
        x=(ci-.5_dp+g%cell_centroid_x(ci,cj,ck))/8
        z=(ck-.5_dp+g%cell_centroid_z(ci,cj,ck))/8
        halo(1,i,j,k)=2+.1_dp*x+.03_dp*z
      endif
    enddo
    enddo
    enddo
    faces='outflow';faces(3)='inflow'
    before=sum(g%volume_fraction*u(1,:,:,:))
    do order=2,4,2
      call reactive_eb_weighted_state_redistribute_3d(g,u,truncated,good,max_order=order)
      call require(good,'truncated boundary comparison')
      call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order, &
        boundary_state=halo,boundary_geometry=hg,boundary_faces=faces)
      call require(good,'physical boundary Qhat halo reconstruction')
      call require(maxval(abs(v-truncated))>1.e-10_dp,'boundary halo affects cut-cell redistribution')
      after=sum(g%volume_fraction*v(1,:,:,:))
      call assert_close(after,before,2.e-12_dp,'physical boundary halo conservative scatter')
      rhs=.1_dp*u
      call reactive_eb_weighted_state_redistribute_3d(g,u+.25_dp*rhs,v,good,max_order=order, &
        boundary_state=halo,boundary_geometry=hg,boundary_faces=faces)
      call require(good,'boundary halo provisional reference')
      call redistribute_sdc_eb_rate_3d(g,u,rhs,.25_dp,'state_redist',rate,good,max_order=order, &
        boundary_state=halo,boundary_geometry=hg,boundary_faces=faces)
      call require(good.and.maxval(abs(u+.25_dp*rate-v))<1.e-14_dp,'SDC boundary Qhat rate forwarding')

    enddo
    call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=2,boundary_state=halo)
    call require(.not.good.and.all(v==0),'incomplete boundary halo rejection')
  end subroutine check_boundary_state_halos

  subroutine check_boundary_graph_stencil()
    real(dp) :: q(1,5),centers(3,5),centroids(3,5),widths(3,5),fractions(5),slopes(3,1)
    integer :: bounds(2,3)
    logical :: faces(6),good
    centers=0;centers(1,:)=[0._dp,-2._dp,-1._dp,1._dp,2._dp]
    centroids=centers;centroids(1,2:3)=-.5_dp
    q(1,:)=[2._dp,1._dp,1._dp,4._dp,7._dp];widths=1;fractions=1
    faces=.false.;faces(1)=.true.;bounds(1,:)=-10;bounds(2,:)=10;bounds(1,1)=0
    call state_redistribution_graph_slopes_3d(q,centroids,centers,widths,fractions,4,slopes,good, &
      extdir_faces=faces,extdir_bounds=bounds)
    call require(good,'face-sampled graph boundary reconstruction')
    call require(maxval(abs(slopes(:,1)-[1.9_dp,0._dp,0._dp]))<2.e-13_dp, &
      'graph dispatch uses extdir regular slope')
    call state_redistribution_graph_slopes_3d(q,centroids,centers,widths,fractions,4,slopes,good, &
      extdir_faces=faces)
    call require(.not.good.and.all(slopes==0),'graph boundary metadata pair rejection')
    bounds(1,1)=1
    call state_redistribution_graph_slopes_3d(q,centroids,centers,widths,fractions,4,slopes,good, &
      extdir_faces=faces,extdir_bounds=bounds)
    call require(.not.good.and.all(slopes==0),'graph boundary domain rejection')
  end subroutine check_boundary_graph_stencil

  subroutine check_grown_graph_stencil()
    real(dp) :: q(1,8),centers(3,8),widths(3,8),fractions(8),slopes(3,1)
    logical :: connected(8),good
    integer :: order,a
    widths=1;fractions=.5_dp;centers=0;q=0
    centers(:,2)=[0._dp,-1._dp,0._dp]
    centers(:,3)=[0._dp,1._dp,0._dp]
    centers(:,4)=[2._dp,0._dp,0._dp]
    q(1,:4)=[0._dp,-1._dp,1._dp,4._dp]
    connected=.false.;connected(:3)=.true.
    do order=2,4,2
      call state_redistribution_graph_slopes_3d(q(:,:4),centers(:,:4),centers(:,:4), &
        widths(:,:4),fractions(:4),order,slopes,good,connected(:4))
      call require(good,'grown nonconnected graph stencil')
      call require(maxval(abs(slopes(:,1)-[2._dp,1._dp,0._dp]))<2.e-13_dp, &
        'grown stencil uses noncovered disconnected donor')
    enddo
    centers=0;q=0
    do a=1,3
      centers(a,2*a)=-1;centers(a,2*a+1)=1
      q(1,2*a)=-real(a,dp);q(1,2*a+1)=real(a,dp)
    enddo
    centers(:,8)=1;q(1,8)=100
    connected=.true.;connected(8)=.false.
    call state_redistribution_graph_slopes_3d(q,centers,centers,widths,fractions,2,slopes,good,connected)
    call require(good,'normal connected graph stencil')
    call require(maxval(abs(slopes(:,1)-[1._dp,2._dp,3._dp]))<2.e-13_dp, &
      'normal stencil excludes disconnected donor')
    call state_redistribution_graph_slopes_3d(q,centers,centers,widths,fractions,2,slopes,good,connected(:1))
    call require(.not.good.and.all(slopes==0),'invalid graph connection mask rejection')
  end subroutine check_grown_graph_stencil

  subroutine check_high_order_state()
    type(eb_geometry_3d) :: g
    real(dp) :: u(3,8,8,8),v(3,8,8,8),zero(3,8,8,8),phi(0:8,0:8,0:8),position(3),before,after
    integer :: i,j,k,c,order,mode
    logical :: good
    do mode=1,2
      if(mode==1) then
        call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,g,good)
      else
        do k=0,8
        do j=0,8
        do i=0,8
          phi(i,j,k)=real(i,dp)/8+.23_dp*real(j,dp)/8+.17_dp*real(k,dp)/8-.43_dp
        enddo
        enddo
        enddo
        call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],g,good)
      endif
      call require(good,'higher-order centroid geometry')
      u=0
      do k=1,8
      do j=1,8
      do i=1,8
        if(g%volume_fraction(i,j,k)<=0) cycle
        position=real([i,j,k],dp)+[g%cell_centroid_x(i,j,k),g%cell_centroid_y(i,j,k),g%cell_centroid_z(i,j,k)]
        u(:,i,j,k)=[1+.1_dp*position(1)+.03_dp*position(2)-.02_dp*position(3), &
          2-.04_dp*position(1)+.01_dp*position(2)+.05_dp*position(3), &
          .5_dp+.02_dp*position(1)-.01_dp*position(2)+.03_dp*position(3)]
      enddo
      enddo
      enddo
      call reactive_eb_weighted_state_redistribute_3d(g,u,zero,good)
      call require(good,'zeroth-order centroid control')
      call require(maxval(abs(zero-u))>1.e-6_dp,'centroid affine fixture requires slope reconstruction')
      do order=2,4,2
        call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=order)
        call require(good,'higher-order weighted centroid reconstruction')
        call require(maxval(abs(v-u))<2.e-12_dp,'higher-order affine field is preserved at cell centroids')
        do c=1,3
          before=sum(u(c,:,:,:)*g%volume_fraction);after=sum(v(c,:,:,:)*g%volume_fraction)
          call assert_close(after,before,2.e-12_dp*max(1._dp,abs(before)),'higher-order component inventory')
        enddo
      enddo
      call reactive_eb_weighted_state_redistribute_3d(g,u,v,good,max_order=1)
      call require(.not.good.and.all(v==0),'unsupported reconstruction order rejects output')
    enddo
  end subroutine check_high_order_state

  subroutine check_sdc_rates()
    type(eb_geometry_3d) :: geometry
    real(dp) :: state(2,8,1,1),rhs(2,8,1,1),rate(2,8,1,1),expected(2,8,1,1),provisional(2,8,1,1)
    real(dp) :: clean(2,8,1,1),saved(2,8,1,1),before,after,dt,nan_value
    character(len=12), parameter :: methods(3)=[character(len=12)::'none','flux_redist','state_redist']
    integer :: i,c,mode,order
    logical :: ok
    call build_axis_plane_eb_geometry_3d(8,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.37_dp,geometry,ok)
    call require(ok,'SDC redistribution geometry')
    dt=.125_dp;nan_value=ieee_value(0._dp,ieee_quiet_nan)
    do i=1,8
      state(:,i,1,1)=[real(i,dp),real(9-i,dp)]
      rhs(:,i,1,1)=[(-1._dp)**i,2*(-1._dp)**(i+1)]
      if(geometry%cell_type(i,1,1)==eb_covered_cell_3d) then
        state(:,i,1,1)=nan_value;rhs(:,i,1,1)=nan_value
      end if
    end do
    do mode=1,3
      clean=state;provisional=rhs
      do i=1,8
        if(geometry%cell_type(i,1,1)/=eb_covered_cell_3d) cycle
        clean(:,i,1,1)=0;provisional(:,i,1,1)=0
      end do
      call redistribute_sdc_eb_rate_3d(geometry,state,rhs,dt,trim(methods(mode)),rate,ok)
      call require(ok,'SDC redistribution source stage')
      do c=1,2
        before=sum(provisional(c,:,:,:)*geometry%volume_fraction)
        after=sum(rate(c,:,:,:)*geometry%volume_fraction)
        call assert_close(after,before,2e-12_dp*max(1._dp,abs(before)),'SDC signed rate conservation')
      end do
      if(mode==3) then
        provisional=clean+dt*provisional
        call reactive_eb_weighted_state_redistribute_3d(geometry,provisional,expected,ok)
        call require(ok.and.maxval(abs(clean+dt*rate-expected))<1e-13_dp,'SDC rate reconstructs redistributed state')
      end if
      do i=1,8
        if(geometry%cell_type(i,1,1)==eb_covered_cell_3d) &
          call require(all(rate(:,i,1,1)==0),'covered SDC rates zero')
      end do
    end do
    ! A one-dimensional stencil exercises the rank-deficient centroid solve.
    clean=0;rhs=0
    do i=1,8
      if(geometry%volume_fraction(i,1,1)<=0) cycle
      clean(:,i,1,1)=[1._dp+.1_dp*(i+geometry%cell_centroid_x(i,1,1)),2._dp]
      rhs(:,i,1,1)=[.03_dp*(i+geometry%cell_centroid_x(i,1,1)),.1_dp]
    enddo
    do order=2,4,2
      call redistribute_sdc_eb_rate_3d(geometry,clean,rhs,dt,'state_redist',rate,ok,max_order=order)
      call require(ok.and.maxval(abs(rate-rhs))<2.e-12_dp,'high-order SDC affine rate on 1D stencil')
      provisional=clean+dt*rhs
      call reactive_eb_weighted_state_redistribute_3d(geometry,provisional,expected,ok,max_order=order)
      call require(ok.and.maxval(abs(clean+dt*rate-expected))<2.e-12_dp,'high-order SDC order propagation')
    enddo
    saved=rate
    call redistribute_sdc_eb_rate_3d(geometry,clean,rhs,dt,'state_redist',rate,ok,max_order=1)
    call require(.not.ok.and.all(rate==saved),'invalid SDC reconstruction order rollback')
    saved=rate;rhs(1,8,1,1)=nan_value
    call redistribute_sdc_eb_rate_3d(geometry,state,rhs,dt,'state_redist',rate,ok)
    call require(.not.ok.and.all(rate==saved),'invalid active SDC rate rollback')
  end subroutine check_sdc_rates


  subroutine check_weighted_boundary_pair()
    type(eb_geometry_3d) :: geometry
    real(dp) :: state(2,2,1,1),result(2,2,1,1),expected(2,2,1,1)
    logical :: good
    integer :: n
    call build_axis_plane_eb_geometry_3d(2,1,1,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'y',.66_dp,geometry,good)
    call require(good,'weighted physical-boundary pair geometry')
    state(:,1,1,1)=[1._dp,-3._dp];state(:,2,1,1)=[3._dp,7._dp]
    ! Native StateRedistUtils: kappa=17/50, alpha_neighbor=8/17,
    ! each cell belongs to two neighborhoods, alpha_self=13/17.
    ! Applying the weighted average and scatter gives these rational weights.
    expected(:,1,1,1)=(185._dp*state(:,1,1,1)+104._dp*state(:,2,1,1))/289._dp
    expected(:,2,1,1)=(104._dp*state(:,1,1,1)+185._dp*state(:,2,1,1))/289._dp
    call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,good)
    call require(good,'weighted physical-boundary pair redistribution')
    call require(maxval(abs(result-expected))<2.e-13_dp,'native rational overlap weights at physical boundary')
    do n=1,2
      call assert_close(sum(result(n,:,:,:)*geometry%volume_fraction), &
        sum(state(n,:,:,:)*geometry%volume_fraction),2.e-13_dp,'weighted boundary pair inventory')
    enddo
  end subroutine check_weighted_boundary_pair
  subroutine check_connected_clusters()
    type(eb_geometry_3d) :: geometry
    type(eb_geometry_2d) :: section
    real(dp), allocatable :: state(:,:,:,:),result(:,:,:,:)
    real(dp) :: phi(0:4,0:4),before,after
    integer :: axis,i,j,k,component,dims(3),order
    character(len=1), parameter :: axes(3)=["x","y","z"]
    logical :: ok

    do axis=1,3
      dims=[8,8,8]; dims(axis)=1
      call build_axis_plane_eb_geometry_3d(dims(1),dims(2),dims(3), &
        0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,axes(axis),.96_dp,geometry,ok)
      call require(ok,'thin boundary sheet geometry')
      call require(maxval(geometry%volume_fraction)<.5_dp,'no large receiver exists')
      allocate(state(2,dims(1),dims(2),dims(3)),result(2,dims(1),dims(2),dims(3)))
      do k=1,dims(3)
        do j=1,dims(2)
          do i=1,dims(1)
            state(:,i,j,k)=[real(i+3*j+7*k,dp),real(2*i-j+k,dp)]
          end do
        end do
      end do
      call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok)
      call require(ok,'connected small-cell fallback')
      do component=1,2
        before=sum(state(component,:,:,:)*geometry%volume_fraction)
        after=sum(result(component,:,:,:)*geometry%volume_fraction)
        call assert_close(after,before,2.e-12_dp*max(1._dp,abs(before)),'cluster conservation')
        call require(minval(result(component,:,:,:))>=minval(state(component,:,:,:))-1.e-12_dp, &
          'cluster convex lower bound')
        call require(maxval(result(component,:,:,:))<=maxval(state(component,:,:,:))+1.e-12_dp, &
          'cluster convex upper bound')
      end do
      call require(maxval(abs(result-state))>1.e-4_dp,'cluster redistribution active')
      do order=2,4,2
        call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok,max_order=order)
        call require(ok.and.maxval(abs(result-state))<2.e-12_dp,'high-order affine connected sheet')
      enddo
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        state(:,i,j,k)=[sin(real(i+3*j+7*k,dp)),cos(real(2*i-j+k,dp))]
      enddo
      enddo
      enddo
      do order=2,4,2
        call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok,max_order=order)
        call require(ok,'high-order nonlinear connected sheet')
        do component=1,2
          before=sum(state(component,:,:,:)*geometry%volume_fraction)
          after=sum(result(component,:,:,:)*geometry%volume_fraction)
          call assert_close(after,before,2.e-12_dp*max(1._dp,abs(before)),'high-order nonlinear cluster inventory')
          call require(minval(result(component,:,:,:))>=minval(state(component,:,:,:))-2.e-12_dp.and. &
            maxval(result(component,:,:,:))<=maxval(state(component,:,:,:))+2.e-12_dp,'high-order cluster bounds')
        enddo
      enddo
      state=7._dp
      do order=0,4,2
        call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok,max_order=order)
        call require(ok.and.maxval(abs(result-state))<1.e-12_dp,'cluster constants preserved')
      enddo
      deallocate(state,result)
    end do

    ! Two thin fluid strips separated by covered cells: no cross-solid mixing.
    do j=0,4
      do i=0,4
        phi(i,j)=abs(real(i,dp)/4._dp-.5_dp)-.45_dp
      end do
    end do
    call build_eb_geometry_2d(phi,0._dp,1._dp,0._dp,1._dp,section,ok)
    call require(ok,'separated strips section')
    call extrude_eb_geometry_3d(section,3,1,0._dp,1._dp,geometry,ok)
    call require(ok,'separated strips geometry')
    allocate(state(1,4,4,1),result(1,4,4,1))
    state=99._dp; state(:,1,:,:)=2._dp; state(:,4,:,:)=8._dp
    call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok)
    call require(ok,'separate connected components')
    call require(maxval(abs(result(:,1,:,:)-2._dp))<1.e-12_dp.and. &
      maxval(abs(result(:,4,:,:)-8._dp))<1.e-12_dp,'solid prevents neighborhood mixing')
    call require(all(result(:,2:3,:,:)==0._dp),'cluster covered cells zero')
    ! Each connected component has volume .8, below a requested target of 1.
    call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok,1._dp)
    call require(.not.ok.and.all(result==0._dp),'insufficient connected volume fails transactionally')
  end subroutine check_connected_clusters

  subroutine check_curved_state()
    type(eb_geometry_2d) :: section
    type(eb_geometry_3d) :: geometry
    real(dp) :: phi(0:10,0:8),before,after
    real(dp), allocatable :: state(:,:,:,:),result(:,:,:,:)
    integer :: i,j,k,axis,component
    logical :: ok
    do j=0,8
      do i=0,10
        phi(i,j)=(real(i,dp)/10-.5_dp)**2+(real(j,dp)/8-.5_dp)**2-.29_dp**2
      end do
    end do
    call build_eb_geometry_2d(phi,0._dp,1._dp,0._dp,1._dp,section,ok)
    call require(ok,'curved section for StateRedist')
    do axis=1,3
      call extrude_eb_geometry_3d(section,axis,3,0._dp,1._dp,geometry,ok)
      call require(ok,'curved StateRedist geometry')
      allocate(state(2,geometry%nx,geometry%ny,geometry%nz),result(2,geometry%nx,geometry%ny,geometry%nz))
      do k=1,geometry%nz
        do j=1,geometry%ny
          do i=1,geometry%nx
            state(:,i,j,k)=[real(i+3*j+7*k,dp),real(2*i-j+k,dp)]
          end do
        end do
      end do
      call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok)
      call require(ok,'nonplanar StateRedist')
      do component=1,2
        before=sum(state(component,:,:,:)*geometry%volume_fraction)
        after=sum(result(component,:,:,:)*geometry%volume_fraction)
        call assert_close(after,before,2.e-12_dp*max(1._dp,abs(before)),'curved weighted conservation')
        call require(minval(result(component,:,:,:),mask=geometry%volume_fraction>0)>= &
          minval(state(component,:,:,:),mask=geometry%volume_fraction>0)-1.e-12_dp,'curved convex lower bound')
        call require(maxval(result(component,:,:,:))<=maxval(state(component,:,:,:))+1.e-12_dp, &
          'curved convex upper bound')
      end do
      call require(maxval(abs(result-state),mask=spread(geometry%volume_fraction>0,1,2))>1.e-4_dp, &
        'curved cells actually redistribute')
      state=7._dp
      call reactive_eb_weighted_state_redistribute_3d(geometry,state,result,ok)
      call require(ok,'curved uniform StateRedist')
      call require(maxval(abs(result-7._dp),mask=spread(geometry%volume_fraction>0,1,2))<1.e-12_dp, &
        'curved uniform state preserved')
      call require(maxval(abs(result),mask=spread(geometry%volume_fraction==0,1,2))==0._dp, &
        'curved covered output zero')
      deallocate(state,result)
    end do
  end subroutine check_curved_state


  subroutine check_regular_and_covered()
    integer, parameter :: nx = 3
    integer, parameter :: ny = 2
    integer, parameter :: nz = 2
    type(eb_geometry_3d) :: geometry
    real(dp) :: rhs(2, nx, ny, nz), redistributed(2, nx, ny, nz)
    logical :: ok
    integer :: i, j, k

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rhs(:, i, j, k) = &
            [real(i + 3 * j + 7 * k, dp), real(2 * i - j + k, dp)]
        end do
      end do
    end do
    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", -0.1_dp, geometry, ok)
    call require(ok, "fully regular geometry")
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok .and. all(redistributed == rhs), &
      "fully regular redistribution identity")

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 1.1_dp, geometry, ok)
    call require(ok, "fully covered geometry")
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok .and. maxval(abs(redistributed)) == 0.0_dp, &
      "fully covered redistribution is zero")

    geometry%volume_fraction(1, 1, 1) = 2.0_dp
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
      "invalid geometry rejection")
  end subroutine check_regular_and_covered

  subroutine check_axis(axis, nx, ny, nz, cut, neighbor)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: nx, ny, nz, cut(3), neighbor(3)

    integer, parameter :: ncomp = 5
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: rhs(:, :, :, :), redistributed(:, :, :, :)
    real(dp) :: reference(ncomp), neighborhood_value
    real(dp) :: expected_cut, expected_neighbor
    real(dp) :: original_integral, redistributed_integral
    real(dp) :: kappa, tolerance
    logical :: ok
    integer :: component, i, j, k

    allocate(rhs(ncomp, nx, ny, nz), redistributed(ncomp, nx, ny, nz))
    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.36875_dp, geometry, ok)
    call require(ok .and. geometry%is_valid(), &
      trim(axis)//"-axis small-cell geometry")
    call require(geometry%cell_type(cut(1), cut(2), cut(3)) == &
      eb_cut_cell_3d, trim(axis)//"-axis cut-cell location")
    kappa = geometry%volume_fraction(cut(1), cut(2), cut(3))
    call assert_close(kappa, 0.05_dp, 2.0e-13_dp, &
      trim(axis)//"-axis small volume fraction")

    reference = [1000.0_dp, -40.0_dp, 7.0_dp, 13.0_dp, -2800.0_dp]
    rhs = 0.0_dp
    rhs(:, cut(1), cut(2), cut(3)) = reference
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok, trim(axis)//"-axis flux redistribution")
    do component = 1, ncomp
      neighborhood_value = &
        kappa * reference(component) / (1.0_dp + kappa)
      expected_cut = kappa * reference(component) + &
        (1.0_dp - kappa) * neighborhood_value
      expected_neighbor = kappa * (1.0_dp - kappa) * &
        (reference(component) - neighborhood_value)
      tolerance = 3.0e-13_dp * max(1.0_dp, abs(reference(component)))
      call assert_close( &
        redistributed(component, cut(1), cut(2), cut(3)), &
        expected_cut, tolerance, trim(axis)//"-axis stabilized cut rhs")
      call assert_close( &
        redistributed(component, neighbor(1), neighbor(2), neighbor(3)), &
        expected_neighbor, tolerance, &
        trim(axis)//"-axis receiving-cell rhs")
      original_integral = sum( &
        geometry%volume_fraction * rhs(component, :, :, :))
      redistributed_integral = sum( &
        geometry%volume_fraction * redistributed(component, :, :, :))
      call assert_close(redistributed_integral, original_integral, &
        tolerance, trim(axis)//"-axis volume-weighted conservation")
    end do
    call require(abs(redistributed(1, cut(1), cut(2), cut(3))) < &
      0.11_dp * abs(rhs(1, cut(1), cut(2), cut(3))), &
      trim(axis)//"-axis small-cell stiffness reduction")
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (all([i, j, k] == cut) .or. &
              all([i, j, k] == neighbor)) cycle
          call require(maxval(abs(redistributed(:, i, j, k))) == 0.0_dp, &
            trim(axis)//"-axis compact support")
        end do
      end do
    end do

    rhs = 0.0_dp
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (geometry%cell_type(i, j, k) /= eb_covered_cell_3d) &
            rhs(:, i, j, k) = reference
        end do
      end do
    end do
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok, trim(axis)//"-axis uniform redistribution")
    call require(maxval(abs(redistributed - rhs)) <= &
      3.0e-13_dp * maxval(abs(rhs)), &
      trim(axis)//"-axis uniform active rhs preservation")

    if (axis == "x") then
      rhs = 0.0_dp
      rhs(1, cut(1), cut(2), cut(3)) = &
        ieee_value(0.0_dp, ieee_quiet_nan)
      call reactive_eb_flux_redistribute_3d( &
        geometry, rhs, redistributed, ok)
      call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
        "nonfinite rhs rejection")
    end if
  end subroutine check_axis

  subroutine check_state_axis(axis, nx, ny, nz, cut, neighbor)
    character(len=*), intent(in) :: axis
    integer, intent(in) :: nx, ny, nz, cut(3), neighbor(3)

    integer, parameter :: ncomp = 3
    type(eb_geometry_3d) :: geometry
    real(dp), allocatable :: provisional(:, :, :, :)
    real(dp), allocatable :: redistributed(:, :, :, :)
    real(dp) :: reference(ncomp), kappa, alpha
    real(dp) :: expected_cut, expected_neighbor
    real(dp) :: original_integral, redistributed_integral, tolerance
    logical :: ok
    integer :: component, i, j, k

    allocate(provisional(ncomp, nx, ny, nz))
    allocate(redistributed(ncomp, nx, ny, nz))
    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      axis, 0.36875_dp, geometry, ok)
    call require(ok, trim(axis)//"-axis weighted-state geometry")
    kappa = geometry%volume_fraction(cut(1), cut(2), cut(3))
    alpha = (0.5_dp - kappa) / &
      geometry%volume_fraction(neighbor(1), neighbor(2), neighbor(3))
    expected_cut = (-kappa + 0.5_dp * alpha) / &
      (kappa + 0.5_dp * alpha)
    expected_neighbor = 1.0_dp - 0.5_dp * alpha + &
      0.5_dp * alpha * expected_cut

    reference = [0.31_dp, -1.7_dp, 2.4_dp]
    do component = 1, ncomp
      provisional(component, :, :, :) = reference(component)
    end do
    provisional(:, cut(1), cut(2), cut(3)) = -reference
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) &
            provisional(:, i, j, k) = 9.0_dp * reference
        end do
      end do
    end do
    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok)
    call require(ok, trim(axis)//"-axis weighted StateRedist")
    do component = 1, ncomp
      tolerance = 5.0e-13_dp * max(1.0_dp, abs(reference(component)))
      call assert_close( &
        redistributed(component, cut(1), cut(2), cut(3)), &
        expected_cut * reference(component), tolerance, &
        trim(axis)//"-axis weighted cut state")
      call assert_close( &
        redistributed(component, neighbor(1), neighbor(2), neighbor(3)), &
        expected_neighbor * reference(component), tolerance, &
        trim(axis)//"-axis weighted receiving state")
      original_integral = sum( &
        geometry%volume_fraction * provisional(component, :, :, :))
      redistributed_integral = sum( &
        geometry%volume_fraction * redistributed(component, :, :, :))
      call assert_close(redistributed_integral, original_integral, &
        tolerance, trim(axis)//"-axis weighted-state conservation")
    end do
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) then
            call require(maxval(abs(redistributed(:, i, j, k))) == &
              0.0_dp, trim(axis)//"-axis covered state is zero")
          else if (.not. all([i, j, k] == cut) .and. &
                   .not. all([i, j, k] == neighbor)) then
            call require(all(redistributed(:, i, j, k) == reference), &
              trim(axis)//"-axis weighted compact support")
          end if
        end do
      end do
    end do

    do component = 1, ncomp
      provisional(component, :, :, :) = reference(component)
    end do
    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok)
    call require(ok, trim(axis)//"-axis uniform weighted StateRedist")
    call require(maxval(abs(redistributed - provisional), &
      mask=spread(geometry%cell_type /= eb_covered_cell_3d, 1, ncomp)) <= &
      8.0e-13_dp, trim(axis)//"-axis uniform state preservation")
    call require(maxval(abs(redistributed), &
      mask=spread(geometry%cell_type == eb_covered_cell_3d, 1, ncomp)) == &
      0.0_dp, trim(axis)//"-axis uniform covered state is zero")

    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok, kappa)
    call require(ok .and. maxval(abs(redistributed - provisional), &
      mask=spread(geometry%cell_type /= eb_covered_cell_3d, 1, ncomp)) == &
      0.0_dp, trim(axis)//"-axis custom target identity")

    call reactive_eb_weighted_state_redistribute_3d( &
      geometry, provisional, redistributed, ok, 0.0_dp)
    call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
      trim(axis)//"-axis invalid weighted target transaction")
    if (axis == "x") then
      provisional(1, cut(1), cut(2), cut(3)) = &
        ieee_value(0.0_dp, ieee_quiet_nan)
      call reactive_eb_weighted_state_redistribute_3d( &
        geometry, provisional, redistributed, ok)
      call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
        "nonfinite weighted state transaction")
      call build_axis_plane_eb_geometry_3d( &
        nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
        axis, 0.99375_dp, geometry, ok)
      call require(ok, "boundary small-cell geometry")
      provisional = 1.0_dp
      call reactive_eb_weighted_state_redistribute_3d( &
        geometry, provisional, redistributed, ok)
      call require(.not. ok .and. maxval(abs(redistributed)) == 0.0_dp, &
        "missing regular receiver transaction")
    end if
  end subroutine check_state_axis

  subroutine check_transverse_neighbors()
    integer, parameter :: nx = 8
    integer, parameter :: ny = 3
    integer, parameter :: nz = 3
    type(eb_geometry_3d) :: geometry
    real(dp) :: rhs(1, nx, ny, nz), redistributed(1, nx, ny, nz)
    real(dp) :: original_integral, redistributed_integral
    logical :: ok

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 0.36875_dp, geometry, ok)
    call require(ok, "transverse-neighbor geometry")
    rhs = 0.0_dp
    rhs(1, 3, 2, 2) = 1.0_dp
    call reactive_eb_flux_redistribute_3d( &
      geometry, rhs, redistributed, ok)
    call require(ok, "transverse-neighbor redistribution")
    call require(redistributed(1, 4, 2, 2) > 0.0_dp, &
      "positive-x neighbor receives excess")
    call require(redistributed(1, 3, 1, 2) > 0.0_dp .and. &
      redistributed(1, 3, 3, 2) > 0.0_dp, &
      "both y neighbors receive excess")
    call require(redistributed(1, 3, 2, 1) > 0.0_dp .and. &
      redistributed(1, 3, 2, 3) > 0.0_dp, &
      "both z neighbors receive excess")
    original_integral = sum(geometry%volume_fraction * rhs(1, :, :, :))
    redistributed_integral = sum( &
      geometry%volume_fraction * redistributed(1, :, :, :))
    call assert_close(redistributed_integral, original_integral, &
      5.0e-14_dp, "overlapping-neighborhood conservation")
  end subroutine check_transverse_neighbors

  subroutine check_reactive_advance()
    integer, parameter :: nx = 8
    integer, parameter :: ny = 1
    integer, parameter :: nz = 1
    integer, parameter :: cut_i = 3
    integer, parameter :: neighbor_i = 4
    type(nasa7_species), allocatable :: species(:)
    type(eb_geometry_3d) :: geometry
    type(state_redist_register_3d) :: registers
    integer :: flags(nx,ny,nz)
    logical :: ghost(nx,ny,nz)
    real(dp), allocatable :: primitive(:), state_cell(:)
    real(dp), allocatable :: mass_fractions(:)
    real(dp), allocatable :: state(:, :, :, :), rhs(:, :, :, :)
    real(dp), allocatable :: new_state(:, :, :, :)
    real(dp), allocatable :: temperature(:, :, :)
    real(dp), allocatable :: new_temperature(:, :, :)
    real(dp) :: mole_fractions(7), cell_temperature, sound_speed
    real(dp) :: kappa, alpha, expected_cut_scale, expected_neighbor_scale
    real(dp) :: original_integral, advanced_integral, tolerance
    logical :: ok
    integer :: component, nvar, species_index,order,i

    call load_h2o2_elementary_thermo(species, ok)
    call require(ok, "reactive thermodynamic database load")
    nvar = reactive_nvar(size(species))
    allocate(primitive(reactive_nprim(size(species))))
    allocate(state_cell(nvar), mass_fractions(size(species)))
    allocate(state(nvar, nx, ny, nz), rhs(nvar, nx, ny, nz))
    allocate(new_state(nvar, nx, ny, nz))
    allocate(temperature(nx, ny, nz), new_temperature(nx, ny, nz))

    mole_fractions = [0.29570_dp, 1.0e-5_dp, 1.0e-5_dp, 0.14784_dp, &
      1.0e-5_dp, 0.0_dp, 0.55643_dp]
    call mass_fractions_from_mole_fractions( &
      species, mole_fractions, mass_fractions, ok)
    call require(ok, "reactive composition conversion")
    primitive(1:5) = [0.31_dp, 0.0_dp, 0.0_dp, 0.0_dp, 135000.0_dp]
    do species_index = 1, size(species)
      primitive(reactive_mass_fraction_component(species_index)) = &
        mass_fractions(species_index)
    end do
    call reactive_primitive_to_conserved( &
      species, primitive, state_cell, cell_temperature, sound_speed, ok)
    call require(ok, "reactive reference-state construction")
    do component = 1, nvar
      state(component, :, :, :) = state_cell(component)
    end do
    temperature = cell_temperature

    call build_axis_plane_eb_geometry_3d( &
      nx, ny, nz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      "x", 0.36875_dp, geometry, ok)
    call require(ok, "reactive small-cell geometry")
    kappa = geometry%volume_fraction(cut_i, 1, 1)
    rhs = 0.0_dp
    rhs(:, cut_i, 1, 1) = -2.0_dp * state_cell
    call require(state(irho, cut_i, 1, 1) + &
      rhs(irho, cut_i, 1, 1) < 0.0_dp, &
      "unredistributed update is nonphysical")
    call advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(ok, "positive redistributed reactive advance")
    expected_cut_scale = 1.0_dp - 4.0_dp * kappa / (1.0_dp + kappa)
    expected_neighbor_scale = 1.0_dp - &
      2.0_dp * kappa * (1.0_dp - kappa) / (1.0_dp + kappa)
    call assert_close(new_state(irho, cut_i, 1, 1), &
      expected_cut_scale * state_cell(irho), 3.0e-13_dp, &
      "positive redistributed cut density")
    call assert_close(new_state(irho, neighbor_i, 1, 1), &
      expected_neighbor_scale * state_cell(irho), 3.0e-13_dp, &
      "positive redistributed neighbor density")
    call assert_close(new_temperature(cut_i, 1, 1), cell_temperature, &
      3.0e-9_dp, "redistributed cut temperature recovery")
    do component = 1, nvar
      original_integral = sum(geometry%volume_fraction * &
        rhs(component, :, :, :))
      advanced_integral = sum(geometry%volume_fraction * &
        (new_state(component, :, :, :) - state(component, :, :, :)))
      tolerance = 5.0e-13_dp * max(1.0_dp, abs(original_integral))
      call assert_close(advanced_integral, original_integral, tolerance, &
        "reactive advanced-state conservation")
    end do

    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(ok, "positive weighted-state reactive advance")
    flags=srd_coarse_regular;flags(cut_i,1,1)=srd_coarse_boundary
    flags(neighbor_i,1,1)=srd_coarse_fine_covered;ghost=.false.;ghost(neighbor_i,1,1)=.true.
    call initialize_state_redist_register_3d(nvar,[nx,ny,nz],registers,ok,flags,ghost)
    call require(ok,'reactive ML register initialization')
    temperature(cut_i,1,1)=0
    call advance_reactive_eb_state_redistributed_3d(species,state,temperature,geometry,rhs,1._dp, &
      new_state,new_temperature,ok,registers=registers)
    call require(.not.ok.and.all(registers%coarse_delta==0).and.all(registers%fine_delta==0), &
      'reactive EOS failure rolls back ML registers')
    temperature(cut_i,1,1)=cell_temperature
    call advance_reactive_eb_state_redistributed_3d(species,state,temperature,geometry,rhs,1._dp, &
      new_state,new_temperature,ok,registers=registers)
    call require(ok.and.any(registers%coarse_delta/=0),'reactive EOS success publishes ML registers')
    alpha = (0.5_dp - kappa) / &
      geometry%volume_fraction(neighbor_i, 1, 1)
    expected_cut_scale = (-kappa + 0.5_dp * alpha) / &
      (kappa + 0.5_dp * alpha)
    expected_neighbor_scale = 1.0_dp - 0.5_dp * alpha + &
      0.5_dp * alpha * expected_cut_scale
    call assert_close(new_state(irho, cut_i, 1, 1), &
      expected_cut_scale * state_cell(irho), 3.0e-13_dp, &
      "positive weighted-state cut density")
    call assert_close(new_state(irho, neighbor_i, 1, 1), &
      expected_neighbor_scale * state_cell(irho), 3.0e-13_dp, &
      "positive weighted-state neighbor density")
    call assert_close(new_temperature(cut_i, 1, 1), cell_temperature, &
      3.0e-9_dp, "weighted-state cut temperature recovery")
    do component = 1, nvar
      original_integral = sum(geometry%volume_fraction * &
        (state(component, :, :, :) + rhs(component, :, :, :)))
      advanced_integral = sum(geometry%volume_fraction * &
        new_state(component, :, :, :))
      tolerance = 5.0e-13_dp * max(1.0_dp, abs(original_integral))
      call assert_close(advanced_integral, original_integral, tolerance, &
        "weighted-state reactive advance conservation")
    end do

    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok, kappa)
    call require(.not. ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "custom weighted target reactive rollback")

    call advance_reactive_eb_state_redistributed_3d( &
      species,state,temperature,geometry,rhs,1._dp,new_state,new_temperature,ok,max_order=1)
    call require(.not.ok.and.all(new_state==state).and.all(new_temperature==temperature), &
      'invalid reactive reconstruction order rollback')
    ! Affine density scaling retains composition, velocity and temperature.
    rhs=0
    do i=1,nx
      state(:,i,1,1)=state_cell*(1._dp+.05_dp*(i+geometry%cell_centroid_x(i,1,1)))
    enddo
    do order=2,4,2
      call advance_reactive_eb_state_redistributed_3d( &
        species,state,temperature,geometry,rhs,1._dp,new_state,new_temperature,ok,max_order=order)
      call require(ok,'high-order reactive affine advance')
      do i=1,nx
        if(geometry%volume_fraction(i,1,1)<=0) cycle
        call require(maxval(abs(new_state(:,i,1,1)-state(:,i,1,1))/max(1._dp,abs(state(:,i,1,1)))) &
          <2.e-12_dp,'high-order reactive affine state preserved')
        call assert_close(new_temperature(i,1,1),cell_temperature,3.e-9_dp,'high-order EOS temperature')
      enddo
    enddo
    do component=1,nvar
      state(component,:,:,:)=state_cell(component)
    enddo
    rhs(:, cut_i, 1, 1) = -50.0_dp * state_cell
    call advance_reactive_eb_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(.not. ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "nonphysical redistributed advance rollback")
    call advance_reactive_eb_state_redistributed_3d( &
      species, state, temperature, geometry, rhs, 1.0_dp, &
      new_state, new_temperature, ok)
    call require(.not. ok .and. all(new_state == state) .and. &
      all(new_temperature == temperature), &
      "nonphysical weighted-state advance rollback")
  end subroutine check_reactive_advance

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2(1x,es24.16))') trim(message), actual, expected
      error stop message
    end if
  end subroutine assert_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine require

end program test_eb_reactive_redistribution_3d
