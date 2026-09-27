program test_rotating_eb_tree_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved, &
    reactive_conserved_to_primitive
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d,fill_reactive_eb_tree_halo_3d,compute_mol_eb_tree_hydro_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: boxes(2)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,extended(3),g
  type(eb_reactive_tree_3d) :: tree,saved
  type(amr_tree_state_3d) :: rates,old_rates
  type(eb_integrated_flux_3d), allocatable :: flux(:),old_flux(:)
  type(rotating_frame_3d) :: frame
  real(dp), allocatable :: q(:),back(:),root(:,:,:,:),temp(:,:,:),ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),inlet(:,:)
  real(dp) :: t,c,potential,position(3),lo(3),hi(3),spacing(3),expected_t,error,tin(6),mach(6),pin(6)
  integer :: nv,i,j,k,node,order,solver,checks,width
  logical :: ok
  character(len=16) :: faces(6)
  checks=0
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  nv=reactive_nvar(size(species))
  allocate(q(reactive_nprim(size(species))),back(reactive_nprim(size(species))), &
    root(nv,4,4,4),temp(4,4,4),ghost(nv,-2:7,-2:7,-2:7),gt(-2:7,-2:7,-2:7),vf(-2:7,-2:7,-2:7))
  frame%angular_speed=100;frame%axis_point=[10._dp,10._dp,0._dp]
  q=0;q(1)=1;q(2:4)=[1._dp,2._dp,-1._dp];q(5)=3e5_dp
  do i=1,size(species)
    if(species(i)%name=='N2') q(5+i)=1
  end do
  do k=1,4
  do j=1,4
  do i=1,4
    position=(real([i,j,k],dp)-.5_dp)/4
    call rotating_potential_3d(frame,position,potential,ok);call require(ok,'root potential')
    call reactive_primitive_to_conserved(species,q,root(:,i,j,k),t,c,ok,potential)
    call require(ok,'root EOS');temp(i,j,k)=t
  end do
  end do
  end do
  call require(maxval(root(5,:,:,:))<0,'negative rotating energy')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,2,3,2,boxes(1),ok);call require(ok,'left box')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,2,3,2,boxes(2),ok);call require(ok,'right box')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1],boxes,topology,ok);call require(ok,'topology')
  call build_axis_plane_eb_geometry_3d(8,8,8,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    'z',-10._dp,finest,ok);call require(ok,'geometry')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,temp,tree,ok,frame)
  call require(ok,'rotating conservative tree initialization')
  saved=tree
  call average_down_reactive_eb_tree_3d(species,tree,ok,frame);call require(ok,'rotating restriction')
  call require(all(tree%nodes(1)%state==root),'exact PCM/restriction conserved round trip')
  call require(maxval(abs(tree%nodes(1)%temperature/temp-1))<1e-12_dp,'root temperature round trip')
  do node=1,3
    g=tree%nodes(node)%geometry
    spacing=[g%dx,g%dy,g%dz];lo=[g%x_lower,g%y_lower,g%z_lower]-3*spacing
    hi=[g%x_upper,g%y_upper,g%z_upper]+3*spacing
    call build_axis_plane_eb_geometry_3d(10,10,10,lo(1),hi(1),lo(2),hi(2),lo(3),hi(3), &
      'z',-10._dp,extended(node),ok);call require(ok,'extended geometry')
  end do
  allocate(inlet(nv,6))
  faces='outflow';error=0
  do node=1,3
    call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,ok,ghost_volume_fraction=vf,halo_width=3, &
      target_volume_fraction=extended(node)%volume_fraction,boundary_faces=faces,frame=frame)
    call require(ok,'rotating tree halo')
    g=tree%nodes(node)%geometry
    do k=-2,7
    do j=-2,7
    do i=-2,7
      position=[g%x_lower,g%y_lower,g%z_lower]+(real([i,j,k],dp)-.5_dp)*[g%dx,g%dy,g%dz]
      call rotating_potential_3d(frame,position,potential,ok);call require(ok,'halo potential')
      call reactive_conserved_to_primitive(species,ghost(:,i,j,k),gt(i,j,k),back,expected_t,c,ok,potential)
      call require(ok,'halo EOS')
      error=max(error,abs(expected_t/gt(i,j,k)-1))
    end do
    end do
    end do
    if(node==2) then
      call require(all(ghost(:,5,2,2)==tree%nodes(3)%state(:,1,2,2)),'same-level peer conserved halo')
      call require(all(ghost(:,2,5,2)==tree%nodes(1)%state(:,1,3,2)),'ancestor conserved halo')
    end if
  end do
  call require(error<1e-12_dp,'receiver-position EOS consistency')
  call reactive_primitive_to_conserved(species,q,inlet(:,1),t,c,ok)
  call require(ok,'ordinary prescribed inlet');inlet=spread(inlet(:,1),2,6)
  tin=t
  do width=1,3,2
    faces='outflow';faces(1)='inflow'
    call fill_reactive_eb_tree_halo_3d(species,tree,1,ghost(:,1-width:4+width,1-width:4+width,1-width:4+width), &
      gt(1-width:4+width,1-width:4+width,1-width:4+width),ok, &
      ghost_volume_fraction=vf(1-width:4+width,1-width:4+width,1-width:4+width),halo_width=width, &
      target_volume_fraction=extended(1)%volume_fraction(4-width:7+width,4-width:7+width,4-width:7+width), &
      boundary_faces=faces,inflow_state=inlet,inflow_temperature=tin,frame=frame)
    call require(ok,'prescribed rotating inlet')
    position=(real([0,2,2],dp)-.5_dp)/4
    call rotating_potential_3d(frame,position,potential,ok);call require(ok,'inlet potential')
    call require(ghost(5,0,2,2)==inlet(5,1)-inlet(1,1)*potential,'inlet energy conversion')
    call require(abs(gt(0,2,2)/t-1)<1e-12_dp,'inlet temperature')
    mach=-1;mach(1)=.2_dp;pin=0;pin(1)=3e5_dp;tin=900
    call fill_reactive_eb_tree_halo_3d(species,tree,1,ghost(:,1-width:4+width,1-width:4+width,1-width:4+width), &
      gt(1-width:4+width,1-width:4+width,1-width:4+width),ok, &
      ghost_volume_fraction=vf(1-width:4+width,1-width:4+width,1-width:4+width),halo_width=width, &
      target_volume_fraction=extended(1)%volume_fraction(4-width:7+width,4-width:7+width,4-width:7+width), &
      inflow_mach=mach,inflow_pressure=pin,inflow_temperature=tin,frame=frame)
    call require(ok,'Mach rotating inlet')
    call reactive_conserved_to_primitive(species,ghost(:,0,2,2),gt(0,2,2),back,expected_t,c,ok,potential)
    call require(ok,'Mach inlet EOS')
    call require(abs(expected_t/900-1)<1e-12_dp.and.abs(back(5)/3e5_dp-1)<1e-12_dp.and. &
      abs(back(2)/c-.2_dp)<1e-12_dp,'Mach inlet imposed state')
    tin=t
  end do
  faces='outflow'
  do solver=0,1
  do order=1,2
    call compute_mol_eb_tree_hydro_3d(species,tree,extended,1e-9_dp,order,rates,flux,ok, &
      boundary_faces=faces,lax_friedrichs=solver==1,frame=frame)
    call require(ok,'rotating tree MOL spatial operator')
  end do
  end do
  old_rates=rates;old_flux=flux;frame%axis=0
  call compute_mol_eb_tree_hydro_3d(species,tree,extended,1e-9_dp,2,rates,flux,ok,frame=frame)
  call require(.not.ok,'invalid frame rejected')
  do node=1,3
    call require(all(rates%fields(node)%values==old_rates%fields(node)%values),'rate rollback')
    call require(all(flux(node)%x==old_flux(node)%x).and.all(flux(node)%y==old_flux(node)%y).and. &
      all(flux(node)%z==old_flux(node)%z),'flux rollback')
  end do
  call average_down_reactive_eb_tree_3d(species,tree,ok,frame)
  call require(.not.ok,'invalid restriction frame')
  do node=1,3
    call require(all(tree%nodes(node)%state==saved%nodes(node)%state),'restriction rollback')
  end do
  print *, 'PASS: rotating tree transfer and MOL; checks, halo EOS error:',checks,error
contains
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if(.not.condition) then
      print *, 'FAIL: ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine
end program
