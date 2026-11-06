module eb_reactive_tree_checkpoint_3d_mod
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d,amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,integrate_reactive_eb_tree_3d, &
    average_down_reactive_eb_tree_3d
  use reactive_eb_3d_driver_mod, only: reactive_eb_integrals_3d,reactive_eb_element_integrals_3d, &
    reactive_eb_composition_integrals_3d
  use reactive_eb_3d_checkpoint_mod, only: write_reactive_eb_3d_checkpoint,read_reactive_eb_3d_checkpoint
  implicit none
  private
  type, public :: eb_tree_checkpoint_3d
    real(dp) :: time=0._dp,minimum_dt=0._dp,minimum_theta=1._dp
    integer :: steps=0
    real(dp), allocatable :: initial(:),scale(:),boundary(:),wall(:),chemical(:),forcing(:)
    type(amr_tree_state_3d) :: reaction_source,heat_release,internal_energy
    type(dual_energy_policy) :: energy_policy
    real(dp) :: energy_correction=0._dp,internal_energy_pressure_work=0._dp
    logical :: has_internal_energy_pressure_work=.false.
  end type
  public :: reactive_eb_tree_checkpoint_record_3d
  public :: checkpoint_held_restriction_matches_3d
contains
  ! Check immediate-child restriction, optionally followed by the saved
  ! dual-energy reset. This validates a copy; never repairs checkpoint data.
  logical function checkpoint_held_restriction_matches_3d(species,topology,node,geometry, &
      actual,held,restricted,restricted_held,temperature,policy,frame) result(matches)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: topology
    integer, intent(in) :: node
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: actual(:,:,:,:),held(:,:,:,:),restricted(:,:,:,:),restricted_held(:,:,:,:)
    real(dp), intent(in) :: temperature(:,:,:)
    type(dual_energy_policy), intent(in) :: policy
    type(rotating_frame_3d), optional, intent(in) :: frame
    logical :: covered(size(temperature,1),size(temperature,2),size(temperature,3)),valid
    real(dp) :: gas(size(actual,1)),h,t,change,phi,position(3)
    integer :: child,lo(3),hi(3),i,j,k
    matches=.false.;covered=.false.
    do child=node+1,size(topology%nodes)
      if(topology%nodes(child)%parent/=node) cycle
      lo=[topology%nodes(child)%patch%coarse_i_lower,topology%nodes(child)%patch%coarse_j_lower, &
        topology%nodes(child)%patch%coarse_k_lower]
      hi=[topology%nodes(child)%patch%coarse_i_upper,topology%nodes(child)%patch%coarse_j_upper, &
        topology%nodes(child)%patch%coarse_k_upper]
      covered(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
    enddo
    do k=1,size(temperature,3)
    do j=1,size(temperature,2)
    do i=1,size(temperature,1)
      gas=restricted(:,i,j,k);h=restricted_held(1,i,j,k)
      if(cell_matches()) cycle
      if(.not.covered(i,j,k).or.geometry%volume_fraction(i,j,k)<=0) return
      t=temperature(i,j,k);phi=0
      if(present(frame)) then
        position=[geometry%x_lower,geometry%y_lower,geometry%z_lower]+ &
          (real([i,j,k],dp)-.5_dp)*[geometry%dx,geometry%dy,geometry%dz]
        call rotating_potential_3d(frame,position,phi,valid)
        if(.not.valid) return
      endif
      change=0
      call reset_nasa7_dual_energy_cell(species,gas,h,t,policy,change,valid,phi)
      if(.not.valid) return
      if(.not.cell_matches()) return
    enddo
    enddo
    enddo
    matches=.true.
  contains
    logical function cell_matches()
      cell_matches=all(abs(gas-actual(:,i,j,k))<=2.e-12_dp*max(1._dp,abs(gas))).and. &
        abs(h-held(1,i,j,k))<=2.e-12_dp*max(1._dp,abs(h))
    end function
  end function checkpoint_held_restriction_matches_3d
  ! Caller owns the formatted stream and stages file publication. The tree
  ! topology/geometry must already match the expected restart configuration.
  ! On reads, no state or metadata publishes until the entire envelope passes.
  subroutine reactive_eb_tree_checkpoint_record_3d(unit,reading,species,reactions,transport,config, &
      tree,record,ok,message,bundle_sha256,chemistry_integrator,base_mole_fractions,element_names,element_atoms,frame)
    type(rotating_frame_3d), optional, intent(in) :: frame
    integer, intent(in) :: unit
    logical, intent(in) :: reading
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_tree_checkpoint_3d), intent(inout) :: record
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    character(len=*), optional, intent(in) :: bundle_sha256,chemistry_integrator
    real(dp), optional, intent(in) :: base_mole_fractions(:),element_atoms(:,:)
    character(len=*), optional, intent(in) :: element_names(:)
    type(eb_reactive_tree_3d) :: candidate,averaged
    type(eb_tree_checkpoint_3d) :: saved
    type(reactive_eb_3d_config) :: cfg
    type(amr_tree_state_3d) :: averaged_held
    type(dual_energy_policy), allocatable :: patch_policy
    real(dp), allocatable :: patch_held(:,:,:),patch_prefix
    real(dp), allocatable :: integral(:),l1(:),elements(:),expected(:),expected_l1(:),expected_elements(:),current(:)
    real(dp) :: got_time,got_dt,got_theta,diffusivity,stored_force(3)
    integer :: nv,ne,node,count,ios,descriptor(12),stored(12),got_step,i,j,k,n(3),controls(2)
    logical :: valid,indexed,sdc,with_heat,with_held,file_held,file_pressure_work
    integer :: held_flags(2),energy_flags(3)
    real(dp) :: stored_eta2
    integer(int64), allocatable :: patch_bytes(:)
    integer(int64) :: patch_start,patch_end
    character(len=64) :: label
    character(len=*), parameter :: magic='PELEF_EB_TREE_CHECKPOINT_V1',terminator='END_EB_TREE_CHECKPOINT'
    ok=.false.;message='Invalid EB tree checkpoint';indexed=.false.
    sdc=config%time_integrator=='sdc'
    if(.not.sdc.and.config%time_integrator/='split') return
    if(.not.all(ieee_is_finite(config%external_forcing))) return
    if(sdc.and.config%sdc_iterations<1) return
    if(.not.sdc.and.any(config%external_forcing/=0)) return
    nv=reactive_nvar(size(species));ne=3
    if (present(element_names)) ne=size(element_names)
    allocate(integral(nv),l1(nv),elements(ne),expected(nv),expected_l1(nv),expected_elements(ne),current(nv))
    call integrate_reactive_eb_tree_3d(tree,current,valid)
    if (.not.valid) return
    candidate=tree;saved=record
    with_held=allocated(record%internal_energy%fields)
    if(with_held) then
      if(.not.record%internal_energy%is_valid(tree%topology)) return
      do node=1,size(tree%nodes)
        if(size(record%internal_energy%fields(node)%values,1)/=1) return
      enddo
    endif
    with_heat=allocated(saved%heat_release%fields)
    count=size(tree%nodes)
    if (reading) then
      read(unit,'(a)',iostat=ios) label
      if (ios/=0) return
      file_pressure_work=label=='PELEF_EB_TREE_CHECKPOINT_V9'.or.label=='PELEF_EB_TREE_CHECKPOINT_V10'
      file_held=label=='PELEF_EB_TREE_CHECKPOINT_V7'.or.label=='PELEF_EB_TREE_CHECKPOINT_V8'.or.file_pressure_work
      saved%has_internal_energy_pressure_work=file_pressure_work
      saved%internal_energy_pressure_work=0
      if(file_held.neqv.with_held) return
      with_heat=label=='PELEF_EB_TREE_CHECKPOINT_V5'.or.label=='PELEF_EB_TREE_CHECKPOINT_V6'
      indexed=label=='PELEF_EB_TREE_CHECKPOINT_V2'.or.label=='PELEF_EB_TREE_CHECKPOINT_V4'.or. &
        label=='PELEF_EB_TREE_CHECKPOINT_V6'.or.label=='PELEF_EB_TREE_CHECKPOINT_V8'.or. &
        label=='PELEF_EB_TREE_CHECKPOINT_V10'
      if(with_held) then
        read(unit,*,iostat=ios) held_flags
        if(ios/=0.or.any(held_flags<0).or.any(held_flags>1)) return
        if(held_flags(1)/=merge(1,0,sdc)) return
        with_heat=held_flags(2)==1
        if(with_heat.and..not.sdc) return
      else if(sdc) then
        if(label/='PELEF_EB_TREE_CHECKPOINT_V3'.and.label/='PELEF_EB_TREE_CHECKPOINT_V4'.and..not.with_heat) return
      else
        if(label/=magic.and.label/='PELEF_EB_TREE_CHECKPOINT_V2') return
      end if
      read(unit,*,iostat=ios) got_step
      if (ios/=0) return
      if (got_step/=count) return
    else
      if (.not.valid_record()) return
      if(with_held) then
        write(unit,'(a)',iostat=ios) 'PELEF_EB_TREE_CHECKPOINT_V9'
        if(ios/=0) return
        write(unit,*,iostat=ios) merge(1,0,sdc),merge(1,0,with_heat)
      else if(sdc) then
        if(with_heat) then
          write(unit,'(a)',iostat=ios) 'PELEF_EB_TREE_CHECKPOINT_V5'
        else
          write(unit,'(a)',iostat=ios) 'PELEF_EB_TREE_CHECKPOINT_V3'
        endif
      else
        write(unit,'(a)',iostat=ios) magic
      end if
      if (ios/=0) return
      write(unit,*,iostat=ios) count
      if (ios/=0) return
    end if
    do node=1,count
      associate(n=>tree%topology%nodes(node),p=>tree%topology%nodes(node)%patch)
        descriptor=[n%parent,n%level,n%shape,p%refinement_ratio,p%coarse_i_lower,p%coarse_i_upper, &
          p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper]
      end associate
      if (reading) then
        read(unit,*,iostat=ios) stored
        if (ios/=0) return
        if (any(stored/=descriptor)) return
      else
        write(unit,*,iostat=ios) descriptor
        if (ios/=0) return
      end if
    end do
    if (reading) then
      if (allocated(saved%initial)) deallocate(saved%initial)
      if (allocated(saved%scale)) deallocate(saved%scale)
      if (allocated(saved%boundary)) deallocate(saved%boundary)
      if (allocated(saved%wall)) deallocate(saved%wall)
      if (allocated(saved%chemical)) deallocate(saved%chemical)
      allocate(saved%initial(nv),saved%scale(nv),saved%boundary(nv),saved%wall(nv),saved%chemical(nv))
      read(unit,*,iostat=ios) saved%time,saved%steps,saved%minimum_dt,saved%minimum_theta
      if (ios/=0) return
      read(unit,*,iostat=ios) saved%initial,saved%scale,saved%boundary,saved%wall,saved%chemical
    else
      write(unit,'(es27.18e3,1x,i0,2(1x,es27.18e3))',iostat=ios) &
        saved%time,saved%steps,saved%minimum_dt,saved%minimum_theta
      if (ios/=0) return
      write(unit,'(*(es27.18e3,1x))',iostat=ios) saved%initial,saved%scale,saved%boundary,saved%wall,saved%chemical
    end if
    if (ios/=0) return
    if(with_held) then
      if(reading) then
        if(file_pressure_work) then
          read(unit,*,iostat=ios) energy_flags,stored_eta2,saved%energy_correction, &
          saved%internal_energy_pressure_work
        else
          ! Previous held formats have no pressure-work history. The caller
          ! may resume MOL (zero pdivu); Godunov must require the new envelope.
          read(unit,*,iostat=ios) energy_flags,stored_eta2,saved%energy_correction
        endif
        if(ios/=0) return
        if(.not.all(ieee_is_finite([stored_eta2,saved%energy_correction,record%energy_policy%eta2]))) return
        if(any(energy_flags/=[merge(1,0,record%energy_policy%allow_small_energy), &
          merge(1,0,record%energy_policy%allow_negative_energy),merge(1,0,record%energy_policy%update_total_from_internal)])) return
        if(stored_eta2/=record%energy_policy%eta2) return
      else
        write(unit,'(3(i0,1x),3(es27.18e3,1x))',iostat=ios) &
          merge(1,0,saved%energy_policy%allow_small_energy),merge(1,0,saved%energy_policy%allow_negative_energy), &
          merge(1,0,saved%energy_policy%update_total_from_internal),saved%energy_policy%eta2,saved%energy_correction, &
          saved%internal_energy_pressure_work
        if(ios/=0) return
      endif
      if(.not.ieee_is_finite(saved%internal_energy_pressure_work)) return
      if(abs(saved%internal_energy_pressure_work)>huge(1._dp)/8) return
      if(abs(saved%energy_correction)>huge(1._dp)/8) return
      allocate(patch_policy,patch_prefix);patch_policy=saved%energy_policy;patch_prefix=saved%energy_correction
    else
      saved%energy_correction=0;saved%internal_energy_pressure_work=0
    endif
    if(sdc) then
      if(reading) then
        read(unit,*,iostat=ios) controls
        if(ios/=0) return
        if(any(controls/=[config%sdc_iterations, &
          merge(1,0,config%chemistry_enabled)+merge(2,0,config%sdc_reactor=="cvode")+ &
          merge(4,0,config%sdc_reactor=="cvode_full")])) return
        read(unit,*,iostat=ios) stored_force
        if(ios/=0) return
        if(.not.all(ieee_is_finite(stored_force))) return
        if(any(stored_force/=config%external_forcing)) return
        if(allocated(saved%forcing)) deallocate(saved%forcing)
        allocate(saved%forcing(nv))
        read(unit,*,iostat=ios) saved%forcing
        if(allocated(saved%reaction_source%fields)) deallocate(saved%reaction_source%fields)
        if(config%chemistry_enabled) then
          allocate(saved%reaction_source%fields(count))
          do node=1,count
            n=tree%topology%nodes(node)%shape
            allocate(saved%reaction_source%fields(node)%values(nv-4,n(1),n(2),n(3)))
          end do
        end if
      else
        write(unit,*,iostat=ios) config%sdc_iterations, &
          merge(1,0,config%chemistry_enabled)+merge(2,0,config%sdc_reactor=="cvode")+ &
          merge(4,0,config%sdc_reactor=="cvode_full")
        if(ios/=0) return
        write(unit,'(*(es27.18e3,1x))',iostat=ios) config%external_forcing
        if(ios/=0) return
        write(unit,'(*(es27.18e3,1x))',iostat=ios) saved%forcing
      end if
      if(ios/=0) return
    end if
    if(reading) then
      if(allocated(saved%heat_release%fields)) deallocate(saved%heat_release%fields)
      if(with_heat) then
        allocate(saved%heat_release%fields(count))
        do node=1,count
          n=tree%topology%nodes(node)%shape
          allocate(saved%heat_release%fields(node)%values(1,n(1),n(2),n(3)))
        enddo
      endif
    endif
    if (indexed) then
      allocate(patch_bytes(count))
      do node=1,count
        read(unit,*,iostat=ios) patch_bytes(node)
        if (ios/=0) return
        if (patch_bytes(node)<=0) return
      end do
    end if
    do node=1,count
      cfg=config
      if(with_held) patch_held=saved%internal_energy%fields(node)%values(1,:,:,:)
      ! SDC controls and sources are bound by the outer V3/V4 envelope. Embedded
      ! gas records retain their legacy format for reuse of mechanism checks.
      if(sdc) then
        cfg%time_integrator='split';cfg%external_forcing=0;cfg%sdc_reactor='native'
      end if
      associate(g=>candidate%nodes(node)%geometry,u=>candidate%nodes(node)%state,t=>candidate%nodes(node)%temperature)
      cfg%nx=g%nx;cfg%ny=g%ny;cfg%nz=g%nz
      cfg%x_lower=g%x_lower;cfg%x_upper=g%x_upper
      cfg%y_lower=g%y_lower;cfg%y_upper=g%y_upper
      cfg%z_lower=g%z_lower;cfg%z_upper=g%z_upper
      if (reading) then
        if (indexed) then
          inquire(unit,pos=patch_start,iostat=ios)
          if (ios/=0.or.patch_start<1) return
        end if
        call read_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,g,u,t,got_time,got_step, &
          integral,l1,elements,got_dt,diffusivity,got_theta,valid,message,bundle_sha256,chemistry_integrator, &
          base_mole_fractions,unit,element_names,element_atoms,tree_patch=.true.,frame=frame, &
          internal_energy_density=patch_held,energy_correction=patch_prefix,energy_policy=patch_policy)
        if (.not.valid) return
        if(with_held) then
          ! Embedded copies bind each patch to the composite prefix; they
          ! are compared with the envelope and never added per patch.
          if(patch_prefix/=saved%energy_correction) return
          saved%internal_energy%fields(node)%values(1,:,:,:)=patch_held
        endif
        if(indexed.and.sdc.and.config%chemistry_enabled) then
          read(unit,'(a)',iostat=ios) label
          if(ios/=0.or.label/='SDC_REACTION_SOURCE') return
          n=tree%topology%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            read(unit,*,iostat=ios) saved%reaction_source%fields(node)%values(:,i,j,k)
            if(ios/=0) return
          end do
          end do
          end do
        end if
        if(indexed.and.with_heat) then
          read(unit,'(a)',iostat=ios) label
          if(ios/=0.or.label/='SDC_HEAT_RELEASE') return
          n=tree%topology%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            read(unit,*,iostat=ios) saved%heat_release%fields(node)%values(1,i,j,k)
            if(ios/=0) return
          enddo
          enddo
          enddo
        endif
        if (indexed) then
          inquire(unit,pos=patch_end,iostat=ios)
          if (ios/=0) return
          if (patch_end-patch_start/=patch_bytes(node)) return
        end if
        message='Inconsistent EB tree patch metadata'
        if (got_time/=saved%time.or.got_step/=saved%steps.or.got_dt/=saved%minimum_dt.or. &
          got_theta/=saved%minimum_theta.or.diffusivity/=0._dp) return
      end if
      call reactive_eb_integrals_3d(u,g,expected,valid,expected_l1)
      if (.not.valid) return
      if (present(element_atoms)) then
        call reactive_eb_composition_integrals_3d(species,u,g,element_atoms,expected_elements,valid)
      else
        call reactive_eb_element_integrals_3d(species,u,g,expected_elements,valid)
      end if
      if (.not.valid) return
      if (reading) then
        if (any(abs(integral-expected)>2.e-12_dp*max(1._dp,abs(expected))).or. &
          any(abs(l1-expected_l1)>2.e-12_dp*max(1._dp,expected_l1)).or. &
          any(abs(elements-expected_elements)>2.e-12_dp*max(1._dp,expected_elements))) return
      else
        call write_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,g,u,t,saved%time,saved%steps, &
          expected,expected_l1,expected_elements,saved%minimum_dt,0._dp,saved%minimum_theta,valid,message, &
          bundle_sha256,chemistry_integrator,base_mole_fractions,unit,element_names,element_atoms,tree_patch=.true.,frame=frame, &
          internal_energy_density=patch_held,energy_correction=patch_prefix,energy_policy=patch_policy)
        if (.not.valid) return
      end if
      end associate
    end do
    if(sdc.and.config%chemistry_enabled.and..not.indexed) then
      if(reading) then
        read(unit,'(a)',iostat=ios) label
        if(ios/=0) return
        if(label/='SDC_REACTION_SOURCES') return
      else
        write(unit,'(a)',iostat=ios) 'SDC_REACTION_SOURCES'
        if(ios/=0) return
      end if
      do node=1,count
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(reading) then
            read(unit,*,iostat=ios) saved%reaction_source%fields(node)%values(:,i,j,k)
          else
            write(unit,'(*(es27.18e3,1x))',iostat=ios) saved%reaction_source%fields(node)%values(:,i,j,k)
          end if
          if(ios/=0) return
        end do
        end do
        end do
      end do
    end if
    if(with_heat.and..not.indexed) then
      if(reading) then
        read(unit,'(a)',iostat=ios) label
        if(ios/=0.or.label/='SDC_HEAT_RELEASES') return
      else
        write(unit,'(a)',iostat=ios) 'SDC_HEAT_RELEASES'
        if(ios/=0) return
      endif
      do node=1,count
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(reading) then
            read(unit,*,iostat=ios) saved%heat_release%fields(node)%values(1,i,j,k)
          else
            write(unit,'(es27.18e3)',iostat=ios) saved%heat_release%fields(node)%values(1,i,j,k)
          endif
          if(ios/=0) return
        enddo
        enddo
        enddo
      enddo
    endif
    message='Invalid EB tree checkpoint envelope or composite inventory'
    if (reading) then
      read(unit,'(a)',iostat=ios) label
      if (ios/=0) return
      if (label/=terminator) return
    end if
    call integrate_reactive_eb_tree_3d(candidate,current,valid)
    if (.not.valid.or..not.valid_record()) return
    averaged=candidate
    if(with_held) then
      averaged_held=saved%internal_energy
      call average_down_reactive_eb_tree_3d(species,averaged,valid,frame,averaged_held,direct_children=.true.)
    else
      call average_down_reactive_eb_tree_3d(species,averaged,valid,frame)
    endif
    if (.not.valid) return
    do node=1,count
      if(with_held) then
        if(.not.checkpoint_held_restriction_matches_3d(species,candidate%topology,node, &
          candidate%nodes(node)%geometry,candidate%nodes(node)%state,saved%internal_energy%fields(node)%values, &
          averaged%nodes(node)%state,averaged_held%fields(node)%values,candidate%nodes(node)%temperature, &
          saved%energy_policy,frame)) return
        cycle
      endif
      if (any(abs(averaged%nodes(node)%state-candidate%nodes(node)%state)> &
        2.e-12_dp*max(1._dp,abs(averaged%nodes(node)%state)))) return
      if(with_held) then
        if(any(abs(averaged_held%fields(node)%values-saved%internal_energy%fields(node)%values)> &
          2.e-12_dp*max(1._dp,abs(averaged_held%fields(node)%values)))) return
      endif
    end do
    if (reading) then
      tree=candidate;record=saved
    else
      write(unit,'(a)',iostat=ios) terminator
      if (ios/=0) return
    end if
    message='';ok=.true.
  contains
    logical function valid_record() result(good)
      real(dp) :: forcing(nv)
      integer :: p
      good=.false.;forcing=0
      if(.not.ieee_is_finite(saved%internal_energy_pressure_work)) return
      if(abs(saved%internal_energy_pressure_work)>huge(1._dp)/8) return
      if(.not.with_held.and.saved%internal_energy_pressure_work/=0) return
      if(with_held) then
        if(.not.saved%internal_energy%is_valid(tree%topology)) return
        if(.not.all(ieee_is_finite([saved%energy_policy%eta2,saved%energy_correction]))) return
        if(abs(saved%energy_correction)>huge(1._dp)/8) return
        do p=1,count
          if(size(saved%internal_energy%fields(p)%values,1)/=1) return
          if(.not.all(ieee_is_finite(saved%internal_energy%fields(p)%values))) return
        enddo
      endif
      if(allocated(saved%heat_release%fields)) then
        if(.not.sdc.or..not.saved%heat_release%is_valid(tree%topology)) return
        if(size(saved%heat_release%fields(1)%values,1)/=1) return
      endif
      if(sdc) then
        if(.not.allocated(saved%forcing)) return
        if(size(saved%forcing)/=nv) return
        if(.not.all(ieee_is_finite(saved%forcing))) return
        forcing=saved%forcing
        if(config%chemistry_enabled) then
          if(.not.saved%reaction_source%is_valid(tree%topology)) return
          if(size(saved%reaction_source%fields(1)%values,1)/=nv-4) return
        else
          if(allocated(saved%reaction_source%fields)) return
        end if
      else
        if(allocated(saved%forcing).or.allocated(saved%reaction_source%fields)) return
      end if
      if(with_held) then
        if(forcing(5)>0.and.saved%energy_correction>0) then
          if(forcing(5)>huge(1._dp)-saved%energy_correction) return
        else if(forcing(5)<0.and.saved%energy_correction<0) then
          if(forcing(5)<-huge(1._dp)-saved%energy_correction) return
        endif
        forcing(5)=forcing(5)+saved%energy_correction
      endif
      if (.not.allocated(saved%initial).or..not.allocated(saved%scale).or..not.allocated(saved%boundary).or. &
        .not.allocated(saved%wall).or..not.allocated(saved%chemical)) return
      if (any([size(saved%initial),size(saved%scale),size(saved%boundary),size(saved%wall),size(saved%chemical)]/=nv)) return
      if (.not.all(ieee_is_finite([saved%time,saved%minimum_dt,saved%minimum_theta,saved%initial,saved%scale, &
        saved%boundary,saved%wall,saved%chemical]))) return
      if (saved%time<0._dp.or.saved%steps<0.or.saved%minimum_dt<0._dp.or. &
        saved%minimum_theta<0._dp.or.saved%minimum_theta>1._dp.or.any(saved%scale<0._dp)) return
      if (any(abs(current-saved%initial-saved%boundary-saved%wall-saved%chemical-forcing)> &
        5.e-10_dp*max(1._dp,saved%scale,abs(current),abs(saved%boundary)+abs(saved%wall)+abs(saved%chemical)+abs(forcing)))) return
      good=.true.
    end function valid_record
  end subroutine reactive_eb_tree_checkpoint_record_3d
end module eb_reactive_tree_checkpoint_3d_mod
