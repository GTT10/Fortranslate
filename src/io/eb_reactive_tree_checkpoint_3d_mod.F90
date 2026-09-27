module eb_reactive_tree_checkpoint_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_tree_state_3d
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
    type(amr_tree_state_3d) :: reaction_source
  end type
  public :: reactive_eb_tree_checkpoint_record_3d
contains
  ! Caller owns the formatted stream and stages file publication. The tree
  ! topology/geometry must already match the expected restart configuration.
  ! On reads, no state or metadata publishes until the entire envelope passes.
  subroutine reactive_eb_tree_checkpoint_record_3d(unit,reading,species,reactions,transport,config, &
      tree,record,ok,message,bundle_sha256,chemistry_integrator,base_mole_fractions,element_names,element_atoms)
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
    real(dp), allocatable :: integral(:),l1(:),elements(:),expected(:),expected_l1(:),expected_elements(:),current(:)
    real(dp) :: got_time,got_dt,got_theta,diffusivity,stored_force(3)
    integer :: nv,ne,node,count,ios,descriptor(12),stored(12),got_step,i,j,k,n(3),controls(2)
    logical :: valid,indexed,sdc
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
    count=size(tree%nodes)
    if (reading) then
      read(unit,'(a)',iostat=ios) label
      if (ios/=0) return
      indexed=label=='PELEF_EB_TREE_CHECKPOINT_V2'.or.label=='PELEF_EB_TREE_CHECKPOINT_V4'
      if(sdc) then
        if(label/='PELEF_EB_TREE_CHECKPOINT_V3'.and.label/='PELEF_EB_TREE_CHECKPOINT_V4') return
      else
        if(label/=magic.and.label/='PELEF_EB_TREE_CHECKPOINT_V2') return
      end if
      read(unit,*,iostat=ios) got_step
      if (ios/=0) return
      if (got_step/=count) return
    else
      if (.not.valid_record()) return
      if(sdc) then
        write(unit,'(a)',iostat=ios) 'PELEF_EB_TREE_CHECKPOINT_V3'
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
    if(sdc) then
      if(reading) then
        read(unit,*,iostat=ios) controls
        if(ios/=0) return
        if(any(controls/=[config%sdc_iterations, &
          merge(1,0,config%chemistry_enabled)+merge(2,0,config%sdc_reactor=="cvode")])) return
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
          merge(1,0,config%chemistry_enabled)+merge(2,0,config%sdc_reactor=="cvode")
        if(ios/=0) return
        write(unit,'(*(es27.18e3,1x))',iostat=ios) config%external_forcing
        if(ios/=0) return
        write(unit,'(*(es27.18e3,1x))',iostat=ios) saved%forcing
      end if
      if(ios/=0) return
    end if
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
          base_mole_fractions,unit,element_names,element_atoms,tree_patch=.true.)
        if (.not.valid) return
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
          bundle_sha256,chemistry_integrator,base_mole_fractions,unit,element_names,element_atoms,tree_patch=.true.)
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
    message='Invalid EB tree checkpoint envelope or composite inventory'
    if (reading) then
      read(unit,'(a)',iostat=ios) label
      if (ios/=0) return
      if (label/=terminator) return
    end if
    call integrate_reactive_eb_tree_3d(candidate,current,valid)
    if (.not.valid.or..not.valid_record()) return
    averaged=candidate
    call average_down_reactive_eb_tree_3d(species,averaged,valid)
    if (.not.valid) return
    do node=1,count
      if (any(abs(averaged%nodes(node)%state-candidate%nodes(node)%state)> &
        2.e-12_dp*max(1._dp,abs(averaged%nodes(node)%state)))) return
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
      good=.false.;forcing=0
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
