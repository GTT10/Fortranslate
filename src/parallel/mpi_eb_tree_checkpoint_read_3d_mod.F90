module mpi_eb_tree_checkpoint_read_3d_mod
  use dual_energy_mod, only: dual_energy_policy
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: iostat_end
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use eb_reactive_tree_checkpoint_3d_mod, only: eb_tree_checkpoint_3d,checkpoint_held_restriction_matches_3d
  use reactive_eb_3d_checkpoint_mod, only: read_reactive_eb_3d_checkpoint
  use reactive_eb_3d_driver_mod, only: reactive_eb_integrals_3d,reactive_eb_element_integrals_3d, &
    reactive_eb_composition_integrals_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,integrate_mpi_eb_tree_3d
  implicit none
  private
  public :: read_mpi_eb_tree_checkpoint_3d
contains
  ! Borrow a formatted stream positioned after the application policy header.
  ! V2/V4 lengths give direct patch seeks; V1/V3 require a root boundary scan.
  subroutine read_mpi_eb_tree_checkpoint_3d(path,unit,species,reactions,transport,config,tree,owners, &
      geometry,state,temperature,record,comm,ok,bundle_sha256,chemistry_integrator,base_mole_fractions, &
      element_names,element_atoms,frame)
    type(rotating_frame_3d), optional, intent(in) :: frame
    character(len=*), intent(in) :: path
    integer, intent(in) :: unit,owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(reactive_eb_3d_config), intent(in) :: config
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    type(eb_tree_checkpoint_3d), intent(inout) :: record
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: bundle_sha256,chemistry_integrator,element_names(:)
    real(dp), optional, intent(in) :: base_mole_fractions(:),element_atoms(:,:)
    type(amr_tree_state_3d) :: candidate,ct,averaged,at,averaged_held
    type(dual_energy_policy), allocatable :: patch_policy
    real(dp), allocatable :: patch_held(:,:,:),patch_prefix
    type(eb_tree_checkpoint_3d) :: saved
    type(reactive_eb_3d_config) :: cfg
    type(eb_geometry_3d) :: g
    integer(MPI_OFFSET_KIND), allocatable :: starts(:),ends(:),source_starts(:),heat_starts(:)
    integer(MPI_OFFSET_KIND) :: position,file_bytes,patch_bytes
    real(dp), allocatable :: current(:),integral(:),l1(:),elements(:),expected(:),expected_l1(:),expected_elements(:)
    real(dp) :: got_time,got_dt,got_theta,diffusivity,stored_force(3)
    integer :: rank,status,ios,node,nv,ne,count,descriptor(12),stored(12),got_steps,controls(2),i,j,k,n(3)
    character(len=64) :: label
    character(len=512) :: message
    integer :: redist_low,redist_high
    logical :: valid,indexed,sdc,with_heat,with_held,file_held,file_pressure_work
    integer :: held_flags(2),energy_flags(3)
    real(dp) :: stored_eta2
    ok=.false.;nv=reactive_nvar(size(species));ne=3
    valid=config%state_redist_max_order==0.or.config%state_redist_max_order==2.or.config%state_redist_max_order==4
    call agree(valid)
    if(.not.valid) return
    call MPI_Allreduce(config%state_redist_max_order,redist_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(config%state_redist_max_order,redist_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.redist_low/=redist_high) return
    sdc=config%time_integrator=='sdc'
    with_held=allocated(record%internal_energy%fields)
    saved%energy_policy=record%energy_policy
    if (present(element_names)) ne=size(element_names)
    candidate=state;ct=temperature
    if(with_held) then
      saved%internal_energy=record%internal_energy;averaged_held=record%internal_energy
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,valid,frame,averaged_held)
    else
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,candidate,ct,comm,valid,frame)
    endif
    if (.not.valid) return
    candidate=state;ct=temperature
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call same_string(trim(path),valid)
    if (.not.valid) return
    allocate(current(nv),integral(nv),l1(nv),elements(ne),expected(nv),expected_l1(nv),expected_elements(ne), &
      saved%initial(nv),saved%scale(nv),saved%boundary(nv),saved%wall(nv),saved%chemical(nv), &
      starts(size(owners)),ends(size(owners)),source_starts(size(owners)),heat_starts(size(owners)))
    read(unit,'(a)',iostat=ios) label
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    call same_string(label,valid)
    if(.not.valid) return
    file_pressure_work=label=='PELEF_EB_TREE_CHECKPOINT_V9'.or.label=='PELEF_EB_TREE_CHECKPOINT_V10'
    file_held=label=='PELEF_EB_TREE_CHECKPOINT_V7'.or.label=='PELEF_EB_TREE_CHECKPOINT_V8'.or.file_pressure_work
    saved%has_internal_energy_pressure_work=file_pressure_work
    saved%internal_energy_pressure_work=0
    valid=file_held.eqv.with_held
    call agree(valid)
    if(.not.valid) return
    with_heat=label=='PELEF_EB_TREE_CHECKPOINT_V5'.or.label=='PELEF_EB_TREE_CHECKPOINT_V6'
    indexed=label=='PELEF_EB_TREE_CHECKPOINT_V2'.or.label=='PELEF_EB_TREE_CHECKPOINT_V4'.or. &
      label=='PELEF_EB_TREE_CHECKPOINT_V6'.or.label=='PELEF_EB_TREE_CHECKPOINT_V8'.or. &
      label=='PELEF_EB_TREE_CHECKPOINT_V10'
    if(with_held) then
      read(unit,*,iostat=ios) held_flags
      valid=ios==0
      if(valid) valid=all(held_flags>=0).and.all(held_flags<=1).and.held_flags(1)==merge(1,0,sdc)
      if(valid) then
        with_heat=held_flags(2)==1
        valid=sdc.or..not.with_heat
      endif
    else if(sdc) then
      valid=ios==0.and.(label=='PELEF_EB_TREE_CHECKPOINT_V3'.or.label=='PELEF_EB_TREE_CHECKPOINT_V4'.or.with_heat)
    else
      valid=ios==0.and.config%time_integrator=='split'.and. &
        (label=='PELEF_EB_TREE_CHECKPOINT_V1'.or.label=='PELEF_EB_TREE_CHECKPOINT_V2')
    end if
    if (valid) then
      read(unit,*,iostat=ios) count
      valid=ios==0
      if (valid) valid=count==size(owners)
    end if
    do node=1,size(owners)
      if (.not.valid) exit
      associate(n=>tree%nodes(node),p=>tree%nodes(node)%patch)
      descriptor=[n%parent,n%level,n%shape,p%refinement_ratio,p%coarse_i_lower,p%coarse_i_upper, &
        p%coarse_j_lower,p%coarse_j_upper,p%coarse_k_lower,p%coarse_k_upper]
      end associate
      read(unit,*,iostat=ios) stored
      valid=ios==0
      if (valid) valid=all(stored==descriptor)
    end do
    if (valid) then
      read(unit,*,iostat=ios) saved%time,saved%steps,saved%minimum_dt,saved%minimum_theta
      valid=ios==0
    end if
    if (valid) then
      read(unit,*,iostat=ios) saved%initial,saved%scale,saved%boundary,saved%wall,saved%chemical
      valid=ios==0
    end if
    if(valid.and.with_held) then
      if(file_pressure_work) then
        read(unit,*,iostat=ios) energy_flags,stored_eta2,saved%energy_correction,saved%internal_energy_pressure_work
      else
        read(unit,*,iostat=ios) energy_flags,stored_eta2,saved%energy_correction
      endif
      valid=ios==0
      if(valid) valid=all(ieee_is_finite([stored_eta2,saved%energy_correction, &
        record%energy_policy%eta2,saved%internal_energy_pressure_work]))
      if(valid) valid=all(energy_flags==[merge(1,0,record%energy_policy%allow_small_energy), &
        merge(1,0,record%energy_policy%allow_negative_energy),merge(1,0,record%energy_policy%update_total_from_internal)])
      if(valid) valid=stored_eta2==record%energy_policy%eta2.and.abs(saved%energy_correction)<=huge(1._dp)/8
      if(valid) valid=abs(saved%internal_energy_pressure_work)<=huge(1._dp)/8
      if(valid) then
        allocate(patch_policy,patch_prefix);patch_policy=record%energy_policy;patch_prefix=saved%energy_correction
      endif
    endif
    if(valid.and.sdc) then
      read(unit,*,iostat=ios) controls
      valid=ios==0
      if(valid) valid=all(controls==[config%sdc_iterations, &
          merge(1,0,config%chemistry_enabled)+merge(2,0,config%sdc_reactor=="cvode")+ &
          merge(4,0,config%sdc_reactor=="cvode_full")])
      if(valid) then
        read(unit,*,iostat=ios) stored_force
        valid=ios==0
        if(valid) valid=all(ieee_is_finite(stored_force))
        if(valid) valid=all(stored_force==config%external_forcing)
      end if
      if(valid) then
        allocate(saved%forcing(nv))
        read(unit,*,iostat=ios) saved%forcing
        valid=ios==0
        if(config%chemistry_enabled) allocate(saved%reaction_source%fields(size(owners)))
        if(with_heat) allocate(saved%heat_release%fields(size(owners)))
      end if
    end if
    if (valid) valid=valid_record(.false.)
    call agree(valid)
    if (.not.valid) return
    if (rank==0) then
      if (indexed) then
        do node=1,size(owners)
          read(unit,*,iostat=ios) ends(node)
          valid=ios==0
          if (valid) valid=ends(node)>0
          if (.not.valid) exit
        end do
        if (valid) then
          inquire(unit,pos=position,size=file_bytes,iostat=ios)
          valid=ios==0.and.position>=1.and.file_bytes<huge(file_bytes)
        end if
        do node=1,size(owners)
          if (.not.valid) exit
          patch_bytes=ends(node)
          valid=position<=file_bytes
          if (valid) valid=patch_bytes<=file_bytes-position
          if (.not.valid) exit
          starts(node)=position;position=position+patch_bytes;ends(node)=position
        end do
        if (valid) then
          read(unit,'()',pos=position,advance='no',iostat=ios)
          valid=ios==0
        end if
      else
        do node=1,size(owners)
          inquire(unit,pos=starts(node),iostat=ios)
          if (ios/=0) then
            valid=.false.;exit
          end if
          read(unit,'(a)',iostat=ios) label
          valid=ios==0.and.label=='PELEF_REACTIVE_EB_3D_CHECKPOINT'
          if (.not.valid) exit
          do
            read(unit,'(a)',iostat=ios) label
            if (ios/=0) exit
            if (label=='END_CHECKPOINT') exit
          end do
          valid=ios==0
          if (.not.valid) exit
          inquire(unit,pos=ends(node),iostat=ios)
          valid=ios==0
          if (.not.valid) exit
        end do
      end if
      if(valid.and.sdc.and..not.indexed.and.config%chemistry_enabled) then
        read(unit,'(a)',iostat=ios) label
        valid=ios==0.and.label=='SDC_REACTION_SOURCES'
        do node=1,size(owners)
          if(.not.valid) exit
          inquire(unit,pos=source_starts(node),iostat=ios)
          n=tree%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(ios==0) read(unit,'(a)',iostat=ios) label
          end do
          end do
          end do
          valid=ios==0
          if(.not.valid) exit
          inquire(unit,pos=ends(node),iostat=ios)
          valid=ios==0
        end do
      end if
      if(valid.and.with_heat.and..not.indexed) then
        read(unit,'(a)',iostat=ios) label
        valid=ios==0.and.label=='SDC_HEAT_RELEASES'
        do node=1,size(owners)
          if(.not.valid) exit
          inquire(unit,pos=heat_starts(node),iostat=ios)
          n=tree%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(ios==0) read(unit,'(a)',iostat=ios) label
          enddo
          enddo
          enddo
          valid=ios==0
          if(.not.valid) exit
          inquire(unit,pos=ends(node),iostat=ios)
          valid=ios==0
        enddo
      endif
      if (valid) then
        read(unit,'(a)',iostat=ios) label
        valid=ios==0.and.label=='END_EB_TREE_CHECKPOINT'
      end if
      if (valid) then
        read(unit,'(a)',iostat=ios) label
        valid=ios==iostat_end
      end if
    end if
    call agree(valid)
    if (.not.valid) return
    call MPI_Bcast(starts,size(owners),MPI_OFFSET,0,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Bcast(ends,size(owners),MPI_OFFSET,0,comm,status)
    if (status/=MPI_SUCCESS) return
    if(sdc.and..not.indexed.and.config%chemistry_enabled) then
      call MPI_Bcast(source_starts,size(owners),MPI_OFFSET,0,comm,status)
      if(status/=MPI_SUCCESS) return
    end if
    if(with_heat.and..not.indexed) then
      call MPI_Bcast(heat_starts,size(owners),MPI_OFFSET,0,comm,status)
      if(status/=MPI_SUCCESS) return
    endif
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if(with_held) patch_held=saved%internal_energy%fields(node)%values(1,:,:,:)
      read(unit,'()',pos=starts(node),advance='no',iostat=ios)
      valid=ios==0
      if (.not.valid) exit
      g=geometry(node);cfg=config
      if(sdc) then
        cfg%time_integrator='split';cfg%external_forcing=0;cfg%sdc_reactor='native'
      end if
      cfg%nx=g%nx;cfg%ny=g%ny;cfg%nz=g%nz
      cfg%x_lower=g%x_lower;cfg%x_upper=g%x_upper
      cfg%y_lower=g%y_lower;cfg%y_upper=g%y_upper
      cfg%z_lower=g%z_lower;cfg%z_upper=g%z_upper
      got_time=0._dp;got_steps=0;got_dt=0._dp;got_theta=1._dp;diffusivity=0._dp
      call read_reactive_eb_3d_checkpoint('',species,reactions,transport,cfg,g,candidate%fields(node)%values, &
        ct%fields(node)%values(1,:,:,:),got_time,got_steps,integral,l1,elements,got_dt,diffusivity,got_theta, &
        valid,message,bundle_sha256,chemistry_integrator,base_mole_fractions,unit,element_names,element_atoms, &
        tree_patch=.true.,frame=frame,internal_energy_density=patch_held,energy_correction=patch_prefix,energy_policy=patch_policy)
      if (.not.valid) exit
      if(with_held) then
        valid=patch_prefix==saved%energy_correction
        if(.not.valid) exit
        saved%internal_energy%fields(node)%values(1,:,:,:)=patch_held
      endif
      if(sdc.and.config%chemistry_enabled) then
        if(indexed) then
          read(unit,'(a)',iostat=ios) label
          valid=ios==0.and.label=='SDC_REACTION_SOURCE'
        else
          read(unit,'()',pos=source_starts(node),advance='no',iostat=ios)
          valid=ios==0
        end if
        if(.not.valid) exit
        n=tree%nodes(node)%shape
        allocate(saved%reaction_source%fields(node)%values(nv-4,n(1),n(2),n(3)))
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(ios==0) read(unit,*,iostat=ios) saved%reaction_source%fields(node)%values(:,i,j,k)
        end do
        end do
        end do
        valid=ios==0
        if(valid) valid=all(ieee_is_finite(saved%reaction_source%fields(node)%values))
        if(.not.valid) exit
      end if
      if(with_heat) then
        if(indexed) then
          read(unit,'(a)',iostat=ios) label
          valid=ios==0.and.label=='SDC_HEAT_RELEASE'
        else
          read(unit,'()',pos=heat_starts(node),advance='no',iostat=ios)
          valid=ios==0
        endif
        if(.not.valid) exit
        n=tree%nodes(node)%shape
        allocate(saved%heat_release%fields(node)%values(1,n(1),n(2),n(3)))
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(ios==0) read(unit,*,iostat=ios) saved%heat_release%fields(node)%values(1,i,j,k)
        enddo
        enddo
        enddo
        valid=ios==0
        if(valid) valid=all(ieee_is_finite(saved%heat_release%fields(node)%values))
        if(.not.valid) exit
      endif
      inquire(unit,pos=position,iostat=ios)
      valid=ios==0.and.position==ends(node).and.matching_eb_metrics_3d(g,geometry(node))
      if (.not.valid) exit
      valid=got_time==saved%time.and.got_steps==saved%steps.and.got_dt==saved%minimum_dt.and. &
        got_theta==saved%minimum_theta.and.diffusivity==0._dp
      if (.not.valid) exit
      call reactive_eb_integrals_3d(candidate%fields(node)%values,g,expected,valid,expected_l1)
      if (.not.valid) exit
      if (present(element_atoms)) then
        call reactive_eb_composition_integrals_3d(species,candidate%fields(node)%values,g,element_atoms,expected_elements,valid)
      else
        call reactive_eb_element_integrals_3d(species,candidate%fields(node)%values,g,expected_elements,valid)
      end if
      if (.not.valid) exit
      valid=all(abs(integral-expected)<=2.e-12_dp*max(1._dp,abs(expected))).and. &
        all(abs(l1-expected_l1)<=2.e-12_dp*max(1._dp,expected_l1)).and. &
        all(abs(elements-expected_elements)<=2.e-12_dp*max(1._dp,expected_elements))
      if (.not.valid) exit
    end do
    call agree(valid)
    if (.not.valid) return
    call integrate_mpi_eb_tree_3d(tree,owners,geometry,candidate,current,comm,valid)
    if (.not.valid) return
    valid=valid_record(.true.);call agree(valid)
    if (.not.valid) return
    averaged=candidate;at=ct
    if(with_held) then
      averaged_held=saved%internal_energy
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,averaged,at,comm,valid,frame,averaged_held, &
        direct_children=.true.)
    else
      call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,averaged,at,comm,valid,frame)
    endif
    if (.not.valid) return
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if(with_held) then
        if(.not.checkpoint_held_restriction_matches_3d(species,tree,node,geometry(node), &
          candidate%fields(node)%values,saved%internal_energy%fields(node)%values, &
          averaged%fields(node)%values,averaged_held%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
          saved%energy_policy,frame)) valid=.false.
        cycle
      endif
      if (any(abs(averaged%fields(node)%values-candidate%fields(node)%values)> &
        2.e-12_dp*max(1._dp,abs(averaged%fields(node)%values)))) valid=.false.
      if(with_held) then
        if(any(abs(averaged_held%fields(node)%values-saved%internal_energy%fields(node)%values)> &
          2.e-12_dp*max(1._dp,abs(averaged_held%fields(node)%values)))) valid=.false.
      endif
    end do
    call agree(valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    record=saved;ok=.true.
  contains
    logical function valid_record(check_inventory) result(good)
      logical, intent(in) :: check_inventory
      real(dp) :: forcing(nv)
      good=.false.;forcing=0
      if(.not.ieee_is_finite(saved%internal_energy_pressure_work)) return
      if(abs(saved%internal_energy_pressure_work)>huge(1._dp)/8) return
      if(.not.with_held.and.saved%internal_energy_pressure_work/=0) return
      if(with_held) then
        if(.not.all(ieee_is_finite([saved%energy_policy%eta2,saved%energy_correction]))) return
        if(abs(saved%energy_correction)>huge(1._dp)/8) return
      endif
      if(sdc) then
        if(.not.allocated(saved%forcing).or.config%sdc_iterations<1) return
        if(size(saved%forcing)/=nv) return
        if(.not.all(ieee_is_finite(saved%forcing))) return
        forcing=saved%forcing
      else
        if(any(config%external_forcing/=0)) return
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
      if (.not.check_inventory) then
        good=.true.;return
      end if
      if (any(abs(current-saved%initial-saved%boundary-saved%wall-saved%chemical-forcing)> &
        5.e-10_dp*max(1._dp,saved%scale,abs(current),abs(saved%boundary)+abs(saved%wall)+abs(saved%chemical)+abs(forcing)))) return
      good=.true.
    end function valid_record
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
    subroutine same_string(value,success)
      character(len=*), intent(in) :: value
      logical, intent(out) :: success
      character(len=:), allocatable :: reference
      integer :: minimum,maximum,ierr,allocation_status
      success=.false.
      call MPI_Allreduce(len(value),minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
      if (ierr/=MPI_SUCCESS) return
      call MPI_Allreduce(len(value),maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if (ierr/=MPI_SUCCESS.or.minimum/=maximum) return
      allocate(character(len=minimum)::reference,stat=allocation_status)
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      reference=value
      call MPI_Bcast(reference,minimum,MPI_CHARACTER,0,comm,ierr)
      success=ierr==MPI_SUCCESS.and.reference==value
      call agree(success)
    end subroutine
  end subroutine
end module mpi_eb_tree_checkpoint_read_3d_mod
