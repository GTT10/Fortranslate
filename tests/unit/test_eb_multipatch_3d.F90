program test_eb_multipatch_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,advance_reactive_eb_multipatch_euler_3d, &
    advance_reactive_eb_two_level_euler_3d,advance_reactive_eb_multipatch_transport_3d, &
    advance_reactive_eb_two_level_transport_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(eb_reactive_patch_3d) :: patches(2),saved(2),reordered(2)
  type(eb_geometry_3d) :: full,coarse,joined
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),sc(:,:,:,:),st(:,:,:),rc(:,:,:,:),rt(:,:,:)
  real(dp), allocatable :: primitive(:)
  real(dp), allocatable :: joined_u(:,:,:,:),joined_t(:,:,:)
  real(dp), allocatable :: boundary(:),wall(:),before(:),after(:),rb(:),rw(:),sb(:),sw(:)
  integer :: lower(3,2),extent(3,2),upper(3,2),bad_lower(3,2),operation,ratio,mode,p,nvar,nf(3),i,j,k,c,step,lo(3),hi(3)
  real(dp) :: amount,reference_amount,interval,theta,accepted_theta,saved_theta,reference_theta,recovered,speed
  logical :: valid,refined(4,4,4)
  call load_h2o2_elementary_thermo(species,valid)
  call require(valid,'thermo')
  call load_h2o2_elementary_transport(transport,valid)
  call require(valid,'transport')
  nvar=reactive_nvar(size(species))
  allocate(primitive(reactive_nprim(size(species))))
  allocate(uc(nvar,4,4,4),tc(4,4,4),boundary(nvar),wall(nvar),before(nvar),after(nvar),rb(nvar),rw(nvar))
  do operation=1,2
  interval=1.e-8_dp
  if (operation==2) interval=1.e-5_dp
  do ratio=2,3
    call build_axis_plane_eb_geometry_3d(4*ratio,4*ratio,4*ratio,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
      'z',.31_dp,full,valid)
    call require(valid,'full geometry')
    call coarsen_eb_geometry_3d(full,ratio,coarse,valid)
    call require(valid,'coarse geometry')
    do mode=1,3
      lower(:,1)=[2,2,2];lower(:,2)=[3,2,2];extent(:,1)=[1,2,2];extent(:,2)=[1,2,2]
      if (mode==2) then
        lower(:,2)=[3,3,2];extent(:,2)=[1,1,2]
      else if (mode==3) then
        lower(:,1)=[1,1,1];lower(:,2)=[3,3,3];extent=2
      end if
      upper=lower+extent-1;refined=.false.
      config%nx=4;config%ny=4;config%nz=4
      call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,valid)
      call require(valid,'coarse initialization')
      do p=1,2
        refined(lower(1,p):upper(1,p),lower(2,p):upper(2,p),lower(3,p):upper(3,p))=.true.
        call extract_eb_geometry_3d(full,(lower(:,p)-1)*ratio+1,upper(:,p)*ratio,patches(p)%geometry,valid)
        call require(valid,'patch extraction')
        nf=extent(:,p)*ratio
        if (allocated(patches(p)%state)) deallocate(patches(p)%state,patches(p)%temperature)
        allocate(patches(p)%state(nvar,nf(1),nf(2),nf(3)),patches(p)%temperature(nf(1),nf(2),nf(3)))
        config%nx=nf(1);config%ny=nf(2);config%nz=nf(3)
        call initialize_reactive_eb_density_sheet_3d(species,config,patches(p)%geometry, &
          patches(p)%state,patches(p)%temperature,valid)
        call require(valid,'patch initialization')
        do k=1,nf(3)
          do j=1,nf(2)
            do i=1,nf(1)
              lo=(lower(:,p)-1)*ratio+[i,j,k]
              patches(p)%state(:,i,j,k)=patches(p)%state(:,i,j,k)*(1._dp+.001_dp*(lo(1)+10*lo(2)+100*lo(3)))
              if (operation==2.and.mode==2.and.patches(p)%geometry%volume_fraction(i,j,k)>0._dp) then
                call reactive_conserved_to_primitive(species,patches(p)%state(:,i,j,k), &
                  patches(p)%temperature(i,j,k),primitive,recovered,speed,valid)
                call require(valid,'viscous fixture primitive state')
                primitive(2:4)=.01_dp*real([lo(2),lo(3),lo(1)],dp)**2
                call reactive_primitive_to_conserved(species,primitive,patches(p)%state(:,i,j,k), &
                  patches(p)%temperature(i,j,k),speed,valid)
                call require(valid,'transverse velocity fixture')
              end if
            end do
          end do
        end do
      end do
      boundary=17._dp;wall=19._dp
      do step=1,3
        call inventory(before)
        sc=uc;st=tc;saved=patches
        call advance(species,coarse,lower,ratio,'rusanov',interval, &
          uc,tc,patches,boundary,wall,valid)
        call require(valid,'multipatch update')
        accepted_theta=theta
        if (operation==2) then
          call require(theta>=0._dp.and.theta<=1._dp,'multipatch transport limiter range')
          call require(all(boundary==0._dp).and.all(wall==0._dp),'closed physical transport boundary')
        end if
        call require(maxval(abs(patches(1)%state-saved(1)%state))>1.e-12_dp,'fine evolution is nonzero')
        call inventory(after)
        call require(all(abs(after-before-boundary-wall)<=3.e-12_dp*max(1._dp,abs(before))), &
          'composite inventory equals boundary and wall changes')
        do p=1,2
          do k=lower(3,p),upper(3,p)
            do j=lower(2,p),upper(2,p)
              do i=lower(1,p),upper(1,p)
                lo=([i,j,k]-lower(:,p))*ratio+1;hi=lo+ratio-1
                do c=1,nvar
                  amount=uc(c,i,j,k)*coarse%volume_fraction(i,j,k)*coarse%dx*coarse%dy*coarse%dz
                  reference_amount=sum(patches(p)%state(c,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))* &
                    patches(p)%geometry%volume_fraction(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))* &
                    patches(p)%geometry%dx*patches(p)%geometry%dy*patches(p)%geometry%dz
                  call require(abs(amount-reference_amount)<=2.e-12_dp*max(1._dp,abs(reference_amount)), &
                    'all refined parents match fine inventory')
                end do
              end do
            end do
          end do
        end do
        if (step==1) then
          rc=sc;rt=st;reordered(1)=saved(2);reordered(2)=saved(1)
          call advance(species,coarse,lower(:,2:1:-1),ratio,'rusanov',interval, &
            rc,rt,reordered,rb,rw,valid)
          call require(theta==accepted_theta,'patch permutation leaves limiter exact')
          call require(valid.and.all(rc==uc).and.all(rt==tc).and.all(rb==boundary).and.all(rw==wall), &
            'patch permutation leaves coarse state and ledgers exact')
          do p=1,2
            call require(all(reordered(3-p)%state==patches(p)%state).and. &
              all(reordered(3-p)%temperature==patches(p)%temperature),'patch permutation leaves fine state exact')
          end do
        end if
        if (step==1.and.mode==1) then
          ! A rectangular refined region should evolve the same way whether
          ! stored as one patch or split across this same-level interface.
          call extract_eb_geometry_3d(full,[ratio+1,ratio+1,ratio+1],[3*ratio,3*ratio,3*ratio],joined,valid)
          call require(valid,'joined reference geometry')
          allocate(joined_u(nvar,2*ratio,2*ratio,2*ratio),joined_t(2*ratio,2*ratio,2*ratio))
          joined_u(:,1:ratio,:,:)=saved(1)%state;joined_u(:,ratio+1:,:,:)=saved(2)%state
          joined_t(1:ratio,:,:)=saved(1)%temperature;joined_t(ratio+1:,:,:)=saved(2)%temperature
          rc=sc;rt=st
          if (operation==1) then
            call advance_reactive_eb_two_level_euler_3d(species,coarse,joined,[2,2,2],ratio,'rusanov',interval, &
              rc,rt,joined_u,joined_t,rb,rw,valid)
          else
            call advance_reactive_eb_two_level_transport_3d(species,transport,coarse,joined,[2,2,2],ratio,interval, &
              .true.,.true.,.true.,.true.,rc,rt,joined_u,joined_t,rb,reference_theta,valid)
            rw=0._dp
          end if
          call require(valid,'one-patch reference step')
          call require(all(abs(rc-uc)<=3.e-12_dp*max(1._dp,abs(rc))).and. &
            all(abs(rt-tc)<=3.e-12_dp*max(1._dp,abs(rt))),'split and joined coarse solutions agree')
          call require(all(abs(rb-boundary)<=3.e-12_dp*max(1._dp,abs(rb))).and. &
            all(abs(rw-wall)<=3.e-12_dp*max(1._dp,abs(rw))),'split and joined ledgers agree')
          do p=1,2
            call require(all(abs(joined_u(:,(p-1)*ratio+1:p*ratio,:,:)-patches(p)%state)<= &
              3.e-12_dp*max(1._dp,abs(patches(p)%state))),'split and joined fine states agree')
          end do
          deallocate(joined_u,joined_t)
        end if
      end do
      sc=uc;st=tc;saved=patches;sb=boundary;sw=wall;saved_theta=theta
      bad_lower=lower;bad_lower(:,2)=lower(:,1)
      call advance(species,coarse,bad_lower,ratio,'rusanov',interval, &
        uc,tc,patches,boundary,wall,valid)
      call require(.not.valid,'overlapping patches rejected')
      call unchanged()
      call advance(species,coarse,lower,ratio,'rusanov',1.e8_dp, &
        uc,tc,patches,boundary,wall,valid)
      call require(.not.valid,'nonphysical multipatch interval rejected')
      call unchanged()
    end do
  end do
  end do
  write(*,'(a)') 'test_eb_multipatch_3d: PASS'
contains
  subroutine advance(gas,geometry,boxes,refinement,solver,dt,u,t,children,b,w,success)
    type(nasa7_species), intent(in) :: gas(:)
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: boxes(:,:),refinement
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: u(:,:,:,:),t(:,:,:),b(:),w(:)
    type(eb_reactive_patch_3d), intent(inout) :: children(:)
    logical, intent(out) :: success
    if (operation==1) then
      call advance_reactive_eb_multipatch_euler_3d(gas,geometry,boxes,refinement,solver,dt,u,t,children,b,w,success)
      if (success) theta=1._dp
    else
      call advance_reactive_eb_multipatch_transport_3d(gas,transport,geometry,boxes,refinement,dt, &
        .true.,.true.,.true.,.true.,u,t,children,b,theta,success)
      if (success) w=0._dp
    end if
  end subroutine advance
  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: component,patch_index
    do component=1,nvar
      values(component)=coarse%dx*coarse%dy*coarse%dz* &
        sum(uc(component,:,:,:)*coarse%volume_fraction,mask=.not.refined)
      do patch_index=1,2
        values(component)=values(component)+patches(patch_index)%geometry%dx*patches(patch_index)%geometry%dy* &
          patches(patch_index)%geometry%dz*sum(patches(patch_index)%state(component,:,:,:)* &
          patches(patch_index)%geometry%volume_fraction)
      end do
    end do
  end subroutine inventory
  subroutine unchanged()
    integer :: patch_index
    call require(theta==saved_theta,'limiter rollback')
    call require(all(uc==sc).and.all(tc==st).and.all(boundary==sb).and.all(wall==sw),'coarse and ledger rollback')
    do patch_index=1,2
      call require(all(patches(patch_index)%state==saved(patch_index)%state).and. &
        all(patches(patch_index)%temperature==saved(patch_index)%temperature),'all fine patches roll back')
    end do
  end subroutine unchanged
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_multipatch_3d
