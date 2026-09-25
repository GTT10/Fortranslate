program test_eb_patch_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: fill_reactive_eb_patch_halo_3d,overlay_reactive_eb_peer_halo_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: full,coarse,patch,saved_geometry
  type(reactive_eb_3d_config) :: config
  type(eb_integrated_flux_3d) :: flux
  real(dp), allocatable :: parent_u(:,:,:,:),parent_t(:,:,:),u(:,:,:,:),t(:,:,:),next(:,:,:,:),nt(:,:,:)
  real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),saved(:,:,:,:),saved_t(:,:,:)
  real(dp) :: area,expected_t
  integer :: nvar,i,j,k,axis,side,idx(3),inside(3),donor(3),lo(3),hi(3),first(3),last(3),mode
  logical :: ok,external
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nvar=reactive_nvar(size(species))
  allocate(parent_u(nvar,6,6,6),parent_t(6,6,6),u(nvar,6,6,6),t(6,6,6),next(nvar,6,6,6),nt(6,6,6))
  allocate(ghost(nvar,0:7,0:7,0:7),gt(0:7,0:7,0:7),saved(nvar,0:7,0:7,0:7),saved_t(0:7,0:7,0:7))
  call build_axis_plane_eb_geometry_3d(12,12,12,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',.17_dp,full,ok)
  call require(ok,'full fine geometry')
  call coarsen_eb_geometry_3d(full,2,coarse,ok)
  call require(ok,'coarse geometry')
  config%nx=6;config%ny=6;config%nz=6
  call initialize_reactive_eb_density_sheet_3d(species,config,coarse,parent_u,parent_t,ok)
  call require(ok,'parent state')
  do k=1,6
    do j=1,6
      do i=1,6
        parent_u(:,i,j,k)=parent_u(:,i,j,k)*(1._dp+.001_dp*(i+10*j+100*k))
      end do
    end do
  end do
  do mode=1,2
    lo=2
    if(mode==2) lo=1
    hi=lo+2
    call extract_eb_geometry_3d(full,(lo-1)*2+1,hi*2,patch,ok)
    call require(ok.and.patch%is_valid(),'patch extraction')
    call require(all(patch%volume_fraction==full%volume_fraction((lo(1)-1)*2+1:hi(1)*2, &
      (lo(2)-1)*2+1:hi(2)*2,(lo(3)-1)*2+1:hi(3)*2)),'exact volume crop')
    call require(all(patch%x_face_fraction==full%x_face_fraction((lo(1)-1)*2:hi(1)*2, &
      (lo(2)-1)*2+1:hi(2)*2,(lo(3)-1)*2+1:hi(3)*2)),'normal face bounds crop')
    call initialize_reactive_eb_density_sheet_3d(species,config,patch,u,t,ok)
    call require(ok,'patch state')
    call fill_reactive_eb_patch_halo_3d(species,coarse,patch,lo,2,parent_u,parent_t,u,t,ghost,gt,ok)
    call require(ok,'parent PCM halo')
    do axis=1,3
      do side=0,1
        first=1;last=6;first(axis)=7*side;last(axis)=first(axis)
        do k=first(3),last(3)
          do j=first(2),last(2)
            do i=first(1),last(1)
              idx=[i,j,k];inside=max(1,min(6,idx))
              select case(axis)
              case(1)
                area=patch%x_face_fraction(6*side,j,k)
              case(2)
                area=patch%y_face_fraction(i,6*side,k)
              case(3)
                area=patch%z_face_fraction(i,j,6*side)
              end select
              if(area<=0._dp) then
                call require(all(ghost(:,i,j,k)==0._dp).and.gt(i,j,k)==0._dp,'closed face unused')
                cycle
              end if
              donor=lo+(inside-1)/2
              donor(axis)=lo(axis)-1
              if(side==1) donor(axis)=hi(axis)+1
              external=donor(axis)<1.or.donor(axis)>6
              if(external) then
                call require(all(ghost(:,i,j,k)==u(:,inside(1),inside(2),inside(3))),'physical outflow donor')
                expected_t=t(inside(1),inside(2),inside(3))
              else
                call require(all(ghost(:,i,j,k)==parent_u(:,donor(1),donor(2),donor(3))),'interior parent donor')
                expected_t=parent_t(donor(1),donor(2),donor(3))
              end if
              call require(abs(gt(i,j,k)-expected_t)<=1.e-8_dp,'donor EOS temperature')
            end do
          end do
        end do
      end do
    end do
    call require(all(ghost(:,0,0,:)==0._dp),'unused corners zero')
    ! This fixture deliberately constructs face-only halos.
    call advance_reactive_eb_euler_3d(species,u,t,patch,'rusanov',1.e-10_dp,next,nt,ok,ghost,gt,flux,difmag=0._dp)
    call require(ok.and.maxval(abs(next-u))>1.e-10_dp,'patch update uses actual coarse donors')
    saved=ghost;saved_t=gt;saved_geometry=patch
    call fill_reactive_eb_patch_halo_3d(species,coarse,patch,lo+1,2,parent_u,parent_t,u,t,ghost,gt,ok)
    call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'misaligned patch halo rollback')
    call extract_eb_geometry_3d(full,[0,1,1],[6,6,6],patch,ok)
    call require(.not.ok.and.all(patch%volume_fraction==saved_geometry%volume_fraction),'invalid crop rollback')
  end do
  parent_u(1,4,3,3)=-1._dp
  call fill_reactive_eb_patch_halo_3d(species,coarse,patch,lo,2,parent_u,parent_t,u,t,ghost,gt,ok)
  call require(.not.ok.and.all(ghost==saved).and.all(gt==saved_t),'bad parent donor halo rollback')
  call check_peers()
  write(*,'(a)') 'test_eb_patch_halo_3d: PASS'
contains
  subroutine check_peers()
    type(eb_geometry_3d) :: target,peer,bad_peer
    real(dp), allocatable :: pu(:,:,:,:),pt(:,:,:),gh(:,:,:,:),ght(:,:,:),keep(:,:,:,:),keep_t(:,:,:)
    integer :: a,s,tl(3),th(3),pl(3),ph(3),offset(3),pn(3),ii,jj,kk,q(3),source(3),bad_index(3),changed
    real(dp) :: aperture
    logical :: valid,replaced
    allocate(gh(nvar,0:7,0:7,0:7),ght(0:7,0:7,0:7))
    do a=1,3
      do s=0,1
        tl=1;th=6;pl=2;ph=5
        if (s==1) then
          pl(a)=7;ph(a)=12
        else
          tl(a)=7;th(a)=12;pl(a)=1;ph(a)=6
        end if
        offset=pl-tl;pn=ph-pl+1
        call extract_eb_geometry_3d(full,tl,th,target,valid)
        call require(valid,'peer target extraction')
        call extract_eb_geometry_3d(full,pl,ph,peer,valid)
        call require(valid,'partial-face peer extraction')
        allocate(pu(nvar,pn(1),pn(2),pn(3)),pt(pn(1),pn(2),pn(3)))
        config%nx=pn(1);config%ny=pn(2);config%nz=pn(3)
        call initialize_reactive_eb_density_sheet_3d(species,config,peer,pu,pt,valid)
        call require(valid,'peer state initialization')
        do kk=1,pn(3)
          do jj=1,pn(2)
            do ii=1,pn(1)
              pu(:,ii,jj,kk)=pu(:,ii,jj,kk)*(1._dp+.001_dp*(ii+10*jj+100*kk))
            end do
          end do
        end do
        gh=-77._dp;ght=-88._dp
        call overlay_reactive_eb_peer_halo_3d(species,target,peer,offset,pu,pt,gh,ght,valid)
        call require(valid,'same-level peer overlay')
        changed=0
        do kk=0,7
          do jj=0,7
            do ii=0,7
              q=[ii,jj,kk];source=q-offset
              replaced=q(a)==s*7.and.all(source>=1).and.all(source<=pn)
              aperture=0._dp
              if (replaced) then
                select case(a)
                case(1);aperture=target%x_face_fraction(s*6,jj,kk)
                case(2);aperture=target%y_face_fraction(ii,s*6,kk)
                case(3);aperture=target%z_face_fraction(ii,jj,s*6)
                end select
              end if
              if (replaced.and.aperture>0._dp) then
                changed=changed+1
                call require(all(gh(:,ii,jj,kk)==pu(:,source(1),source(2),source(3))),'exact peer conserved donor')
                call require(abs(ght(ii,jj,kk)-pt(source(1),source(2),source(3)))< &
                  1.e-8_dp*pt(source(1),source(2),source(3)),'peer EOS temperature')
              else
                call require(all(gh(:,ii,jj,kk)==-77._dp).and.ght(ii,jj,kk)==-88._dp,'other halo slots preserved')
              end if
            end do
          end do
        end do
        call require(changed>0,'peer face actually exercised')
        keep=gh;keep_t=ght
        bad_index=pn;bad_index(a)=pn(a)
        if (s==1) bad_index(a)=1
        pu(1,bad_index(1),bad_index(2),bad_index(3))=-1._dp
        call overlay_reactive_eb_peer_halo_3d(species,target,peer,offset,pu,pt,gh,ght,valid)
        call require(.not.valid.and.all(gh==keep).and.all(ght==keep_t),'invalid peer donor rolls back entire overlay')
        call overlay_reactive_eb_peer_halo_3d(species,target,peer,[0,0,0],pu,pt,gh,ght,valid)
        call require(.not.valid.and.all(gh==keep).and.all(ght==keep_t),'overlapping boxes rejected')
        if (a/=1) then
          call build_axis_plane_eb_geometry_3d(pn(1),pn(2),pn(3),peer%x_lower,peer%x_upper, &
            peer%y_lower,peer%y_upper,peer%z_lower,peer%z_upper,'x',.171_dp,bad_peer,valid)
          call require(valid,'individually valid mismatched interface geometry')
          call initialize_reactive_eb_density_sheet_3d(species,config,bad_peer,pu,pt,valid)
          call require(valid,'mismatched peer valid state')
          call overlay_reactive_eb_peer_halo_3d(species,target,bad_peer,offset,pu,pt,gh,ght,valid)
          call require(.not.valid.and.all(gh==keep).and.all(ght==keep_t),'shared aperture mismatch rejected')
        end if
        deallocate(pu,pt)
      end do
    end do
  end subroutine check_peers

  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) error stop message
  end subroutine require
end program test_eb_patch_halo_3d
