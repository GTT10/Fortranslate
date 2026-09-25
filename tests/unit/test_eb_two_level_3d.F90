program test_eb_two_level_3d
  use eb_reactive_transfer_3d_mod, only: fill_reactive_eb_patch_halo_3d
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d, &
    coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_two_level_3d_mod, only: advance_reactive_eb_two_level_euler_3d, &
    advance_reactive_eb_two_level_subcycled_euler_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: full,coarse,fine,bad
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: uc(:,:,:,:),tc(:,:,:),uf(:,:,:,:),tf(:,:,:),sc(:,:,:,:),stc(:,:,:),sf(:,:,:,:),stf(:,:,:)
  real(dp), allocatable :: boundary(:),wall(:),before(:),after(:),saved_b(:),saved_w(:)
  procedure(advance_reactive_eb_two_level_euler_3d), pointer :: advance
  logical :: refined(4,4,4),ok
  integer :: scheme,ratio,mode,n,nvar,lo(3),hi(3),step,c,i,j,k,fl(3),fh(3)
  character(len=32) :: argument
  real(dp) :: coefficient=0
  real(dp) :: coarse_amount,fine_amount,interval
  call get_command_argument(1,argument)
  if(argument=="viscosity") coefficient=.037_dp
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  nvar=reactive_nvar(size(species))
  allocate(uc(nvar,4,4,4),tc(4,4,4),sc(nvar,4,4,4),stc(4,4,4))
  allocate(boundary(nvar),wall(nvar),before(nvar),after(nvar),saved_b(nvar),saved_w(nvar))
  do ratio=2,3
    n=2*ratio
    allocate(uf(nvar,n,n,n),tf(n,n,n),sf(nvar,n,n,n),stf(n,n,n))
    call build_axis_plane_eb_geometry_3d(4*ratio,4*ratio,4*ratio,0._dp,1._dp,0._dp,1._dp, &
      0._dp,1._dp,'z',.31_dp,full,ok)
    call require(ok,'full fine geometry')
    call coarsen_eb_geometry_3d(full,ratio,coarse,ok)
    call require(ok,'coarse geometry')
    do scheme=1,2
    advance=>advance_reactive_eb_two_level_euler_3d
    interval=1.e-10_dp
    if (scheme==2) then
      advance=>advance_reactive_eb_two_level_subcycled_euler_3d
      interval=1.e-7_dp
    end if
    do mode=1,2
      lo=2
      if(mode==2) lo=1
      hi=lo+1
      refined=.false.;refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      call extract_eb_geometry_3d(full,(lo-1)*ratio+1,hi*ratio,fine,ok)
      call require(ok,'fine patch')
      config%nx=4;config%ny=4;config%nz=4
      call initialize_reactive_eb_density_sheet_3d(species,config,coarse,uc,tc,ok)
      call require(ok,'coarse state')
      config%nx=n;config%ny=n;config%nz=n
      call initialize_reactive_eb_density_sheet_3d(species,config,fine,uf,tf,ok)
      call require(ok,'fine state')
      if(coefficient>0) call check_complete_halo()
      boundary=17._dp;wall=19._dp
      do step=1,3
        call inventory(before)
        sc=uc;stc=tc;sf=uf;stf=tf
        call advance(species,coarse,fine,lo,ratio,'rusanov',interval, &
          uc,tc,uf,tf,boundary,wall,ok,difmag=coefficient)
        call require(ok,'synchronized two-level hydro step')
        call inventory(after)
        do c=1,nvar
          call require(abs(after(c)-before(c)-boundary(c)-wall(c))<=2.e-12_dp*max(1._dp,abs(before(c))), &
            'composite state equals boundary plus wall ledger')
        end do
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            do i=lo(1),hi(1)
              fl=([i,j,k]-lo)*ratio+1;fh=fl+ratio-1
              do c=1,nvar
                coarse_amount=uc(c,i,j,k)*coarse%volume_fraction(i,j,k)*coarse%dx*coarse%dy*coarse%dz
                fine_amount=sum(uf(c,fl(1):fh(1),fl(2):fh(2),fl(3):fh(3))* &
                  fine%volume_fraction(fl(1):fh(1),fl(2):fh(2),fl(3):fh(3)))*fine%dx*fine%dy*fine%dz
                call require(abs(coarse_amount-fine_amount)<=2.e-12_dp*max(1._dp,abs(fine_amount)), &
                  'refined parent synchronized with child inventory')
              end do
            end do
          end do
        end do
        call require(all(tc>0._dp).and.all(tf>0._dp),'positive temperatures')
        call require(maxval(abs(uf-sf))>1.e-12_dp,'fine evolution nonzero')
        if (scheme==2.and.step==1) then
          call advance_reactive_eb_two_level_euler_3d(species,coarse,fine,lo,ratio,'rusanov',interval, &
            sc,stc,sf,stf,saved_b,saved_w,ok,difmag=coefficient)
          call require(ok,'common-step reference for subcycle comparison')
          call require(maxval(abs(uf-sf))>1.e-12_dp,'subcycled fine evolution differs from one fine Euler step')
        end if
        call require(maxval(abs(wall))>0._dp,'embedded wall ledger nonzero')
        call require(all(uc(:,:,:,1)==sc(:,:,:,1)).and.all(tc(:,:,1)==stc(:,:,1)), &
          'covered coarse storage unchanged')
      end do
      sc=uc;stc=tc;sf=uf;stf=tf;saved_b=boundary;saved_w=wall
      call advance(species,coarse,fine,lo,ratio,'rusanov',1.e8_dp, &
        uc,tc,uf,tf,boundary,wall,ok,difmag=coefficient)
      call require(.not.ok,'nonphysical large step rejected')
      call unchanged()
      call build_axis_plane_eb_geometry_3d(n,n,n,fine%x_lower,fine%x_upper,fine%y_lower,fine%y_upper, &
        fine%z_lower,fine%z_upper,'z',.311_dp,bad,ok)
      call require(ok,'valid but inconsistent fine geometry')
      call advance(species,coarse,bad,lo,ratio,'rusanov',1.e-10_dp, &
        uc,tc,uf,tf,boundary,wall,ok,difmag=coefficient)
      call require(.not.ok,'coarse fine geometry mismatch rejected')
      call unchanged()
      call advance(species,coarse,fine,lo,0,'rusanov',1.e-10_dp, &
        uc,tc,uf,tf,boundary,wall,ok,difmag=coefficient)
      call require(.not.ok,'invalid ratio rejected')
      call unchanged()
    end do
    end do
    deallocate(uf,tf,sf,stf)
  end do
  write(*,'(a)') 'test_eb_two_level_3d: PASS'
contains
  subroutine check_complete_halo()
    real(dp) :: ghost(nvar,0:n+1,0:n+1,0:n+1),gt(0:n+1,0:n+1,0:n+1),vol(0:n+1,0:n+1,0:n+1)
    integer :: x,y,z,idx(3),parent(3)
    ghost=-17;gt=-18;vol=-19
    call fill_reactive_eb_patch_halo_3d(species,coarse,fine,lo,ratio,uc,tc,uf,tf,ghost,gt,ok,vol)
    call require(ok,'complete halo')
    do z=0,n+1
      do y=0,n+1
        do x=0,n+1
          idx=[x,y,z]
          where(lo==1) idx=max(1,idx)
          where(hi==4) idx=min(n,idx)
          if(all(idx>=1.and.idx<=n)) then
            call require(vol(x,y,z)==fine%volume_fraction(idx(1),idx(2),idx(3)),'fine halo metric')
            if(vol(x,y,z)>0) call require(all(ghost(:,x,y,z)==uf(:,idx(1),idx(2),idx(3))),'fine donor')
          else
            parent=lo+floor(real(idx-1,dp)/ratio)
            call require(vol(x,y,z)==coarse%volume_fraction(parent(1),parent(2),parent(3)),'parent halo metric')
            if(vol(x,y,z)>0) &
              call require(all(ghost(:,x,y,z)==uc(:,parent(1),parent(2),parent(3))),'parent donor')
          end if
          if(vol(x,y,z)>0) call require(gt(x,y,z)>0,'fluid halo temperature')
        end do
      end do
    end do
    sf=uf;sf(1,n,n,n)=-1
    ghost=-17;gt=-18;vol=-19
    call fill_reactive_eb_patch_halo_3d(species,coarse,fine,lo,ratio,uc,tc,sf,tf,ghost,gt,ok,vol)
    call require(.not.ok.and.all(ghost==-17).and.all(gt==-18).and.all(vol==-19),'complete halo rollback')
  end subroutine
  subroutine inventory(values)
    real(dp), intent(out) :: values(:)
    integer :: component
    do component=1,size(values)
      values(component)=coarse%dx*coarse%dy*coarse%dz* &
        sum(uc(component,:,:,:)*coarse%volume_fraction,mask=.not.refined)+ &
        fine%dx*fine%dy*fine%dz*sum(uf(component,:,:,:)*fine%volume_fraction)
    end do
  end subroutine inventory
  subroutine unchanged()
    call require(all(uc==sc).and.all(tc==stc).and.all(uf==sf).and.all(tf==stf).and. &
      all(boundary==saved_b).and.all(wall==saved_w),'whole hierarchy and ledger rollback')
  end subroutine unchanged
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(.not.condition) error stop message
  end subroutine require
end program test_eb_two_level_3d
