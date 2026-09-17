program test_eb_reactive_transfer_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d,prolong_reactive_eb_state_pcm_3d
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(eb_geometry_3d) :: fine,coarse
  type(reactive_eb_3d_config) :: config
  real(dp), allocatable :: uf(:,:,:,:),tf(:,:,:),uc(:,:,:,:),tc(:,:,:),saved(:,:,:,:),saved_t(:,:,:)
  real(dp), allocatable :: prolonged(:,:,:,:),prolonged_t(:,:,:),primitive(:)
  real(dp) :: phi(0:6,0:6,0:6),temperature,speed,nan_value,before,after,vf,vc
  integer :: mode,ratio,i,j,k,n,nvar,c
  logical :: ok
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo load')
  nvar=reactive_nvar(size(species))
  allocate(uf(nvar,6,6,6),tf(6,6,6),primitive(reactive_nprim(size(species))))
  allocate(prolonged(nvar,6,6,6),prolonged_t(6,6,6))
  nan_value=ieee_value(0._dp,ieee_quiet_nan)
  config%nx=6;config%ny=6;config%nz=6
  do mode=1,2
    do k=0,6
      do j=0,6
        do i=0,6
          if (mode==1) then
            phi(i,j,k)=real(i+j+k,dp)/6-1.231_dp
          else
            phi(i,j,k)=(real(i,dp)/6-.5_dp)**2+(real(j,dp)/6-.5_dp)**2+ &
              (real(k,dp)/6-.5_dp)**2-.29_dp**2
          end if
        end do
      end do
    end do
    call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[1._dp,1._dp,1._dp],fine,ok)
    call require(ok,'fine geometry')
    call initialize_reactive_eb_density_sheet_3d(species,config,fine,uf,tf,ok)
    call require(ok,'initial state')
    do k=1,6
      do j=1,6
        do i=1,6
          if (fine%volume_fraction(i,j,k)<=0._dp) then
            uf(:,i,j,k)=nan_value;tf(i,j,k)=nan_value
            cycle
          end if
          call reactive_conserved_to_primitive(species,uf(:,i,j,k),tf(i,j,k),primitive,temperature,speed,ok)
          call require(ok,'initial primitive')
          primitive(1)=primitive(1)*(1._dp+.02_dp*i)
          primitive(2:4)=[real(i,dp),real(j,dp),real(k,dp)]
          primitive(5)=primitive(5)*(1._dp+.01_dp*(j+k))
          call reactive_primitive_to_conserved(species,primitive,uf(:,i,j,k),tf(i,j,k),speed,ok)
          call require(ok,'variable initial state')
        end do
      end do
    end do
    do ratio=2,3
      n=6/ratio
      call coarsen_eb_geometry_3d(fine,ratio,coarse,ok)
      call require(ok,'coarse geometry')
      allocate(uc(nvar,n,n,n),tc(n,n,n),saved(nvar,n,n,n),saved_t(n,n,n))
      uc=nan_value;tc=nan_value
      call restrict_reactive_eb_state_3d(species,fine,ratio,uf,tf,uc,tc,ok)
      call require(ok,'volume weighted restriction')
      vf=fine%dx*fine%dy*fine%dz;vc=coarse%dx*coarse%dy*coarse%dz
      do c=1,nvar
        before=sum(uf(c,:,:,:)*fine%volume_fraction,mask=fine%volume_fraction>0)*vf
        after=sum(uc(c,:,:,:)*coarse%volume_fraction,mask=coarse%volume_fraction>0)*vc
        call require(abs(after-before)<=2.e-12_dp*max(1._dp,abs(before)),'restriction conserved component')
      end do
      do k=1,n
        do j=1,n
          do i=1,n
            if (coarse%volume_fraction(i,j,k)<=0._dp) then
              call require(all(same_bits(uc(:,i,j,k),nan_value)).and.same_bits(tc(i,j,k),nan_value), &
                'covered coarse storage retained')
              cycle
            end if
            call reactive_conserved_to_primitive(species,uc(:,i,j,k),tc(i,j,k),primitive,temperature,speed,ok)
            call require(ok.and.abs(temperature-tc(i,j,k))<=1.e-8_dp,'restricted EOS temperature')
          end do
        end do
      end do
      saved=uc;saved_t=tc
      prolonged=uf;prolonged_t=tf
      call prolong_reactive_eb_state_pcm_3d(species,fine,ratio,uc,tc,prolonged,prolonged_t,ok)
      call require(ok,'PCM prolongation')
      do c=1,nvar
        before=sum(uc(c,:,:,:)*coarse%volume_fraction,mask=coarse%volume_fraction>0)*vc
        after=sum(prolonged(c,:,:,:)*fine%volume_fraction,mask=fine%volume_fraction>0)*vf
        call require(abs(after-before)<=2.e-12_dp*max(1._dp,abs(before)),'prolongation conserved component')
      end do
      call require(all(same_bits(prolonged,uf).or.spread(fine%volume_fraction>0,1,nvar)).and. &
        all(same_bits(prolonged_t,tf).or.fine%volume_fraction>0),'covered fine storage retained')
      call restrict_reactive_eb_state_3d(species,fine,ratio,prolonged,prolonged_t,uc,tc,ok)
      call require(ok,'PCM restriction roundtrip')
      do c=1,nvar
        call require(maxval(abs(uc(c,:,:,:)-saved(c,:,:,:)),mask=coarse%volume_fraction>0)<= &
          2.e-12_dp*max(1._dp,maxval(abs(saved(c,:,:,:)),mask=coarse%volume_fraction>0)), &
          'roundtrip conserved state')
      end do
      call require(maxval(abs(tc-saved_t),mask=coarse%volume_fraction>0)<=1.e-8_dp,'roundtrip temperature')
      saved=uc;saved_t=tc
      prolonged(1,6,6,6)=-1._dp
      call restrict_reactive_eb_state_3d(species,fine,ratio,prolonged,prolonged_t,uc,tc,ok)
      call require(.not.ok.and.all(same_bits(uc,saved)).and.all(same_bits(tc,saved_t)), &
        'bad late donor restriction rollback')
      uc(1,n,n,n)=-1._dp
      prolonged=uf;prolonged_t=tf
      call prolong_reactive_eb_state_pcm_3d(species,fine,ratio,uc,tc,prolonged,prolonged_t,ok)
      call require(.not.ok.and.all(same_bits(prolonged,uf)).and.all(same_bits(prolonged_t,tf)), &
        'bad late donor prolongation rollback')
      uc=saved
      call restrict_reactive_eb_state_3d(species,fine,0,uf,tf,uc,tc,ok)
      call require(.not.ok.and.all(same_bits(uc,saved)).and.all(same_bits(tc,saved_t)),'bad ratio rollback')
      call prolong_reactive_eb_state_pcm_3d(species,fine,ratio,uc(:,:n-1,:,:),tc,prolonged,prolonged_t,ok)
      call require(.not.ok.and.all(same_bits(prolonged,uf)).and.all(same_bits(prolonged_t,tf)), &
        'bad shape rollback')
      deallocate(uc,tc,saved,saved_t)
    end do
  end do
  write(*,'(a)') 'test_eb_reactive_transfer_3d: PASS'
contains
  pure elemental logical function same_bits(a,b)
    real(dp), intent(in) :: a,b
    same_bits=transfer(a,0_int64)==transfer(b,0_int64)
  end function same_bits
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require
end program test_eb_reactive_transfer_3d
