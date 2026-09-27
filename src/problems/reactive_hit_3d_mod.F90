module reactive_hit_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64,iostat_end
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_conserved_to_primitive
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_boundaries_3d
  implicit none
  private
  public :: initialize_reactive_hit_3d
contains
  ! Optional pinned HIT input: x-fast records x,y,z,u,v,w. Velocities are
  ! normalized by uin_norm and scaled by Mach*sound_speed/sqrt(3).
  subroutine initialize_reactive_hit_3d(input_path,config,species,x,y,z,state,temperature,ok)
    character(len=*), intent(in) :: input_path
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: x(:),y(:),z(:),temperature(:,:,:)
    real(dp), intent(inout) :: state(:,:,:,:)
    logical, intent(out) :: ok
    character(len=1024) :: input_file,line
    logical :: binfmt,valid
    integer :: inres,unit,ios,i,j,k,c,a,b,d,lo(3),hi(3),donor(3),left,right,mid
    integer(int64) :: bytes
    real(dp) :: mach_t0,uin_norm,sound,temp,scale,length,position(3),fraction(3),weight,velocity(3),rho
    real(dp), allocatable :: data(:,:,:,:),knots(:),width(:),candidate(:,:,:,:),primitive(:)
    namelist /reactive_hit/ input_file,inres,binfmt,mach_t0,uin_norm
    ok=.false.;input_file='';inres=0;binfmt=.false.;mach_t0=.1_dp;uin_norm=1
    open(newunit=unit,file=input_path,status='old',action='read',iostat=ios)
    if(ios/=0) return
    read(unit,nml=reactive_hit,iostat=ios)
    close(unit)
    if(ios>0) return
    if(len_trim(input_file)==0) then
      ok=inres==0;return
    end if
    if(config%problem/='uniform_reactor') return
    if(any(resolved_reactive_boundaries_3d(config)/='periodic')) return
    if(inres<2) return
    if(real(inres,dp)**3>real(huge(1)/6,dp)) return
    if(.not.all(ieee_is_finite([mach_t0,uin_norm]))) return
    if(mach_t0<0.or.uin_norm<=0) return
    allocate(data(6,inres,inres,inres),knots(inres),width(inres),stat=ios)
    if(ios/=0) return
    if(binfmt) then
      inquire(file=trim(input_file),size=bytes,iostat=ios)
      if(ios/=0) return
      if(bytes/=6_int64*int(inres,int64)**3*int(storage_size(0._dp)/8,int64)) return
      open(newunit=unit,file=trim(input_file),status='old',action='read',access='stream',form='unformatted',iostat=ios)
      if(ios/=0) return
      read(unit,iostat=ios) data
      close(unit)
      if(ios/=0) return
    else
      open(newunit=unit,file=trim(input_file),status='old',action='read',iostat=ios)
      if(ios/=0) return
      read(unit,'(a)',iostat=ios) line
      if(ios==0) then
        if(index(trim(adjustl(line)),'#')/=1) ios=1
      end if
      do k=1,inres
      do j=1,inres
      do i=1,inres
        if(ios/=0) exit
        read(unit,*,iostat=ios) data(:,i,j,k)
      end do
      end do
      end do
      if(ios==0) then
        do
          read(unit,'(a)',iostat=ios) line
          if(ios/=0.or.len_trim(line)/=0) exit
        end do
        if(ios==iostat_end) then
          ios=0
        else
          ios=1
        end if
      end if
      close(unit)
      if(ios/=0) return
    end if
    if(.not.all(ieee_is_finite(data))) return
    knots=data(1,:,1,1);width(2:)=knots(2:)-knots(:inres-1);width(1)=width(2)
    if(any(width<=0).or.knots(1)<0) return
    length=knots(inres)+.5_dp*width(inres)
    if(.not.ieee_is_finite(length).or.length<=0) return
    do k=1,inres
    do j=1,inres
    do i=1,inres
      if(any(abs(data(:3,i,j,k)-[knots(i),knots(j),knots(k)])>64*epsilon(1._dp)*length)) return
    end do
    end do
    end do
    allocate(primitive(size(state,1)))
    call reactive_conserved_to_primitive(species,state(:,1,1,1),temperature(1,1,1),primitive,temp,sound,valid)
    if(.not.valid) return
    scale=mach_t0*sound/(sqrt(3._dp)*uin_norm)
    candidate=state
    do k=1,size(z)
    do j=1,size(y)
    do i=1,size(x)
      position=mod([x(i),y(j),z(k)],length)
      do c=1,3
        left=1;right=inres
        if(position(c)>=knots(inres)) then
          left=inres
        else if(position(c)>knots(1)) then
          do while(right-left>1)
            mid=(left+right)/2
            if(position(c)>=knots(mid)) then
              left=mid
            else
              right=mid
            end if
          end do
        end if
        lo(c)=left;hi(c)=modulo(left,inres)+1
        fraction(c)=(position(c)-knots(left))/width(left)
      end do
      velocity=0
      do d=0,1
      do b=0,1
      do a=0,1
        donor=merge(hi,lo,[a,b,d]==1)
        weight=product(merge(fraction,1-fraction,[a,b,d]==1))
        velocity=velocity+weight*data(4:6,donor(1),donor(2),donor(3))
      end do
      end do
      end do
      rho=state(1,i,j,k)
      velocity=state(2:4,i,j,k)/rho+scale*velocity
      candidate(2:4,i,j,k)=rho*velocity
      candidate(5,i,j,k)=state(5,i,j,k)+.5_dp*(rho*sum(velocity**2)-sum(state(2:4,i,j,k)**2)/rho)
    end do
    end do
    end do
    if(.not.all(ieee_is_finite(candidate))) return
    state=candidate;ok=.true.
    write(*,'(a,a)') 'HIT velocity input: ',trim(input_file)
    write(*,'(a,es24.16)') 'HIT velocity scale: ',scale
  end subroutine
end module
