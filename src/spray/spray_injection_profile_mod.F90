module spray_injection_profile_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: iostat_end,iostat_eor
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: spray_profile_max_points=16384
  type, public :: spray_rate_profile
    real(dp), allocatable :: times(:),rates(:)
  end type
  integer, parameter, public :: spray_schedule_max_events=65536
  type, public :: spray_event_schedule
    real(dp), allocatable :: times(:),ends(:),scales(:)
    integer, allocatable :: ports(:)
  end type
  public :: read_spray_event_schedule
  public :: integrate_spray_injection_profile,read_spray_injection_profile,read_spray_rate_profile
contains
  ! Four SI columns: begin, end, integer port, scale. Publish only a complete table.
  subroutine read_spray_event_schedule(path,schedule,ok)
    character(len=*), intent(in) :: path
    type(spray_event_schedule), intent(inout) :: schedule
    logical, intent(out) :: ok
    type(spray_event_schedule) :: candidate
    character(len=4096) :: line
    integer :: unit,ios,n,i,words,cut
    logical :: in_word
    ok=.false.;n=0
    allocate(candidate%times(spray_schedule_max_events),candidate%ends(spray_schedule_max_events), &
      candidate%scales(spray_schedule_max_events),candidate%ports(spray_schedule_max_events))
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if(ios/=0) return
    do
      read(unit,'(a)',advance='no',iostat=ios) line
      if(ios==iostat_end) exit
      if(ios/=iostat_eor) go to 910
      cut=scan(line,'#!')
      if(cut>0) line(cut:)=' '
      words=0;in_word=.false.
      do i=1,len_trim(line)
        if(line(i:i)==' '.or.line(i:i)==achar(9)) then
          in_word=.false.
        else
          if(scan('0123456789+-.eEdD',line(i:i))==0) go to 910
          if(.not.in_word) words=words+1
          in_word=.true.
        end if
      end do
      if(words==0) cycle
      if(words/=4.or.n==spray_schedule_max_events) go to 910
      n=n+1
      read(line,*,iostat=ios) candidate%times(n),candidate%ends(n),candidate%ports(n),candidate%scales(n)
      if(ios/=0) go to 910
      if(.not.all(ieee_is_finite([candidate%times(n),candidate%ends(n),candidate%scales(n)]))) go to 910
      if(candidate%times(n)<0.or.candidate%ends(n)<candidate%times(n).or.candidate%ends(n)>1e6_dp) go to 910
      if(candidate%scales(n)<0.or.candidate%scales(n)>1e6_dp) go to 910
      if(candidate%ports(n)<1.or.candidate%ports(n)>16) go to 910
      if(n>1) then
        if(candidate%times(n)<candidate%times(n-1)) go to 910
      end if
    end do
    close(unit,iostat=ios)
    if(ios/=0.or.n<1) return
    schedule%times=candidate%times(:n);schedule%ends=candidate%ends(:n)
    schedule%ports=candidate%ports(:n);schedule%scales=candidate%scales(:n)
    ok=.true.
    return
910 close(unit)
  end subroutine
  ! Two whitespace-separated SI columns; blank lines and #/! comments are allowed.
  ! Publish a zero-padded table only after parsing and validating every record.
  subroutine read_spray_injection_profile(path,times,rates,count,ok)
    character(len=*), intent(in) :: path
    real(dp), intent(inout) :: times(:),rates(:)
    integer, intent(inout) :: count
    logical, intent(out) :: ok
    real(dp) :: candidate_times(size(times)),candidate_rates(size(rates)),mass
    character(len=4096) :: line
    integer :: unit,ios,n,i,words,cut
    logical :: in_word,valid
    ok=.false.;n=0;candidate_times=0;candidate_rates=0
    if(size(times)/=size(rates).or.size(times)<2.or.size(times)>spray_profile_max_points) return
    open(newunit=unit,file=path,status='old',action='read',iostat=ios)
    if (ios/=0) return
    do
      read(unit,'(a)',advance='no',iostat=ios) line
      if (ios==iostat_end) exit
      if (ios/=iostat_eor) go to 900
      cut=scan(line,'#!')
      if (cut>0) line(cut:)=' '
      words=0;in_word=.false.
      do i=1,len_trim(line)
        if (line(i:i)==' '.or.line(i:i)==achar(9)) then
          in_word=.false.
        else
          if (scan('0123456789+-.eEdD',line(i:i))==0) go to 900
          if (.not.in_word) words=words+1
          in_word=.true.
        end if
      end do
      if (words==0) cycle
      if (words/=2.or.n==size(times)) go to 900
      n=n+1
      read(line,*,iostat=ios) candidate_times(n),candidate_rates(n)
      if (ios/=0) go to 900
    end do
    close(unit,iostat=ios)
    if (ios/=0.or.n<2) return
    call integrate_spray_injection_profile(candidate_times(:n),candidate_rates(:n), &
      candidate_times(1),candidate_times(n),mass,valid)
    if (.not.valid) return
    times=candidate_times;rates=candidate_rates;count=n;ok=.true.
    return
900 close(unit)
  end subroutine
  ! Allocatable file entry; retain the caller's previous profile on failure.
  subroutine read_spray_rate_profile(path,profile,ok)
    character(len=*), intent(in) :: path
    type(spray_rate_profile), intent(inout) :: profile
    logical, intent(out) :: ok
    real(dp), allocatable :: times(:),rates(:)
    integer :: count
    allocate(times(spray_profile_max_points),rates(spray_profile_max_points))
    count=0
    call read_spray_injection_profile(path,times,rates,count,ok)
    if(.not.ok) return
    profile%times=times(:count);profile%rates=rates(:count)
  end subroutine
  pure subroutine integrate_spray_injection_profile(times,rates,begin_time,end_time,mass,ok)
    real(dp), intent(in) :: times(:),rates(:),begin_time,end_time
    real(dp), intent(out) :: mass
    logical, intent(out) :: ok
    real(dp) :: lo,hi,flo,fhi,width
    integer :: knot
    mass=0;ok=.false.
    if (size(times)<2.or.size(times)>spray_profile_max_points.or.size(rates)/=size(times)) return
    if (.not.all(ieee_is_finite(times)).or..not.all(ieee_is_finite(rates))) return
    if (.not.all(ieee_is_finite([begin_time,end_time]))) return
    if (any(times<0).or.any(times>1e6_dp).or.any(rates<0).or.any(rates>1)) return
    if (any(times(2:)<=times(:size(times)-1))) return
    if (begin_time<times(1).or.end_time>times(size(times)).or.end_time<begin_time) return
    do knot=1,size(times)-1
      lo=max(begin_time,times(knot));hi=min(end_time,times(knot+1))
      if (hi<=lo) cycle
      width=times(knot+1)-times(knot)
      flo=(lo-times(knot))/width;fhi=(hi-times(knot))/width
      flo=(1-flo)*rates(knot)+flo*rates(knot+1)
      fhi=(1-fhi)*rates(knot)+fhi*rates(knot+1)
      mass=mass+.5_dp*(hi-lo)*(flo+fhi)
    end do
    ok=.true.
  end subroutine
end module
