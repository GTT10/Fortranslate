program test_state_redist_tracker_3d
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,build_level_set_eb_geometry_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_box_3d,select_state_redist_tracker_3d, &
    reactive_eb_weighted_state_redistribute_3d,state_redist_register_3d,initialize_state_redist_register_3d
  implicit none
  type(eb_geometry_3d) :: base,g
  type(state_redist_box_3d) :: box
  type(state_redist_register_3d) :: registers
  real(dp) :: state(1,17,17,17),result(1,17,17,17)
  real(dp) :: phi(0:17,0:17,0:17)
  integer :: coarse_flags(17,17,17)
  integer :: neighbors(3,7),count,id,cell(3),image(3),periodic(3),input,output,status,i,j,k,order
  real(dp) :: normal(3),kappa,neighbor_kappa,target,before,after,volume
  logical :: good
  character(len=4096) :: input_path,output_path
  call build_axis_plane_eb_geometry_3d(17,17,17,-6._dp,11._dp,-6._dp,11._dp,-6._dp,11._dp, &
    'x',-7._dp,base,good)
  if(.not.good) error stop 'Tracker base geometry'
  box%valid_lower=7;box%valid_upper=11
  box%domain_lower=7_int64;box%domain_upper=11_int64
  call get_command_argument(1,input_path)
  if(len_trim(input_path)==0) then
    cell=9;normal=[1._dp,1._dp,1._dp]
    kappa=.1_dp;neighbor_kappa=1._dp;target=.5_dp
    call prepare()
    call select_state_redist_tracker_3d(g,cell,target,box,neighbors,count,good)
    if(.not.good.or.count/=7) error stop 'Tracker three-way symmetry'
    if(any(neighbors(:,1)/=[0,0,1]).or.any(neighbors(:,2)/=[1,0,0]).or. &
      any(neighbors(:,3)/=[1,0,1]).or.any(neighbors(:,5)/=[1,1,0]).or. &
      any(neighbors(:,6)/=[0,1,1]).or.any(neighbors(:,7)/=[1,1,1])) error stop 'Tracker native order'
    normal=[1._dp,0._dp,0._dp]
    call prepare()
    call select_state_redist_tracker_3d(g,cell,target,box,neighbors,count,good)
    if(.not.good.or.count/=1.or.any(neighbors(:,1)/=[1,0,0])) error stop 'Tracker axial merge'
    normal=[1._dp,.2_dp,.1_dp];neighbor_kappa=.1_dp
    call prepare()
    call select_state_redist_tracker_3d(g,cell,target,box,neighbors,count,good)
    if(.not.good.or.count/=7) error stop 'Tracker volume expansion'
    normal=0
    call prepare()
    call select_state_redist_tracker_3d(g,cell,target,box,neighbors,count,good)
    if(good) error stop 'Tracker degenerate normal rejection'
    do k=0,17
    do j=0,17
    do i=0,17
      phi(i,j,k)=real(i+j+k,dp)-27.8_dp
    enddo
    enddo
    enddo
    call build_level_set_eb_geometry_3d(phi,[-6._dp,-6._dp,-6._dp],[11._dp,11._dp,11._dp],g,good)
    if(.not.good) error stop 'Tracker boundary level set setup'
    if(.not.g%is_valid()) error stop 'Tracker boundary geometry validity'
    box%periodic=.true.
    state=2
    call reactive_eb_weighted_state_redistribute_3d(g,state,result,good,box=box)
    if(.not.good) error stop 'Tracker periodic corner control'
    box%periodic=.false.
    state=2;coarse_flags=0
    call initialize_state_redist_register_3d(1,[17,17,17],registers,good,coarse_flags=coarse_flags)
    if(.not.good) error stop 'Tracker boundary register setup'
    registers%coarse_delta=3
    call reactive_eb_weighted_state_redistribute_3d(g,state,result,good,registers=registers,box=box)
    if(.not.good) error stop 'Tracker physical corner constant state'
    do k=7,11
    do j=7,11
    do i=7,11
      if(g%volume_fraction(i,j,k)<=0) cycle
      if(abs(result(1,i,j,k)-2._dp)>3e-12_dp) error stop 'Tracker physical constant preservation'
    end do
    end do
    end do
    if(any(registers%coarse_delta/=3)) error stop 'Tracker physical exterior makes no register transfer'
    do k=1,17
    do j=1,17
    do i=1,17
      state(1,i,j,k)=2._dp+.01_dp*(i-7)+.02_dp*(j-7)-.03_dp*(k-7)
    end do
    end do
    end do
    before=sum(state(1,7:11,7:11,7:11)*g%volume_fraction(7:11,7:11,7:11))
    do order=0,4,2
      call reactive_eb_weighted_state_redistribute_3d(g,state,result,good,box=box,max_order=order)
      if(.not.good) error stop 'Tracker physical corner nonconstant state'
      after=sum(result(1,7:11,7:11,7:11)*g%volume_fraction(7:11,7:11,7:11))
      if(abs(after-before)>3e-12_dp*max(abs(before),abs(after))) error stop 'Tracker physical conservation'
    end do
    ! The outward normal ends at a physical x face. Four tangent donors
    ! contain less than target; the connected grow1 ring can supply it.
    ! Tangent directions are periodic and this plane/state are periodic there.
    call build_axis_plane_eb_geometry_3d(17,17,17,-6._dp,11._dp,-6._dp,11._dp,-6._dp,11._dp, &
      'x',4.92_dp,g,good)
    if(.not.good) error stop 'Tracker physical tangent fallback geometry'
    box%periodic=[.false.,.true.,.true.];state=2
    call reactive_eb_weighted_state_redistribute_3d(g,state,result,good,box=box)
    if(.not.good) error stop 'Tracker physical connected tangent receivers'
    do k=7,11
    do j=7,11
      if(abs(result(1,11,j,k)-2._dp)>3e-12_dp) error stop 'Tracker tangent fallback constant preservation'
    end do
    end do
    do k=1,17
    do j=1,17
    do i=1,17
      state(1,i,j,k)=2._dp+.1_dp*cos(2._dp*acos(-1._dp)*real(j-7,dp)/5._dp)+ &
        .2_dp*sin(2._dp*acos(-1._dp)*real(k-7,dp)/5._dp)
    end do
    end do
    end do
    before=sum(state(1,7:11,7:11,7:11)*g%volume_fraction(7:11,7:11,7:11))
    do order=0,4,2
      call reactive_eb_weighted_state_redistribute_3d(g,state,result,good,box=box,max_order=order)
      if(.not.good) error stop 'Tracker tangent fallback periodic nonconstant state'
      after=sum(result(1,7:11,7:11,7:11)*g%volume_fraction(7:11,7:11,7:11))
      if(abs(after-before)>3e-12_dp*max(abs(before),abs(after))) error stop 'Tracker tangent periodic conservation'
    end do
    box%periodic=.false.
    ! Only a tiny isolated physical corner remains; all positive native
    ! tracker donors are exterior. No target-sized physical group exists.
    do k=0,17
    do j=0,17
    do i=0,17
      phi(i,j,k)=real(i+j+k,dp)-32.8_dp
    end do
    end do
    end do
    call build_level_set_eb_geometry_3d(phi,[-6._dp,-6._dp,-6._dp],[11._dp,11._dp,11._dp],g,good)
    if(.not.good.or..not.g%is_valid()) error stop 'Tracker isolated physical corner geometry'
    volume=sum(g%volume_fraction(7:11,7:11,7:11))
    if(volume<=0.or.volume>=.5_dp) error stop 'Tracker isolated physical volume control'
    state=2
    call reactive_eb_weighted_state_redistribute_3d(g,state,result,good,registers=registers,box=box)
    if(good) error stop 'Tracker insufficient physical volume accepted'
    if(any(registers%coarse_delta/=3)) error stop 'Tracker boundary register rollback'
    print *, 'PASS: native tracker symmetry, physical corner conservation orders0/2/4, isolated-volume rollback'
    stop
  endif
  call get_command_argument(2,output_path)
  open(newunit=input,file=trim(input_path),status='old',action='read')
  open(newunit=output,file=trim(output_path),status='replace',action='write')
  do
    read(input,*,iostat=status) id,image,normal,kappa,neighbor_kappa,target,periodic
    if(status<0) exit
    if(status/=0) error stop 'Tracker fixture read'
    cell=image+7;box%periodic=periodic/=0
    call prepare()
    call select_state_redist_tracker_3d(g,cell,target,box,neighbors,count,good)
    if(.not.good) then
      print *, 'Tracker fixture rejected: ',id
      error stop 'Tracker fixture'
    endif
    write(output,'(23(i0,1x))') id,count,neighbors
  enddo
  close(input);close(output)
contains
  subroutine prepare()
    g=base
    do k=cell(3)-1,cell(3)+1
    do j=cell(2)-1,cell(2)+1
    do i=cell(1)-1,cell(1)+1
      g%volume_fraction(i,j,k)=neighbor_kappa
    enddo
    enddo
    enddo
    g%volume_fraction(cell(1),cell(2),cell(3))=kappa
    g%x_face_fraction(cell(1)-1,cell(2),cell(3))=.5_dp-.25_dp*normal(1)
    g%x_face_fraction(cell(1),cell(2),cell(3))=.5_dp+.25_dp*normal(1)
    g%y_face_fraction(cell(1),cell(2)-1,cell(3))=.5_dp-.25_dp*normal(2)
    g%y_face_fraction(cell(1),cell(2),cell(3))=.5_dp+.25_dp*normal(2)
    g%z_face_fraction(cell(1),cell(2),cell(3)-1)=.5_dp-.25_dp*normal(3)
    g%z_face_fraction(cell(1),cell(2),cell(3))=.5_dp+.25_dp*normal(3)
  end subroutine prepare
end program test_state_redist_tracker_3d
