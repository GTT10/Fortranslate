module reactive_csv_io_3d_mod
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use state_indices_mod, only: imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive
  implicit none
  private

  public :: write_reactive_3d_csv,reactive_3d_csv_header,format_reactive_3d_csv_row

contains

  function reactive_3d_csv_header(species) result(header)
    type(nasa7_species), intent(in) :: species(:)
    character(len=:), allocatable :: header
    integer :: i
    header='time,x,y,z,rho,u,v,w,pressure,temperature,rhoE,rhou,rhov,rhow'
    do i=1,size(species)
      header=header//',Y_'//trim(species(i)%name)
    end do
    do i=1,size(species)
      header=header//',rhoY_'//trim(species(i)%name)
    end do
  end function reactive_3d_csv_header

  subroutine format_reactive_3d_csv_row(species,xyz,time,state,temperature,row,ok,internal_energy_density)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: xyz(3),time,state(:),temperature
    character(len=*), intent(out) :: row
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: internal_energy_density
    real(dp) :: primitive(reactive_nprim(size(species))),local_temperature,sound_speed
    integer :: i,status
    ok=.false.;row=''
    if (size(state)/=reactive_nvar(size(species)).or.len(row)/=26*(14+2*size(species))-1) return
    call reactive_conserved_to_primitive(species,state,temperature,primitive,local_temperature, &
      sound_speed,ok,internal_energy_density=internal_energy_density)
    if (.not.ok) return
    write(row,'(*(es25.16e3,:,","))',iostat=status) time,xyz,primitive(1:5),local_temperature, &
      state(iet),state(imx),state(imy),state(imz), &
      (primitive(reactive_mass_fraction_component(i)),i=1,size(species)), &
      (state(reactive_species_component(i)),i=1,size(species))
    ok=status==0
  end subroutine format_reactive_3d_csv_row


  subroutine write_reactive_3d_csv( &
      path, species, x, y, z, state, temperature, nx, ny, nz, time, &
      ok, message,heat_release,rotating_potential,internal_energy_density)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: time
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), optional, intent(in) :: heat_release(:,:,:),rotating_potential(:,:,:),internal_energy_density(:,:,:)
    real(dp), allocatable :: auxiliary_cell

    real(dp) :: ordinary_cell(reactive_nvar(size(species)))
    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed
    logical :: cell_ok
    integer :: unit, io_status, i, j, k, species_index, nspecies

    ok = .false.
    message = ""
    nspecies = size(species)
    if (size(x) /= nx .or. size(y) /= ny .or. size(z) /= nz .or. &
        size(state, 1) /= reactive_nvar(nspecies) .or. &
        size(state, 2) /= nx .or. size(state, 3) /= ny .or. &
        size(state, 4) /= nz .or. size(temperature, 1) /= nx .or. &
        size(temperature, 2) /= ny .or. size(temperature, 3) /= nz) then
      message = "Reactive 3D CSV array shapes are inconsistent"
      return
    end if

    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=[nx,ny,nz])) then
        message='Reactive 3D CSV internal-energy shape is inconsistent';return
      endif
      if(.not.all(ieee_is_finite(internal_energy_density))) then
        message='Reactive 3D CSV internal energy is nonfinite';return
      endif
      allocate(auxiliary_cell)
    endif
    if(present(heat_release)) then
      if(any(shape(heat_release)/=[nx,ny,nz])) then
        message='Reactive 3D CSV heat shape is inconsistent';return
      endif
      if(.not.all(ieee_is_finite(heat_release))) then
        message='Reactive 3D CSV heat contains nonfinite values';return
      endif
    endif

    if(present(rotating_potential)) then
      if(any(shape(rotating_potential)/=[nx,ny,nz])) then
        message='Reactive 3D CSV rotation shape is inconsistent';return
      endif
      if(.not.all(ieee_is_finite(rotating_potential)).or.any(rotating_potential<0)) then
        message='Reactive 3D CSV rotation is invalid';return
      endif
    endif
    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=io_status)
    if (io_status /= 0) then
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if

    write(unit, '(a)', advance='no') &
      "time,x,y,z,rho,u,v,w,pressure,temperature,rhoE,rhou,rhov,rhow"
    do species_index = 1, nspecies
      write(unit, '(a)', advance='no') &
        ",Y_" // trim(species(species_index)%name)
    end do
    do species_index = 1, nspecies
      write(unit, '(a)', advance='no') &
        ",rhoY_" // trim(species(species_index)%name)
    end do
    if(present(heat_release)) write(unit,'(a)',advance='no') ',heat_release_W_m3'
    write(unit, '(a)') ""

    allocate(primitive(reactive_nprim(nspecies)))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          ordinary_cell=state(:,i,j,k)
          if(present(rotating_potential)) ordinary_cell(iet)=ordinary_cell(iet)+ &
            ordinary_cell(1)*rotating_potential(i,j,k)
          if(allocated(auxiliary_cell)) auxiliary_cell=internal_energy_density(i,j,k)
          call reactive_conserved_to_primitive( &
            species, ordinary_cell, temperature(i, j, k), primitive, &
            local_temperature, sound_speed, cell_ok,internal_energy_density=auxiliary_cell)
          if (.not. cell_ok) then
            close(unit)
            write(message, '(a,i0,a,i0,a,i0,a)') &
              "Non-physical reactive state while writing cell (", i, ",", &
              j, ",", k, ")"
            return
          end if
          write(unit, '(*(es25.16e3,:,","))',advance='no') &
            time, x(i), y(j), z(k), primitive(1), primitive(2), &
            primitive(3), primitive(4), primitive(5), local_temperature, &
            state(iet, i, j, k), state(imx, i, j, k), &
            state(imy, i, j, k), state(imz, i, j, k), &
            (primitive(reactive_mass_fraction_component(species_index)), &
              species_index = 1, nspecies), &
            (state(reactive_species_component(species_index), i, j, k), &
              species_index = 1, nspecies)
          if(present(heat_release)) then
            write(unit,'(a,es25.16e3)') ',',heat_release(i,j,k)
          else
            write(unit,'(a)') ''
          endif
        end do
      end do
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_3d_csv

end module reactive_csv_io_3d_mod
