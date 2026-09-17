program test_eb_composition_integrals_3d
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: reactive_eb_composition_integrals_3d
  implicit none
  type(nasa7_species) :: species(3)
  type(eb_geometry_3d) :: geometry
  real(dp) :: state(8,2,2,2),atoms(4,3),values(4),saved(4),nan
  logical :: ok
  call build_axis_plane_eb_geometry_3d(2,2,2,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp, &
    'z',.6_dp,geometry,ok)
  if (.not.ok) error stop 'geometry'
  ! Deliberately non-formula labels: composition comes only from the table.
  species(1)%name='fuel_alias';species(1)%molecular_weight=16._dp
  species(2)%name='product_alias';species(2)%molecular_weight=44._dp
  species(3)%name='inert_alias';species(3)%molecular_weight=40._dp
  atoms(:,1)=[0._dp,1._dp,4._dp,0._dp] ! Ar,C,H,O for methane
  atoms(:,2)=[0._dp,1._dp,0._dp,2._dp] ! carbon dioxide
  atoms(:,3)=[1._dp,0._dp,0._dp,0._dp] ! argon
  state=0._dp;state(1,:,:,:)=100._dp;state(5,:,:,:)=1.e6_dp
  state(6,:,:,:)=16._dp;state(7,:,:,:)=44._dp;state(8,:,:,:)=40._dp
  nan=ieee_value(0._dp,ieee_quiet_nan)
  state(:,:,:,1)=nan ! Solid half does not own inventory.
  values=-7._dp
  call reactive_eb_composition_integrals_3d(species,state,geometry,atoms,values,ok)
  if (.not.ok.or.any(abs(values-[.4_dp,.8_dp,1.6_dp,.8_dp])>1.e-14_dp)) error stop 'analytical cut volume'
  saved=values
  atoms(1,1)=nan
  call reactive_eb_composition_integrals_3d(species,state,geometry,atoms,values,ok)
  if (ok.or.any(values/=saved)) error stop 'nonfinite metadata rollback'
  atoms(1,1)=0._dp;state(6,1,1,2)=-1._dp
  call reactive_eb_composition_integrals_3d(species,state,geometry,atoms,values,ok)
  if (ok.or.any(values/=saved)) error stop 'negative species rollback'
  state(6,1,1,2)=16._dp;atoms(:,3)=0._dp
  call reactive_eb_composition_integrals_3d(species,state,geometry,atoms,values,ok)
  if (ok.or.any(values/=saved)) error stop 'missing composition rollback'
  write(*,'(a)') 'test_eb_composition_integrals_3d: PASS'
end program test_eb_composition_integrals_3d
