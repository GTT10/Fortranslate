program test_eb_reactive_transfer_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use nasa7_thermo_mod, only: nasa7_species,universal_gas_constant
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive, &
    reactive_primitive_to_conserved
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d,coarsen_eb_geometry_3d, &
    build_axis_plane_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d,prolong_reactive_eb_state_pcm_3d, &
    recover_reactive_eb_cell_3d,plan_eb_linear_halo_3d,interpolate_eb_linear_halo_3d, &
    reconstruct_eb_pressure_state_3d,interpolate_eb_pressure_time_3d
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
  block
    real(dp) :: left(nvar),right(nvar),mixed(nvar),t0,t1,h0,h1,cs,f,pvalue,p0,p1
    integer :: step
    primitive=0;primitive(1)=.4_dp;primitive(5)=60000;primitive(size(primitive))=1
    call reactive_primitive_to_conserved(species,primitive,left,t0,cs,ok)
    call require(ok,'temporal pressure left endpoint')
    h0=left(5);p0=primitive(5)
    primitive(1)=.2_dp;primitive(5)=66000
    call reactive_primitive_to_conserved(species,primitive,right,t1,cs,ok)
    call require(ok,'temporal pressure right endpoint')
    h1=right(5);p1=primitive(5)
    do step=0,4
      f=real(step,dp)/4;pvalue=-777
      call interpolate_eb_pressure_time_3d(species,left,right,h0,h1,t0,t1,f,pvalue,ok)
      call require(ok.and.abs(pvalue-((1-f)*p0+f*p1))<1.e-7_dp,'pressure uses parent endpoint clock')
    enddo
    mixed=.5_dp*(left+right)
    call reactive_conserved_to_primitive(species,mixed,.5_dp*(t0+t1),primitive,temperature,cs,ok, &
      internal_energy_density=.5_dp*(h0+h1))
    call require(ok.and.abs(primitive(5)-.5_dp*(p0+p1))>1._dp,'nonlinear mixing distinguishes pressure interpolation')
    pvalue=-777
    call interpolate_eb_pressure_time_3d(species,left,right,h0,h1,t0,t1,nan_value,pvalue,ok)
    call require(.not.ok.and.pvalue==-777,'temporal invalid fraction rollback')
    call interpolate_eb_pressure_time_3d(species,left,right,h0,nan_value,t0,t1,.5_dp,pvalue,ok)
    call require(.not.ok.and.pvalue==-777,'temporal invalid endpoint rollback')
  end block
  block
    real(dp) :: state(nvar),original(nvar),held,old_held,t,old_t,delta,cs,target,gap
    integer :: refusal
    primitive=0;primitive(1)=1;primitive(2:4)=[10._dp,-2._dp,3._dp]
    primitive(5)=101325;primitive(size(primitive))=1
    call reactive_primitive_to_conserved(species,primitive,state,t,cs,ok)
    call require(ok,'pressure reconstruction initial EOS')
    held=state(5)-.5_dp*sum(state(2:4)**2)/state(1)
    ! A nonzero dual-energy gap must survive the auxiliary representation.
    state(5)=state(5)+17._dp;original=state;old_held=held;old_t=t
    delta=-777;target=1.1_dp*primitive(5);gap=state(5)-held
    call reconstruct_eb_pressure_state_3d(species,state,held,t,target,delta,ok)
    call require(ok,'pressure reconstruction admitted')
    call require(all(same_bits(state(:4),original(:4))).and.all(same_bits(state(6:),original(6:))), &
      'pressure reconstruction preserves material and momentum bits')
    call require(abs(state(5)-held-gap)<1.e-10_dp,'pressure reconstruction preserves energy gap')
    call require(abs(state(5)-original(5)-delta)<1.e-10_dp.and.abs(held-old_held-delta)<1.e-10_dp, &
      'pressure reconstruction returns both energy corrections')
    call reactive_conserved_to_primitive(species,state,t,primitive,temperature,cs,ok,internal_energy_density=held)
    call require(ok.and.abs(primitive(5)-target)<1.e-8_dp,'pressure reconstruction target readback')
    do refusal=1,6
      state=original;held=old_held;t=old_t;delta=-777;target=101325
      select case(refusal)
      case(1)
        target=-1
      case(2)
        target=nan_value
      case(3)
        target=1.e20_dp ! Outside the pinned thermo temperature interval.
      case(4)
        held=nan_value
      case(5)
        t=nan_value
      case(6)
        state(1)=0
      end select
      block
        real(dp) :: kept(nvar),kept_scalars(3)
        kept=state;kept_scalars=[held,t,delta]
        call reconstruct_eb_pressure_state_3d(species,state,held,t,target,delta,ok)
        call require(.not.ok,'pressure reconstruction refuses invalid state/target')
        call require(all(same_bits(state,kept)).and.all(same_bits([held,t,delta],kept_scalars)), &
          'pressure reconstruction refusal is transactional including correction')
      end block
    enddo
  end block
  block
    type(eb_geometry_3d) :: extended
    integer(int64) :: origins(3,2),domains(3,2)
    integer :: shapes(3,2),requests(4,6),signs(3,6)
    real(dp) :: offset(3),value(2),neighbors(2,6),kept(2)
    character(len=16) :: faces(6)
    logical :: linear
    type(eb_geometry_3d) :: corner
    origins(:,1)=0;origins(:,2)=4;domains(:,1)=4;domains(:,2)=8;shapes=4;faces='outflow'
    call build_axis_plane_eb_geometry_3d(8,8,8,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp, &
      'z',.2_dp,extended,ok)
    call require(ok,'linear halo geometry')
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[3,3,3],2,[4_int64,4_int64,4_int64], &
      extended,2,faces,requests,signs,offset,linear,ok)
    call require(ok.and.linear.and.all(offset==-.25_dp),'connected coarse stencil has fine-center offset')
    call require(all(requests(:,1)==[1,2,3,3]).and.all(requests(:,2)==[1,4,3,3]),'axial donors')
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[2,2,1],2,[2_int64,2_int64,0_int64], &
      extended,2,faces,requests,signs,offset,linear,ok)
    call require(ok.and..not.linear,'covered neighbor retains PCM without reading solid state')
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[2,2,1],2,[2_int64,2_int64,0_int64], &
      extended,2,faces,requests,signs,offset,linear,ok,connected_directions=.false.)
    call require(ok.and..not.linear,'explicit false preserves pinned cut stencil')
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[2,2,1],2,[2_int64,2_int64,0_int64], &
      extended,2,faces,requests,signs,offset,linear,ok,connected_directions=.true.)
    call require(ok.and.linear,'directional cut stencil admits tangent reconstruction')
    call require(all(requests(:,1)==[1,1,2,1]).and.all(requests(:,2)==[1,3,2,1]).and. &
      all(requests(:,3)==[1,2,1,1]).and.all(requests(:,4)==[1,2,3,1]),'connected tangent donors retained')
    call require(all(requests(:,5)==[1,2,2,1]).and.all(requests(:,6)==[1,2,2,1]).and.all(signs==1), &
      'blocked normal axis reads donor twice without solid or reflection')
    block
      real(dp) :: center(2),samples(2,6),child(2),mean(2),shift(3)
      integer :: ix,iy,iz
      center=[10._dp,20._dp];samples=spread(center,2,6)
      samples(:,1)=center-[2._dp,4._dp];samples(:,2)=center+[2._dp,4._dp]
      samples(:,3)=center-[4._dp,8._dp];samples(:,4)=center+[4._dp,8._dp]
      mean=0
      do iz=0,1
      do iy=0,1
      do ix=0,1
        shift=.5_dp*real([ix,iy,iz],dp)-.25_dp;child=center
        call interpolate_eb_linear_halo_3d(child,samples,shift,ok)
        call require(ok.and.all(child==center+shift(1)*[2._dp,4._dp]+shift(2)*[4._dp,8._dp]), &
          'tangent affine interpolation ignores blocked normal offset')
        mean=mean+child/8
      enddo
      enddo
      enddo
      call require(all(mean==center),'directional symmetric children preserve donor average')
    end block
    corner=extended
    ! Synthetic planar cut-cell topology with one tangential diagonal sealed,
    ! while the center's immediate x/y neighbors remain connected.
    corner%x_face_fraction(4,5,3)=0
    corner%y_face_fraction(5,4,3)=0
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[2,2,1],2,[2_int64,2_int64,0_int64], &
      corner,2,faces,requests,signs,offset,linear,ok,connected_directions=.true.)
    call require(ok.and..not.linear,'disconnected tangent diagonal retains PCM')
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[2,2,1],2,[2_int64,2_int64,0_int64], &
      extended,2,faces,requests,signs,offset,linear,ok,connected_directions=.true.)
    call require(ok.and.linear,'corner fallback does not mutate source geometry')
    block
      type(eb_geometry_3d) :: shaped
      real(dp) :: nodal(0:8,0:8,0:8),xyz(3)
      integer :: shape_case,ix,iy,iz,q,admitted,constant
      do shape_case=1,2
        do iz=0,8
        do iy=0,8
        do ix=0,8
          xyz=-.5_dp+.25_dp*real([ix,iy,iz],dp)
          if(shape_case==1) then
            nodal(ix,iy,iz)=sum((xyz-.5_dp)**2)-.45_dp**2
          else
            nodal(ix,iy,iz)=min(xyz(1)-.38_dp,xyz(3)-.42_dp)
          endif
        enddo
        enddo
        enddo
        call build_level_set_eb_geometry_3d(nodal,[-.5_dp,-.5_dp,-.5_dp],[1.5_dp,1.5_dp,1.5_dp],shaped,ok)
        call require(ok,'curved/corner stencil geometry builds')
        admitted=0;constant=0
        do iz=1,4
        do iy=1,4
        do ix=1,4
          call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[ix,iy,iz],2, &
            2_int64*int([ix,iy,iz]-1,int64),shaped,2,faces,requests,signs,offset,linear,ok, &
            connected_directions=.true.)
          call require(ok,'curved/corner stencil returns valid plan or PCM')
          if(linear) then
            admitted=admitted+1
            do q=1,6
              call require(requests(1,q)==1.and.all(requests(2:4,q)>=1).and.all(requests(2:4,q)<=4), &
                'curved/corner donor stays on available coarse level')
              call require(shaped%volume_fraction(requests(2,q)+2,requests(3,q)+2,requests(4,q)+2)>0, &
                'curved/corner admitted stencil never reads solid')
            enddo
          else
            constant=constant+1
          endif
        enddo
        enddo
        enddo
        write(*,'(a,3(i0,1x))') 'Directional geometry case/admitted/PCM: ',shape_case,admitted,constant
        call require(admitted>0.and.constant>0,'curved/corner test exercises both admission and fallback')
      enddo
    end block
    faces(1:2)='periodic'
    call plan_eb_linear_halo_3d([0,1],shapes,origins,domains,1,[1,3,3],2,[0_int64,4_int64,4_int64], &
      extended,2,faces,requests,signs,offset,linear,ok)
    call require(ok.and.linear.and.all(requests(:,1)==[1,4,3,3]),'periodic stencil crosses root seam')
    value=[10._dp,20._dp];neighbors=spread(value,2,6)
    neighbors(:,1)=value-[2._dp,4._dp];neighbors(:,2)=value+[2._dp,4._dp]
    call interpolate_eb_linear_halo_3d(value,neighbors,[.25_dp,0._dp,0._dp],ok)
    call require(ok.and.all(value==[10.5_dp,21._dp]),'affine interpolation at fine cell center')
    value=[10._dp,20._dp];neighbors(2,1)=20
    call interpolate_eb_linear_halo_3d(value,neighbors,[.25_dp,0._dp,0._dp],ok)
    call require(ok.and.all(value==[10._dp,20._dp]),'stored energy extremum limits all gas components')
    kept=value;neighbors(2,2)=nan_value
    call interpolate_eb_linear_halo_3d(value,neighbors,[.25_dp,0._dp,0._dp],ok)
    call require(.not.ok.and.all(value==kept),'nonfinite stencil rollback')
  end block
  block
    real(dp) :: center(8),neighbors(8,6),value(8),legacy(8),explicit(8),average(8),short(2)
    integer :: ix,iy,iz
    center=[2._dp,1._dp,0._dp,0._dp,10._dp,.4_dp,1.6_dp,10._dp]
    neighbors=spread(center,2,6)
    neighbors([1,6,7],1)=.9_dp*center([1,6,7])
    neighbors([1,6,7],2)=1.1_dp*center([1,6,7])
    neighbors(2,1)=0;neighbors(2,2)=1.1_dp
    neighbors(3,1)=-1.e-26_dp;neighbors(3,2)=-1.e-31_dp
    legacy=center;explicit=center
    call interpolate_eb_linear_halo_3d(legacy,neighbors,[.25_dp,0._dp,0._dp],ok)
    call require(ok.and.all(legacy==center),'default retains shared momentum-extremum constraint')
    call interpolate_eb_linear_halo_3d(explicit,neighbors,[.25_dp,0._dp,0._dp],ok,separate_momentum=.false.)
    call require(ok.and.all(same_bits(legacy,explicit)),'explicit legacy policy preserves default bits')
    value=center
    call interpolate_eb_linear_halo_3d(value,neighbors,[.25_dp,0._dp,0._dp],ok,separate_momentum=.true.)
    call require(ok.and.abs(value(1)-2.05_dp)<1.e-14_dp,'candidate retains resolved scalar gradient')
    call require(abs(value(2)-1.05_dp)<1.e-14_dp.and.value(2)<=neighbors(2,2), &
      'candidate separately limits momentum instead of removing its bound')
    call require(value(3)==0.and.abs(sum(value(6:7))-value(1))<1.e-14_dp, &
      'candidate preserves species closure and zero momentum extremum')
    average=0
    do iz=-1,1,2
    do iy=-1,1,2
    do ix=-1,1,2
      value=center
      call interpolate_eb_linear_halo_3d(value,neighbors,.25_dp*real([ix,iy,iz],dp),ok,separate_momentum=.true.)
      call require(ok.and.all(value(6:7)>=0),'candidate positive species at every child')
      call require(abs(sum(value(6:7))-value(1))<1.e-14_dp,'candidate child species closure')
      average=average+value/8
    enddo
    enddo
    enddo
    call require(maxval(abs(average-center))<1.e-14_dp,'candidate child average conserves every component')
    value=center;neighbors(2,2)=nan_value
    call interpolate_eb_linear_halo_3d(value,neighbors,[.25_dp,0._dp,0._dp],ok,separate_momentum=.true.)
    call require(.not.ok.and.all(same_bits(value,center)),'candidate nonfinite refusal preserves caller')
    short=1
    call interpolate_eb_linear_halo_3d(short,spread(short,2,6),[0._dp,0._dp,0._dp],ok,separate_momentum=.true.)
    call require(.not.ok.and.all(short==1),'candidate refuses layout without momentum slots')
  end block
  block
    real(dp) :: center(8),neighbors(8,6),value(8),average(8),old(8),explicit(8)
    integer :: ix,iy,iz
    center=[2._dp,20._dp,-4._dp,0._dp,100._dp,.5_dp,1.5_dp,900._dp]
    neighbors=spread(center,2,6)
    neighbors(1:7,1)=.9_dp*center(1:7);neighbors(1:7,2)=1.1_dp*center(1:7)
    ! Tiny temperature extremum suppresses the shared scalar density slope.
    neighbors(8,1)=900._dp-1e-10_dp;neighbors(8,2)=900._dp-2e-10_dp
    old=center;explicit=center
    call interpolate_eb_linear_halo_3d(old,neighbors,[.25_dp,0._dp,0._dp],ok,separate_momentum=.true.)
    call require(ok.and.abs(old(2)/old(1)-10._dp)>.1_dp,'old separate policy breaks constant velocity')
    call interpolate_eb_linear_halo_3d(explicit,neighbors,[.25_dp,0._dp,0._dp],ok, &
      separate_momentum=.true.,consistent_velocity=.false.)
    call require(ok.and.all(same_bits(old,explicit)),'disabled velocity policy preserves bits')
    average=0
    do iz=-1,1,2
    do iy=-1,1,2
    do ix=-1,1,2
      value=center
      call interpolate_eb_linear_halo_3d(value,neighbors,.25_dp*real([ix,iy,iz],dp),ok, &
        separate_momentum=.true.,consistent_velocity=.true.)
      call require(ok.and.maxval(abs(value(2:4)/value(1)-[10._dp,-2._dp,0._dp]))<1e-13_dp, &
        'density-coupled momentum preserves constant velocity')
      average=average+value/8
    enddo
    enddo
    enddo
    call require(maxval(abs(average-center))<1e-13_dp,'density-coupled child average conservation')
    value=center
    call interpolate_eb_linear_halo_3d(value,neighbors,[0._dp,0._dp,0._dp],ok,consistent_velocity=.true.)
    call require(.not.ok.and.all(same_bits(value,center)),'velocity policy requires separate layout')
    neighbors=spread(center,2,6)
    neighbors(2,1)=center(2)-2;neighbors(2,2)=center(2)+2;value=center
    call interpolate_eb_linear_halo_3d(value,neighbors,[.25_dp,0._dp,0._dp],ok, &
      separate_momentum=.true.,consistent_velocity=.true.)
    call require(ok.and.abs(value(2)-20.5_dp)<1e-13_dp,'velocity difference retains resolved slope')
    neighbors(1,1)=0;value=center
    call interpolate_eb_linear_halo_3d(value,neighbors,[0._dp,0._dp,0._dp],ok, &
      separate_momentum=.true.,consistent_velocity=.true.)
    call require(.not.ok.and.all(same_bits(value,center)),'velocity policy refuses zero density transactionally')
  end block
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
  call check_held_recovery()
  call check_held_transfer()
  write(*,'(a)') 'test_eb_reactive_transfer_3d: PASS'
contains
  subroutine check_held_recovery()
    type(nasa7_species) :: constant_cv(1)
    type(rotating_frame_3d) :: frame
    real(dp) :: gas(6),initial(6),q(6),stored,rt,cs,r,phi,xyz(3),expected_cs
    integer :: index(3)
    logical :: good
    constant_cv(1)%name='X';constant_cv(1)%molecular_weight=10
    constant_cv(1)%temperature_min=200;constant_cv(1)%temperature_mid=1000;constant_cv(1)%temperature_max=5000
    r=universal_gas_constant/constant_cv(1)%molecular_weight
    constant_cv(1)%low_coefficients=0;constant_cv(1)%low_coefficients(1)=1+1000._dp/r
    constant_cv(1)%high_coefficients=constant_cv(1)%low_coefficients
    gas=[2._dp,2.e10_dp,0._dp,0._dp,1.e20_dp,2._dp];initial=gas;stored=600000._dp;index=[6,6,6]
    call require(gas(5)-.5_dp*gas(2)**2/gas(1)==0,'held EOS fixture has exact kinetic cancellation')
    call recover_reactive_eb_cell_3d(constant_cv,fine,index,gas,300._dp,q,rt,cs,good)
    call require(.not.good,'ordinary EOS rejects total-minus-kinetic energy below temperature bounds')
    call recover_reactive_eb_cell_3d(constant_cv,fine,index,gas,300._dp,q,rt,cs,good,internal_energy_density=stored)
    expected_cs=sqrt((1+r/1000._dp)*r*300._dp)
    call require(good.and.abs(rt-300._dp)<1.e-10_dp.and.abs(q(5)-2*r*300._dp)<1.e-8_dp.and. &
      abs(cs-expected_cs)<1.e-10_dp.and.q(2)==1.e10_dp.and.all(gas==initial), &
      'EB held EOS recovers analytic temperature/pressure/sound without kinetic subtraction')
    frame%axis=2;frame%angular_speed=1.e10_dp;frame%axis_point=[10._dp,10._dp,10._dp]
    xyz=[fine%x_lower,fine%y_lower,fine%z_lower]+(real(index,dp)-.5_dp)*[fine%dx,fine%dy,fine%dz]
    call rotating_potential_3d(frame,xyz,phi,good)
    call require(good,'held EOS rotating potential')
    gas(5)=gas(5)-gas(1)*phi;initial=gas
    call recover_reactive_eb_cell_3d(constant_cv,fine,index,gas,300._dp,q,rt,cs,good,frame,stored)
    call require(good.and.abs(rt-300._dp)<1.e-10_dp.and.abs(q(5)-2*r*300._dp)<1.e-8_dp.and. &
      abs(cs-expected_cs)<1.e-10_dp.and.all(gas==initial), &
      'rotating EB held EOS reads ordinary UEINT directly despite unresolved modified total energy')
    stored=ieee_value(0._dp,ieee_quiet_nan)
    call recover_reactive_eb_cell_3d(constant_cv,fine,index,gas,300._dp,q,rt,cs,good,frame,stored)
    call require(.not.good.and.all(q==0).and.rt==0.and.cs==0.and.all(gas==initial), &
      'invalid held EOS refuses with cleared outputs and unchanged conserved state')
  end subroutine check_held_recovery
  subroutine check_held_transfer()
    type(nasa7_species) :: held_species(1)
    type(eb_geometry_3d) :: cg
    type(rotating_frame_3d), allocatable :: frame
    real(dp), allocatable :: gas(:,:,:,:),ft(:,:,:),fh(:,:,:),coarse_gas(:,:,:,:),ct(:,:,:),ch(:,:,:)
    real(dp), allocatable :: pg(:,:,:,:),pt(:,:,:),ph(:,:,:),sg(:,:,:,:),st(:,:,:),sh(:,:,:)
    real(dp) :: rho,energy,vfine,vcoarse,original_inventory,new_inventory,r,phi,xyz(3),formation
    integer :: i,j,k,c,mode,last_f(3),last_c(3)
    logical :: good
    held_species(1)%name='X';held_species(1)%molecular_weight=10
    held_species(1)%temperature_min=200;held_species(1)%temperature_mid=1000;held_species(1)%temperature_max=5000
    r=universal_gas_constant/held_species(1)%molecular_weight
    held_species(1)%low_coefficients=0;held_species(1)%low_coefficients(1)=1+1000._dp/r
    held_species(1)%high_coefficients=held_species(1)%low_coefficients
    call coarsen_eb_geometry_3d(fine,2,cg,good)
    call require(good,'held transfer restricted geometry')
    allocate(gas(6,6,6,6),ft(6,6,6),fh(6,6,6),coarse_gas(6,3,3,3),ct(3,3,3),ch(3,3,3))
    vfine=fine%dx*fine%dy*fine%dz;vcoarse=cg%dx*cg%dy*cg%dz
    do mode=1,3
      formation=0
      if(mode==3) then
        if(allocated(frame)) deallocate(frame)
        formation=2.e6_dp
        held_species(1)%low_coefficients(6)=-formation/r
        held_species(1)%high_coefficients=held_species(1)%low_coefficients
      endif
      if(mode==2) then
        allocate(frame);frame%axis=2;frame%angular_speed=1.e10_dp;frame%axis_point=[10._dp,10._dp,10._dp]
      endif
      gas=nan_value;ft=nan_value;fh=nan_value
      do k=1,6
      do j=1,6
      do i=1,6
        if(fine%volume_fraction(i,j,k)<=0) cycle
        rho=1+.01_dp*i;ft(i,j,k)=300+10*i+5*j+2*k;energy=.5_dp*rho*1.e20_dp
        if(allocated(frame)) then
          xyz=[fine%x_lower,fine%y_lower,fine%z_lower]+ &
            (real([i,j,k],dp)-.5_dp)*[fine%dx,fine%dy,fine%dz]
          call rotating_potential_3d(frame,xyz,phi,good)
          call require(good,'held transfer native rotating energy convention')
          energy=energy-rho*phi
        endif
        gas(:,i,j,k)=[rho,rho*1.e10_dp,0._dp,0._dp,energy,rho]
        fh(i,j,k)=rho*(1000._dp*ft(i,j,k)-formation);last_f=[i,j,k]
      enddo
      enddo
      enddo
      coarse_gas=nan_value;ct=nan_value;ch=nan_value
      call restrict_reactive_eb_state_3d(held_species,fine,2,gas,ft,coarse_gas,ct,good,frame,fh,ch)
      call require(good,'paired gas/T/UEINT restriction despite unresolved total thermal energy')
      do k=1,3
      do j=1,3
      do i=1,3
        if(cg%volume_fraction(i,j,k)<=0) then
          call require(same_bits(ch(i,j,k),nan_value),'covered coarse UEINT storage retained')
          cycle
        endif
        call require(abs(ct(i,j,k)-(ch(i,j,k)/coarse_gas(1,i,j,k)+formation)/1000._dp)<1.e-10_dp, &
          'restricted temperature follows independently conserved analytic UEINT')
        last_c=[i,j,k]
      enddo
      enddo
      enddo
      original_inventory=sum(fh*fine%volume_fraction,mask=fine%volume_fraction>0)*vfine
      new_inventory=sum(ch*cg%volume_fraction,mask=cg%volume_fraction>0)*vcoarse
      call require(abs(new_inventory-original_inventory)<3.e-14_dp*max(1._dp,abs(original_inventory)), &
        'independent kappa-volume held inventory survives restriction')
      do c=1,6
        original_inventory=sum(gas(c,:,:,:)*fine%volume_fraction,mask=fine%volume_fraction>0)*vfine
        new_inventory=sum(coarse_gas(c,:,:,:)*cg%volume_fraction,mask=cg%volume_fraction>0)*vcoarse
        call require(abs(new_inventory-original_inventory)<3.e-14_dp*max(1._dp,abs(original_inventory)), &
          'paired restriction preserves each independent gas inventory')
      enddo
      pg=gas;pt=ft;ph=fh
      call prolong_reactive_eb_state_pcm_3d(held_species,fine,2,coarse_gas,ct,pg,pt,good,frame,ch,ph)
      call require(good,'paired UEINT PCM prolongation in ordinary/rotating coordinates')
      call require(all(same_bits(ph,fh).or.fine%volume_fraction>0),'covered fine UEINT storage retained')
      do k=1,6
      do j=1,6
      do i=1,6
        if(fine%volume_fraction(i,j,k)<=0) cycle
        call require(abs(pt(i,j,k)-(ph(i,j,k)/pg(1,i,j,k)+formation)/1000._dp)<1.e-10_dp, &
          'prolonged temperature follows held energy directly at the fine cell center')
      enddo
      enddo
      enddo
      original_inventory=sum(ch*cg%volume_fraction,mask=cg%volume_fraction>0)*vcoarse
      new_inventory=sum(ph*fine%volume_fraction,mask=fine%volume_fraction>0)*vfine
      call require(abs(new_inventory-original_inventory)<3.e-14_dp*max(1._dp,abs(original_inventory)), &
        'independent kappa-volume held inventory survives prolongation')
      sg=coarse_gas;st=ct;sh=ch;fh(last_f(1),last_f(2),last_f(3))=nan_value
      call restrict_reactive_eb_state_3d(held_species,fine,2,gas,ft,coarse_gas,ct,good,frame,fh,ch)
      call require(.not.good.and.all(same_bits(coarse_gas,sg)).and.all(same_bits(ct,st)).and.all(same_bits(ch,sh)), &
        'late held donor failure rolls back gas/T/UEINT restriction together')
      call restrict_reactive_eb_state_3d(held_species,fine,2,gas,ft,coarse_gas,ct,good,frame,coarse_internal_energy=ch)
      call require(.not.good.and.all(same_bits(coarse_gas,sg)).and.all(same_bits(ct,st)).and.all(same_bits(ch,sh)), &
        'incomplete paired held restriction arguments reject atomically')
      sg=pg;st=pt;sh=ph;ch(last_c(1),last_c(2),last_c(3))=nan_value
      call prolong_reactive_eb_state_pcm_3d(held_species,fine,2,coarse_gas,ct,pg,pt,good,frame,ch,ph)
      call require(.not.good.and.all(same_bits(pg,sg)).and.all(same_bits(pt,st)).and.all(same_bits(ph,sh)), &
        'late held donor failure rolls back gas/T/UEINT prolongation together')
    enddo
  end subroutine check_held_transfer

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
