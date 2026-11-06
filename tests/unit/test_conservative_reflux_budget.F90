program test_conservative_reflux_budget
  use precision_mod, only: dp
  use conservative_reflux_budget_mod, only: spread_reflux_return
  implicit none
  real(dp) :: u(9,5),candidate(9,5),saved(9,5),vol(5),seed(5),returned(9),flux(9,5),saved_flux(9,5)
  real(dp) :: before(9),after(9),balance(9,5),reordered(9,5),other_flux(9,5),coarse(9)
  integer :: edges(2,5),other_edges(2,5),i,e
  logical :: ok
  vol=[1._dp,.125_dp,.5_dp,2._dp,.25_dp]
  seed=[1._dp,0._dp,0._dp,0._dp,0._dp]
  edges(:,1)=[1,2];edges(:,2)=[2,3];edges(:,3)=[3,4];edges(:,4)=[4,5];edges(:,5)=[2,5]
  u=0;u(6,:)=[0._dp,0._dp,.2_dp,.1_dp,.05_dp];u(7,:)=.2_dp;u(8,:)=1.78_dp
  u(1,:)=sum(u(6:8,:),dim=1);u(2,:)=.5_dp*u(1,:);u(3,:)=-.2_dp*u(1,:)
  u(5,:)=10000;u(9,:)=1000
  returned=[-.12_dp,.06_dp,-.24_dp,.01_dp,-100._dp,-.12_dp,0._dp,0._dp,-30._dp]
  candidate=-9;flux=77
  call spread_reflux_return(u,vol,edges,seed,returned,6,8,candidate,flux,ok)
  call require(ok,'wide donor budget accepts when the boundary-adjacent cells have no fuel')
  call require(all(candidate(6:8,:)>=0).and.all(candidate(1,:)>0),'strict nonnegative species and positive density')
  call require(all(candidate(:,1:2)==u(:,1:2)),'empty adjacent donors retain their state')
  call require(all(candidate(6,3:5)<u(6,3:5)),'withdrawal reaches multiple remote donors')
  call require(maxval(abs(candidate(1,:)-sum(candidate(6:8,:),dim=1)))<1e-15_dp,'density follows species transfer')
  call inventory(u,before);call inventory(candidate,after)
  call require(maxval(abs(after-before-returned))<2e-11_dp,'all conserved inventories match the returned vector')
  ! An independent finite-volume balance reconstructs every cell from the
  ! external interface integral and the published internal face integrals.
  call audit_flux(edges,flux,candidate)
  coarse=0;coarse(1)=1;coarse(5)=10000;coarse(6:8)=[0._dp,.2_dp,.8_dp]
  call require(maxval(abs((coarse-returned)+(after-before)-coarse))<2e-11_dp,'coarse/fine composite conservation')
  saved=candidate;saved_flux=flux
  do e=1,5
    other_edges(:,e)=edges(2:1:-1,6-e)
  enddo
  reordered=0;other_flux=0
  call spread_reflux_return(u,vol,other_edges,seed,returned,6,8,reordered,other_flux,ok)
  call require(ok.and.all(reordered==candidate),'donor states do not depend on edge order or orientation')
  call audit_flux(other_edges,other_flux,reordered)
  returned(6)=-1;returned(1)=-1
  call spread_reflux_return(u,vol,edges,seed,returned,6,8,candidate,flux,ok)
  call require(.not.ok.and.all(candidate==saved).and.all(flux==saved_flux),'insufficient inventory rolls back state and flux')
  returned=0
  other_edges=edges;other_edges(:,4)=[3,4];other_edges(:,5)=[2,3]
  call spread_reflux_return(u,vol,other_edges,seed,returned,6,8,candidate,flux,ok)
  call require(.not.ok.and.all(candidate==saved).and.all(flux==saved_flux),'disconnected stencil rejects atomically')
  returned(1)=-sum(u(1,:)*vol)
  call spread_reflux_return(u,vol,edges,seed,returned,6,8,candidate,flux,ok)
  call require(.not.ok.and.all(candidate==saved),'cannot empty all donor density')
  returned=0
  call spread_reflux_return(u,vol,edges,seed,returned,6,8,candidate,flux,ok)
  call require(ok.and.all(candidate==u).and.all(flux==0),'zero return is unchanged')
  print *, 'PASS: wide conservative donor budget, mixed volumes, face balance, permutation and rollback'
contains
  subroutine inventory(values,total)
    real(dp), intent(in) :: values(:,:)
    real(dp), intent(out) :: total(:)
    integer :: c
    do c=1,size(values,1)
      total(c)=sum(values(c,:)*vol)
    enddo
  end subroutine
  subroutine audit_flux(graph,integrals,values)
    integer, intent(in) :: graph(:,:)
    real(dp), intent(in) :: integrals(:,:),values(:,:)
    integer :: face,cell
    balance=0;balance(:,1)=returned
    do face=1,size(graph,2)
      balance(:,graph(1,face))=balance(:,graph(1,face))-integrals(:,face)
      balance(:,graph(2,face))=balance(:,graph(2,face))+integrals(:,face)
    enddo
    do cell=1,5
      call require(maxval(abs((values(:,cell)-u(:,cell))*vol(cell)-balance(:,cell)))<2e-12_dp, &
        'published flux divergence equals accepted cell increment')
    enddo
  end subroutine
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(condition) return
    print *, 'FAIL: ',message
    error stop 1
  end subroutine
end program test_conservative_reflux_budget
