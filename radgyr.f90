subroutine radgyr_tensor(rg,r,ID,mC,tot_atom,N,box,binv)
integer, intent(in) :: N,tot_atom
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
integer, intent(in), dimension(0:N-1) :: ID
real(4), intent(in), dimension(0:N-1) :: mC
real(8), intent(in), dimension(0:2) :: box,binv
real(4), intent(inout), dimension(0:2) :: rg
real(4), dimension(0:2) :: rCOM,dr
integer :: i,j,k

dr(:) = 0.
rg(:) = 0.
rCOM(:) = 0.

do i=0, N-1
!call pbcdr(dr,r(ID(i)-1,0:2),(/0.,0.,0./),box,binv)
call pbcdr(dr,r(ID(i)-1,0:2),r(ID(0)-1,0:2),box,binv)
rCOM(0:2) = rCOM(0:2) + dr(0:2)*mC(i)
end do
rCOM(0:2) = rCOM(0:2)/SUM(mC) + r(ID(0)-1,0:2)

rg(:)=0.0
do i=0, N-1
  call pbcdr(dr,r(ID(i)-1,0:2),rCOM(0:2),box,binv)
  !do j=0,2
  rg(:) = rg(:) + dr(:)*dr(:)
  !end do
end do
rg(:) = rg(:)/N
end subroutine radgyr_tensor
