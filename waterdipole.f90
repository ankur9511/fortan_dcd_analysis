subroutine waterdipole(costheta,rO,rH1,rH2,box,binv)
implicit none
real(4), intent(in), dimension(0:2) :: rO,rH1,rH2
real(4), dimension(0:2) :: rOH1,rOH2,rOHH
real(8), intent(in), dimension(0:2) :: box,binv
real(4), intent(inout) :: costheta

call pbcdr(rOH1,rO,rH1,box,binv)
call pbcdr(rOH2,rO,rH2,box,binv)
rOHH(:) = 0.5*(rOH1(:)+rOH2(:))
!!!! a_vec (dot) z = |a_vec_z| = |a_vec| costheta
costheta = (rOHH(2)*rOHH(2))/(rOHH(0)*rOHH(0)+rOHH(1)*rOHH(1)+rOHH(2)*rOHH(2))
costheta = SIGN(sqrt(costheta),rOHH(2))
end subroutine waterdipole
