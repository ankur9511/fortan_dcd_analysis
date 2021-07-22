subroutine pbcdr(dr,r1,r2,box,binv)
real(8), intent(in), dimension(0:2) :: box,binv
real(4), dimension(0:2), intent(in) :: r1,r2
real(4), dimension(0:2), intent(inout) :: dr

dr = r1 - r2
dr = dr - box*anint(dr*binv)

end subroutine pbcdr