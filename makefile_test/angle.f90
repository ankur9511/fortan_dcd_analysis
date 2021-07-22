subroutine cosangle(costheta,rO,rH1,rH2,box,binv)
    real(4), intent(in), dimension(0:2) :: rO,rH1,rH2
    real(4), dimension(0:2) :: rOH1,rOH2,rOHH
    real(8), intent(in), dimension(0:2) :: box,binv
    real(4), intent(inout) :: costheta
    
    call pbcdr(rOH1,rH1,rO,box,binv)
    call pbcdr(rOH2,rH2,rO,box,binv)
    costheta = dot_product(rOH1,rOH2)/sqrt(dot_product(rOH1,rOH1)*dot_product(rOH2,rOH2))
end subroutine cosangle
