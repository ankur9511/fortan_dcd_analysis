subroutine qfourcalc(qt,neighborID,center,centerID,ncenters,shell,nshell,r,tot_atom,box,binv)
    implicit none
    integer, intent(in) :: tot_atom,ncenters,nshell
    real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
    real(8), intent(in), dimension(0:2) :: box,binv
    integer, intent(in) :: center,centerID
    integer, intent(in), dimension(0:nshell-1) :: shell
    integer, intent(in), dimension(0:ncenters-1,0:3) :: neighborID
    real(4), dimension(0:2) :: dr1,dr2
    real(4) :: dr12,dr22,costheta2
    integer :: k,k1,l,l1,lID,l1ID
    real(4), intent(inout) :: qt
    
    qt = 0.0
    do k = 0,2
     do k1 = k+1,3
     l = neighborID(center,k)
     l1 = neighborID(center,k1)
     if (l<0 .OR. l1<0) then
     print *, "Four NNs not found"
     end if
     if (l>=0 .AND. l1>=0) then
      lID = shell(l)-1
      l1ID = shell(l1)-1
      call pbcdr(dr1,r(lID,:),r(centerID,:),box,binv)
      dr12  = dr1(0)*dr1(0)+dr1(1)*dr1(1)+dr1(2)*dr1(2)
      dr12 = sqrt(dr12)
      call pbcdr(dr2,r(l1ID,:),r(centerID,:),box,binv)
      dr22  = dr2(0)*dr2(0)+dr2(1)*dr2(1)+dr2(2)*dr2(2)
      dr22 = sqrt(dr22)
      ! a_vec (dot) b_vec = |a||b|costheta
      costheta2 = dr1(0)*dr2(0)+dr1(1)*dr2(1)+dr1(2)*dr2(2)
      costheta2 = costheta2/dr12/dr22
      costheta2 = costheta2+1.0/3.0
      qt = qt+(costheta2*costheta2)
     end if
     end do
    end do
    qt = 1.0-0.375*qt
end subroutine qfourcalc
