subroutine endtoend_orient(costheta,box,binv,headcoord, &
                          & backbone,allcoord,nbackbone, &
                          & nallcoord)
integer, intent(in) :: nbackbone, nallcoord
real(4), intent(in), dimension(0:nallcoord-1,0:2) :: allcoord
real(4), intent(in), dimension(0:nbackbone-1) :: backbone
real(4), intent(in), dimension(0:2) :: headcoord
real(8), intent(in),dimension(0:2) :: box,binv
real(4), dimension(0:2) :: dra,dr1
integer :: j=0,ja=0
real(4), intent(out) :: costheta

dra(:) = 0.0
do j = 0, nbackbone-1
    ja = backbone(j)-1
    dr1(:) = 0.0
    call pbcdr(dr1,allcoord(ja,:),headcoord(:),box,binv)
    if ( NORM2(dr1) >= NORM2(dra) ) then
        dra(:) = dr1(:)
    end if
end do
costheta = 0.0
costheta = (dra(2)*dra(2))/(dra(0)*dra(0)+dra(1)*dra(1)+dra(2)*dra(2))
costheta = SIGN(sqrt(costheta),dra(2))
if (costheta < -1 .OR. costheta > 1) then
    print *, "Unusual dipole orientation &
    & wrt surface normal calculated for", &
    & "=", costheta
end if
end subroutine endtoend_orient


subroutine average_loc(avgr,allcoord,ridarr,marr,narr,nallcoord)
implicit none
integer, intent(in) :: narr,nallcoord
real(4), intent(in), dimension(0:nallcoord-1,0:2) :: allcoord
integer, intent(in), dimension(0:narr-1) :: ridarr
real(4), intent(in), dimension(0:narr-1) :: marr
real(4), intent(out), dimension(0:2) :: avgr
real(4), dimension(0:2) :: dra
integer :: i = 0, j = 0 , ja = 0
dra(:) = 0.0

do j = 0, narr-1

ja = ridarr(j)-1
avgr(:) = avgr(:) + marr(j)*allcoord(ja,:)

end do

avgr = avgr/(SUM(marr))

end subroutine average_loc
