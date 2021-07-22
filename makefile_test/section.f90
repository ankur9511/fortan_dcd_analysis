subroutine rsection(sectionlist,rmin,rmax,rdirection, &
                    & box,binv, &
                    & rid,rcoord, &
                    & nid,tot_atom)
implicit none
integer, intent(in) :: nid,tot_atom,rdirection
integer, intent(in), dimension(0:nid-1) :: rid
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: rcoord
real(4), intent(in) :: rmin,rmax
real(8), intent(in), dimension(0:2) :: box,binv
integer, intent(out), dimension(0:nid-1) :: sectionlist
!integer, external :: pbcdrdir
integer :: i=0
real(4) :: dr1,dr2,dr3
sectionlist(:) = 0

!do i = 0,nid-1
!dr1 = rcoord(rid(i)-1,rdirection) - rmin
!dr1 = dr1 - box(rdirection)*anint(dr1*binv(rdirection))
!dr2 = rcoord(rid(i)-1,rdirection) - rmax
!dr2 = dr2 - box(rdirection)*anint(dr2*binv(rdirection))
!dr3 = dr1*dr2
!sectionlist(i) = int(SIGN(1.,(rcoord(rid(i)-1,rdirection) - rmin)*(rcoord(rid(i)-1,rdirection) - rmax)))
!end do
sectionlist(:) = int(SIGN(1.,(rcoord(rid(:)-1,rdirection) - rmin)*(rcoord(rid(:)-1,rdirection) - rmax)))
!
end subroutine rsection

!integer function pbcdrdir(dirmin,dirmax,val,box,binv)
!implicit none
!real(4), intent(in) :: dirmin,dirmax,val,box,binv
!real(4) :: dr1=0., dr2=0., dr3=0.

!dr1 = val - dirmin
!dr1 = dr1 - box*anint(dr1*binv)
!dr2 = val - dirmax
!dr2 = dr2 - box*anint(dr2*binv)
!dr3 = dr1*dr2
!pbcdrdir = int(SIGN(1.,dr3))
!return
!end 


