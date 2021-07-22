subroutine voldensity(Arz,r,shell_ID,centers_ID,box,binv,dx,dy,dz,nx,ny,nz,nzdd,n_shell_ID,n_centers_ID,tot_atom)
implicit none
integer, intent(in) :: tot_atom,n_shell_ID,n_centers_ID,nx,ny,nz,nzdd
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(4), intent(in) :: dx,dy,dz
integer, intent(in), dimension(0:n_centers_ID-1) :: centers_ID
integer, intent(in), dimension(0:n_shell_ID-1) :: shell_ID
real(8), intent(in), dimension(0:2) :: box,binv
real(4), intent(inout), dimension(0:nzdd-1) :: Arz
real(4), dimension(0:2) :: minr,tempr
real(4) :: minr2,tempr2,x0,y0,z0,lastOrz
real(4), dimension(0:n_centers_ID-1,0:2) :: dr
real(4), dimension(0:n_centers_ID-1) :: dr2
real :: T1,T2
integer :: i,j,k,l,m,shellelement,ompid
character systime*8
minr(:) = 1000.0
minr2 = 10000.0
tempr(:) = 0.0
x0 = -0.5*(box(0))
y0 = -0.5*(box(1))
z0 = -0.5*(box(2))
lastOrz = nzdd*dz

do l = 0,nx*ny*nz-1
    dr(:,0) = r(centers_ID(:)-1,0) - x0+MOD(int(l/(ny*nz)),nx)*dx
    dr(:,1) = r(centers_ID(:)-1,1) - y0+MOD(int(l/(nz)),ny)*dy
    dr(:,2) = r(centers_ID(:)-1,2) - z0+MOD(int(l),nz)*dz
    dr(:,0) = dr(:,0) - box(0)*anint(dr(:,0)*binv(0))
    dr(:,1) = dr(:,1) - box(1)*anint(dr(:,1)*binv(1))
    dr(:,2) = dr(:,2) - box(2)*anint(dr(:,2)*binv(2))
    dr2(:) = dr(:,0)*dr(:,0)+dr(:,1)*dr(:,1)+dr(:,2)*dr(:,2)
    minr2 = sqrt(minval(dr2))
    m = nint((minr2 - (0))*nzdd/lastOrz)
    if (m < nzdd) then
    Arz(m) = Arz(m)+1.0
    end if
    if (MOD(l,1000000)==0) then
    call time(systime)
    i = int(float(l-0)/float(nx*ny*nz-1-0)*100)
    print *, "AT time:",systime," loop index:",l," Progress of frames:",i,"%"
    end if
end do

end subroutine voldensity
