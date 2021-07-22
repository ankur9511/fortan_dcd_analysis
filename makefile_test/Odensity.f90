subroutine Odensity(Orz,r,shell_ID,centers_ID,box,binv,dz,nzdd,n_shell_ID,n_centers_ID,tot_atom)
implicit none
integer, intent(in) :: tot_atom,n_shell_ID,n_centers_ID,nzdd
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(4), intent(in) :: dz
integer, intent(in), dimension(0:n_centers_ID-1) :: centers_ID
integer, intent(in), dimension(0:n_shell_ID-1) :: shell_ID
real(8), intent(in), dimension(0:2) :: box,binv
real(4), intent(inout), dimension(0:nzdd-1) :: Orz
real(4), dimension(0:2) :: minr,tempr
real(4) :: minr2,tempr2,lastOrz
real(4), dimension(0:n_centers_ID-1,0:2) :: dr
real(4), dimension(0:n_centers_ID-1) :: dr2
integer :: i,j,k,l,m,shellelement
character systime*8
minr(:) = 1000.0
minr2 = 10000.0
tempr(:) = 0.0
lastOrz = nzdd*dz
do k = 0,n_shell_ID-1
  shellelement = shell_ID(k)
  dr(:,0) = r(centers_ID(:)-1,0) - r(shellelement-1,0)
  dr(:,1) = r(centers_ID(:)-1,1) - r(shellelement-1,1)
  dr(:,2) = r(centers_ID(:)-1,2) - r(shellelement-1,2)
  dr(:,0) = dr(:,0) - box(0)*anint(dr(:,0)*binv(0))
  dr(:,1) = dr(:,1) - box(1)*anint(dr(:,1)*binv(1))
  dr(:,2) = dr(:,2) - box(2)*anint(dr(:,2)*binv(2))
  dr2(:) = dr(:,0)*dr(:,0)+dr(:,1)*dr(:,1)+dr(:,2)*dr(:,2)
  minr2 = sqrt(minval(dr2))
 
  j = nint((minr2 - (0))*nzdd/lastOrz)
  if (j < nzdd) then    
    Orz(j) = Orz(j)+1.0
  end if
  !if (MOD(k,2000)==0) then
  !  call time(systime)
  !  i = int(float(k-0)/float(n_shell_ID-0)*100)
  !  print *, "AT time:",systime," loop index:",k," Progress of frames:",i,"%"
  !end if
end do
end subroutine Odensity
