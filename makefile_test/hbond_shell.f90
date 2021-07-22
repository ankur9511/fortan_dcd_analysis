subroutine nhbond_nshell(shell1, npairs, nshell1, &
    & centers,Hcenters1,Hcenters2, &
    & shell,Hshell1,Hshell2, &
    & r,box,binv,dist2,cosangle_crit, &
    & tot_atom,ncenters,nshell)

implicit none

!character*100, intent(in) :: criteria

integer, intent(in) :: tot_atom,ncenters,nshell
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(8), intent(in), dimension(0:2) :: box,binv
integer, intent(in), dimension(0:ncenters-1) :: centers
integer, intent(in), dimension(0:nshell-1) :: shell
real(4), intent(in) :: dist2 ! xl, xh, yl, yh, zl, zh
real(4), intent(in) :: cosangle_crit
integer, intent(in), dimension(0:ncenters-1) :: Hcenters1,Hcenters2
integer, intent(in), dimension(0:nshell-1) :: Hshell1,Hshell2
real(4), dimension(0:2) :: dr
real(4) :: dr2,costheta
integer, dimension(0:ncenters-1,0:1) :: Hcenters
integer, dimension(0:nshell-1,0:1) :: Hshell
integer :: i,j,k,flag
integer, intent(inout) :: npairs,nshell1
integer, intent(inout), dimension(0:nshell-1) :: shell1
npairs = 0
shell1 = -1
nshell1 = 0

Hshell = -1
Hshell(:,0) = Hshell1(:)
Hshell(:,1) = Hshell2(:)

Hcenters = -1
Hcenters(:,0) = Hcenters1(:)
Hcenters(:,1) = Hcenters2(:)

!print*, nshell,ncenters
!print*, shell

do j = 0,ncenters-1        
do i = 0,nshell-1   
flag = -1
if (centers(j)/=shell(i)) then

call pbcdr(dr,r(centers(j)-1,:), &
   & r(shell(i)-1,:),box,binv)
dr2 = dr(0)*dr(0) + dr(1)*dr(1) + dr(2)*dr(2)
costheta = 0.0

if (Hcenters(j,0)/=-1) then
costheta = 0.0
call cosangle(costheta, &
         & r(centers(j)-1,:), &
         & r(Hcenters(j,0)-1,:), &
         & r(shell(i)-1,:),box,binv)

if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
npairs = npairs+1
flag = 1
end if
end if

if (Hcenters(j,1)/=-1) then
costheta = 0.0
call cosangle(costheta, &
         & r(centers(j)-1,:), &
         & r(Hcenters(j,1)-1,:), &
         & r(shell(i)-1,:),box,binv)

if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
npairs = npairs+1
flag = 1
end if
end if

if (Hshell(i,0)/=-1) then
costheta = 0.0
call cosangle(costheta, &
         & r(shell(i)-1,:), &
         & r(Hshell(i,0)-1,:), &
         & r(centers(j)-1,:),box,binv)

if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
npairs = npairs+1
flag = 1
end if
end if

if (Hshell(i,1)/=-1) then
costheta = 0.0
call cosangle(costheta, &
        & r(shell(i)-1,:), &
        & r(Hshell(i,1)-1,:), &
        & r(centers(j)-1,:),box,binv)

if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
npairs = npairs+1
flag = 1
end if
end if

if (flag == 1) then
if ( .NOT. ANY(shell1==i)) then
shell1(nshell1) = i
nshell1 = nshell1+1
end if
end if
!if (j == ncenters1-1) then
!   if ( .NOT. ANY(shellc11==i)) then
!      notshellc11(nnshellc11) = i
!      nnshellc11 = nnshellc11+1
!   end if
!end if

end if
end do
end do
end subroutine nhbond_nshell