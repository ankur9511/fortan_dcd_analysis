subroutine findfournn(neighborID,centers,ncenters,shell,nshell,r,tot_atom,box,binv)
implicit none
integer, intent(in) :: tot_atom,ncenters,nshell
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(8), intent(in), dimension(0:2) :: box,binv
integer, intent(in), dimension(0:ncenters-1) :: centers
integer, intent(in), dimension(0:nshell-1) :: shell
integer, intent(inout), dimension(0:ncenters-1,0:3) :: neighborID
real(4),dimension(0:2) :: dr
real(4) :: dr2
integer :: i,j,k
real(4), dimension(0:ncenters-1,0:3) :: neighbor

neighbor(:,:) = 1000.0
neighborID(:,:) = -1 

do j = 0,ncenters-1 !! j = Water in center 
  do i = 0,nshell-1    !! i = Water in tetrahedral shell
   if (centers(j)/=shell(i)) then
    call pbcdr(dr,r(centers(j)-1,:),r(shell(i)-1,:),box,binv)
    dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2)
    do k =0,3
     if (dr2 < neighbor(j,k)) then
      neighbor(j,k)=dr2
      neighborID(j,k)=i
      EXIT
     else
      continue
     endif
    end do
   end if
  end do
 end do
end subroutine findfournn
