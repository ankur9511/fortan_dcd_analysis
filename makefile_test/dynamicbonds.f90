subroutine dynnc_bond_angle(npairs,centers,ncenters,shell,nshell,r,criteria,tot_atom,box,binv,dist2,cosangle_crit,H)
implicit none
character*100, intent(in) :: criteria
integer, intent(in) :: tot_atom,ncenters,nshell
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(8), intent(in), dimension(0:2) :: box,binv
integer, intent(in), dimension(0:ncenters-1) :: centers
!!!!! If angle criteria is present, center must contain all O's that have a bonded H
integer, intent(in), dimension(0:nshell-1) :: shell
real(4), intent(in) :: dist2
integer, optional, intent(in), dimension(0:ncenters-1,0:1) :: H
!!!!! If angle criteria is present, H must contain atomID of H's bonded to centers. Max is 2. If only 1 present, atomID is -1
real(4), optional, intent(in) :: cosangle_crit
real(4), dimension(0:2) :: dr
real(4) :: dr2,costheta
integer :: i,j,k,index1,index2
integer, intent(inout) :: npairs
npairs = 0
index1=index(criteria,'dist')
index2=index(criteria,'angle')

if (index1 * index2 > 0) then
  do j = 0,ncenters-1 !! j = Water or any other O in center
   do i = 0,nshell-1    !! i = Water in tetrahedral shell
    if (centers(j)/=shell(i)) then
     call pbcdr(dr,r(centers(j)-1,:),r(shell(i)-1,:),box,binv)
     dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2)
     costheta = 0.0
     if (H(j,0)/=-1) then
      costheta = 0.0
      call cosangle(costheta,r(centers(j)-1,:),r(H(j,0)-1,:),r(shell(i)-1,:),box,binv)
      if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
        npairs = npairs+1      
      end if
     end if
     if (H(j,1)/=-1) then
      costheta = 0.0
      call cosangle(costheta,r(centers(j)-1,:),r(H(j,1)-1,:),r(shell(i)-1,:),box,binv)
      if (dr2 <= dist2 .and. costheta >= cosangle_crit) then
        npairs = npairs+1      
      end if
     end if
    end if
   end do
  end do

else if (index1 > 0) then

  do j = 0,ncenters-1 !! j = Water or any other O in center 
    print *, "Inside iteration only dist",j
    do i = 0,nshell-1    !! i = Water in tetrahedral shell
      if (centers(j)/=shell(i)) then
        call pbcdr(dr,r(centers(j)-1,:),r(shell(i)-1,:),box,binv)
        dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2)
        if (dr2 <= dist2) then
           npairs = npairs+1      
        end if
      end if
    end do
    print *, "Iteration only dist ",j,"Ended"
   end do
 
end if

end subroutine dynnc_bond_angle
