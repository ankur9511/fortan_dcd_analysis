subroutine densfluc(occ,occkappa,centers,ncenters,shell,nshell,r,tot_atom,box,binv)
implicit none
integer, intent(in) :: tot_atom,ncenters,nshell
real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
real(4), dimension(0:2) :: rcenter
real(8), intent(in), dimension(0:2) :: box,binv
integer, intent(in), dimension(0:ncenters-1) :: centers
integer, intent(in), dimension(0:nshell-1) :: shell
real(4), intent(inout),dimension(0:2*ncenters-1) :: occ,occkappa
real(4), dimension(0:2*ncenters-1) :: tocc
real(4),dimension(0:2) :: dr
real(4) :: dr2,ctocdist,methsize,halfmethsizesq
integer :: i,j,k

tocc(:) = 0.0

methsize = 3.3 ! Angstrom
halfmethsizesq = 0.25*methsize*methsize
ctocdist = 0.5*(2.96568+methsize) ! 2.96568 is obtained from www.pnas.org/cgi/doi/10.1073/pnas.0809029105
do j = 0,ncenters-1 !! j = Water in center 
  do i = 0,nshell-1    !! i = Water in tetrahedral shell
    rcenter(:) = r(centers(j)-1,:)
    rcenter(2) = rcenter(2) + ctocdist
    call pbcdr(dr,rcenter(:),r(shell(i)-1,:),box,binv)
    dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2)

    if (dr2 <= halfmethsizesq) then
       occ(j) = occ(j) + 1.0
       tocc(j) = tocc(j) + 1.0
    end if

  end do
end do

ctocdist = 25+0.5*(methsize) ! 2.96568 is obtained from www.pnas.org/cgi/doi/10.1073/pnas.0809029105
do j = 0,ncenters-1 !! j = Water in center 
  do i = 0,nshell-1    !! i = Water in tetrahedral shell
    rcenter(:) = r(centers(j)-1,:)
    rcenter(2) = rcenter(2) + ctocdist
    call pbcdr(dr,rcenter(:),r(shell(i)-1,:),box,binv)
    dr2 = dr(0)*dr(0)+dr(1)*dr(1)+dr(2)*dr(2)

    if (dr2 <= halfmethsizesq) then
       occ(j+ncenters) = occ(j+ncenters) + 1.0
       tocc(j+ncenters) = tocc(j+ncenters) + 1.0
    end if

  end do
end do

occkappa(:) = occkappa(:) + tocc(:)*tocc(:)

end subroutine densfluc
