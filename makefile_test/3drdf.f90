subroutine threedrdfsolv_pf(nsolv3drdfarr,tsolv3drdfarr,center_ID,shell_ID,r,lastsolv3drdf,firstsolv3drdf,box,binv,hbox0,ncenter,nshell,tot_atom,lensolv3drdfarr)
    implicit none
    integer, intent(in) :: ncenter,nshell,tot_atom,lensolv3drdfarr
    real(8), intent(in), dimension(0:2) :: box,binv,hbox0
    real(4), intent(in) :: lastsolv3drdf,firstsolv3drdf
    real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
    integer, intent(in), dimension(0:ncenter-1) :: center_ID
    integer, intent(in), dimension(0:nshell-1) :: shell_ID
    real(4), dimension(0:lensolv3drdfarr-1), intent(inout) :: tsolv3drdfarr
    real(4), intent(inout) :: nsolv3drdfarr
    real(4), dimension(0:2) :: dr1,dr2
    real(4) :: dr12,dr22
    integer :: i,j,k,l,index,zindexk,zindexl,flagkl
    !real(4), dimension(0:nbincrdf-1) :: zk,zl
    !print *, "Inside 3D"    
    !binctobin = int(0.5*dbinrdf/dbincrdf)
    do i = 0,ncenter-1
        k = center_ID(i)-1
        !zk(:) = 0.0
        !zindexk = nint((r(k,2)+hbox0(2))*1.0/dbincrdf)
        !zk(zindexk-binctobin:zindexk+binctobin) = 1.0
        do j = 0,nshell-1
            l = shell_ID(j)-1
            if (l .NE. k) then
            !zl(:)=0.0
            !zindexl = nint((r(l,2)+hbox0(2))*1.0/dbincrdf)
            !zl(zindexl-binctobin:zindexl+binctobin) = 1.0
            call pbcdr(dr1,r(l,:),r(k,:),box,binv)
            dr12 = dot_product(dr1,dr1)
            nsolv3drdfarr = nsolv3drdfarr+1.0
            index = nint((sqrt(dr12) - (firstsolv3drdf))*float(lensolv3drdfarr)/(lastsolv3drdf-firstsolv3drdf))
            if ( index <= lensolv3drdfarr-1 ) then
             tsolv3drdfarr(index) = tsolv3drdfarr(index)+1.0
            end if
            end if
        end do
    end do
end subroutine threedrdfsolv_pf
