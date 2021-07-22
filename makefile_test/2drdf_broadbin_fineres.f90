subroutine twodrdfsolv_pf(nsolv2drdfarr,tsolv2drdfarr,center_shell,r,dbinrdf,dbincrdf,lastsolv2drdf,firstsolv2drdf,box,binv,hbox0,ncenter_shell,tot_atom,lensolv2drdfarr,nbincrdf)
    implicit none
    integer, intent(in) :: ncenter_shell,tot_atom,lensolv2drdfarr,nbincrdf
    real(8), intent(in), dimension(0:2) :: box,binv,hbox0
    real(4), intent(in) :: dbincrdf,dbinrdf,lastsolv2drdf,firstsolv2drdf
    real(4), intent(in), dimension(0:tot_atom-1,0:2) :: r
    integer, intent(in), dimension(0:ncenter_shell-1) :: center_shell
    real(4), dimension(0:nbincrdf-1,0:lensolv2drdfarr-1), intent(inout) :: tsolv2drdfarr
    real(4), dimension(0:nbincrdf-1), intent(inout) :: nsolv2drdfarr
    real(4), dimension(0:2) :: dr1,dr2
    real(4) :: dr12,dr22
    integer :: i,j,k,l,index,zindexk,zindexl,flagkl,binctobin
    real(4), dimension(0:nbincrdf-1) :: zk,zl
    print *, "Inside 2D"    
    binctobin = int(0.5*dbinrdf/dbincrdf)
    do i = 0,ncenter_shell-2
        k = center_shell(i)-1
        zk(:) = 0.0
        zindexk = nint((r(k,2)+hbox0(2))*1.0/dbincrdf)
        zk(zindexk-binctobin:zindexk+binctobin) = 1.0
        do j = i+1,ncenter_shell-1
            l = center_shell(j)-1
            zl(:)=0.0
            zindexl = nint((r(l,2)+hbox0(2))*1.0/dbincrdf)
            zl(zindexl-binctobin:zindexl+binctobin) = 1.0
            call pbcdr(dr1,r(l,:),r(k,:),box,binv)
            dr12 = dot_product(dr1(0:1),dr1(0:1))
            nsolv2drdfarr(:) = nsolv2drdfarr(:)+2.0*zl(:)*zk(:)
            index = nint((sqrt(dr12) - (firstsolv2drdf))*float(lensolv2drdfarr)/(lastsolv2drdf-firstsolv2drdf))
            if ( index <= lensolv2drdfarr-1 ) then
             tsolv2drdfarr(:,index) = tsolv2drdfarr(:,index)+2.0*zl(:)*zk(:)
            end if
        end do
    end do
end subroutine twodrdfsolv_pf
