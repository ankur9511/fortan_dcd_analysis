subroutine pbcind(ind,indin,maxind)
integer, intent(in) :: maxind
integer, intent(in) :: indin
integer, intent(inout) :: ind

!dr = r1 - r2
ind = ind - maxind*anint(ind/float(maxind))

end subroutine pbcind
