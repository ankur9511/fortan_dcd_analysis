subroutine pbcwrapxy(r,rcid,tot_atom,box,binv)
implicit none
integer, intent(in) :: tot_atom
real(8), intent(in), dimension(0:2) :: box,binv
real(4), intent(in), dimension(0:2) :: rcid
real(4), dimension(0:tot_atom-1,0:2), intent(inout) :: r

r(:,0) = r(:,0) - rcid(0) ! - box(0)*anint((r(:,0)-rcid(0))*binv(0))
r(:,1) = r(:,1) - rcid(1) ! - box(1)*anint((r(:,1)-rcid(1))*binv(1))
r(:,2) = r(:,2) - rcid(2) ! - box(2)*anint((r(:,2)-rcid(2))*binv(2))
end subroutine pbcwrapxy
