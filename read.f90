subroutine opendcd(filenum,filename)
integer, intent(in) :: filenum
character(LEN=200), intent(in) :: filename
open(unit = filenum, file = trim(filename), form = "UNFORMATTED", status = "OLD")
end subroutine opendcd

subroutine closedcd(filenum)
integer, intent(in) :: filenum
close(filenum)
end subroutine closedcd


subroutine skipframedcd(filenum)
integer, intent(in) :: filenum
read(UNIT = filenum) 
read(UNIT = filenum) 
read(UNIT = filenum) 
read(UNIT = filenum) 
end subroutine skipframedcd

subroutine skipheaderdcd(filenum)
integer, intent(in) :: filenum
read(UNIT = filenum) 
read(UNIT = filenum) 
read(UNIT = filenum) 
end subroutine skipheaderdcd

subroutine readheaderdcd(tot_atom, tot_frame, box, filename)
! Input: 
!  Pathname of dcd file
! Units: Length (nm)
! Read a dcd file and gather/do:
!  Header: Name of simulation
!  Total number of frames in dcd
!  Total number of atoms in dcd
!  Box length vectors in dcd
!  Create a grid with cells of size 0.05 * units along the box vectors
! Output: 
!  Header
!  Number of atoms
!  Number of frames
!  Box vectors
implicit none
character(len=200), intent(in) :: filename
character(len=4) :: head
!f2py intent(out) :: head
integer, dimension(1:18) :: dat
!real(4), dimension(4) :: dat1
real(8), dimension(6) :: box5
real(8), intent(out), dimension(0:2) :: box
real :: delta
!f2py intent(out) :: boxLx,boxLy,boxLz
integer, intent(out) :: tot_atom, tot_frame
!f2py intent(out) :: tot_atom, tot_frame, nx, ny, nz
integer :: i=0
dat(:) = 0.0
print *, "read_header is ",filename
open(unit = 4000, file = trim(filename), form = "UNFORMATTED", status ="OLD")
read(4000) head, tot_frame, (dat(i), i=1,8), delta, (dat(i), i=9,18)
read(4000) 
read(4000) tot_atom
read(4000) box5
box(0) = box5(1)
box(1) = box5(3)
box(2) = box5(6)
close(unit = 4000)
end subroutine readheaderdcd

subroutine readframedcd(rcoord,box,filenum,tot_atom)
integer, intent(in) :: filenum
integer, intent(in) :: tot_atom
real(8), dimension(6) :: box6
real(8), intent(out), dimension(0:2) :: box
real(4), intent(out), dimension(0:tot_atom-1,0:2) :: rcoord
read(UNIT = filenum) box6
read(UNIT = filenum) rcoord(:,0)
read(UNIT = filenum) rcoord(:,1)
read(UNIT = filenum) rcoord(:,2)
box(0) = box6(1)
box(1) = box6(3)
box(2) = box6(6)
end subroutine readframedcd
