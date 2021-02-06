subroutine read_header(tot_atom, tot_frame, box, nx, ny, nz, dx, dy, dz, filename)
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
!  Number of grid points
implicit none
character(len=200), intent(in) :: filename
character(len=4) :: head
!f2py intent(out) :: head
integer, dimension(1:18) :: dat
!real(4), dimension(4) :: dat1
real(8), dimension(6) :: box5
real(8), intent(out), dimension(0:2) :: box
!f2py intent(out) :: boxLx,boxLy,boxLz
real(4), intent(out) :: dx, dy, dz
integer, intent(out) :: tot_atom, tot_frame, nx, ny, nz
!f2py intent(out) :: tot_atom, tot_frame, nx, ny, nz
integer :: i=0, j=0, k=0
real :: a=0.0, b=0.0,delta=0.0
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
!Computing number of volume cells to initialize size apriori
!Each cell is of the size ~ 1/5 Angstrom
dx = 0.2
dy = 0.2
dz = 0.3
nx = anint((box(0))/dx)+2
ny = anint((box(1))/dy)+2
nz = anint((box(2))/dz)+2
print *, nx,ny,nz,dx,dy,dz
close(unit = 4000)
end subroutine read_header
