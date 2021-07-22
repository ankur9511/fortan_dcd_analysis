SUBROUTINE  llst_icell(ibox,rxi,ryi,rzi,icell)

    USE global
    IMPLICIT NONE 
  
    ! Passed
    INTEGER :: ibox, imol, icell
  
    ! Local
    DOUBLE PRECISION :: icellx, icelly, icellz
    DOUBLE PRECISION :: rxi,ryi,rzi
  
    ! Zero icell for good measure
    icell = 0
  
    ! Calculate the cell number in each direction
    icellx = INT((rxi + 0.5d0*box(1,ibox))/dcell(1,ibox))
    icelly = INT((ryi + 0.5d0*box(2,ibox))/dcell(2,ibox))
    icellz = INT((rzi + 0.5d0*box(3,ibox))/dcell(3,ibox))
  
    icell = 1 + icellx + (icelly*n_cells(1,ibox)) + (icellz*n_cells(1,ibox)*n_cells(2,ibox))
  
    IF ((icell .LT. 1) .OR. (icell .GT. n_cells_tot(ibox))) THEN
      WRITE(*,'(A,I4,A)') 'FATAL ERROR IN LLST_ICELL: IN BOX ', ibox, ' IS OUT OF BOUNDS'
      WRITE(*,'(A,I4,3F20.8)') 'ICELL: ', icell,rxi,ryi,rzi
      STOP
    ENDIF 
  
  
    RETURN
  
  END SUBROUTINE llst_icell
  