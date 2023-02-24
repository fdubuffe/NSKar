MODULE mod_bounds
CONTAINS
!

  SUBROUTINE bounds_PR(w)
!=====================================================================
  USE constantes, ONLY : nx,nz,Ray,dz,Tt,Tb,T0
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
!=====================================================================
!    w(-1:nx+3,  -1)=-0.5*Ray*2.*dz+w(-1:nx+3,   2)
!    w(-1:nx+3,   0)=-0.5*Ray*1.*dz+w(-1:nx+3,   1)
!    w(-1:nx+3,nz+1)=-0.5*Ray*1.*dz+w(-1:nx+3,  nz)
!    w(-1:nx+3,nz+2)=-0.5*Ray*2.*dz+w(-1:nx+3,nz-1)
    w(-1:nx+3,  -1)=+Ray*(Tt-T0)*2.*dz+w(-1:nx+3,   2)
    w(-1:nx+3,   0)=+Ray*(Tt-T0)*1.*dz+w(-1:nx+3,   1)
    w(-1:nx+3,nz+1)=-Ray*(Tb-T0)*1.*dz+w(-1:nx+3,  nz)
    w(-1:nx+3,nz+2)=-Ray*(Tb-T0)*2.*dz+w(-1:nx+3,nz-1)
    w(-1,-1:nz+2)=w(nx-1,-1:nz+2)
    w( 0,-1:nz+2)=w(nx  ,-1:nz+2)
    w(nx+1,-1:nz+2)=w(1,-1:nz+2)
    w(nx+2,-1:nz+2)=w(2,-1:nz+2)
    w(nx+3,-1:nz+2)=w(3,-1:nz+2)
  END SUBROUTINE bounds_PR
!
!
  SUBROUTINE bounds_TP(w)
!=====================================================================
  USE constantes, ONLY : nx,nz,Tt,Tb
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
!=====================================================================
    w(1:nx+1,  -1)=2.*Tt-w(1:nx+1,   2)
    w(1:nx+1,   0)=2.*Tt-w(1:nx+1,   1)
    w(1:nx+1,nz+1)=2.*Tb-w(1:nx+1,nz  )
    w(1:nx+1,nz+2)=2.*Tb-w(1:nx+1,nz-1)
    w(-1  ,-1:nz+2)=w(nx-1,-1:nz+2)
    w( 0  ,-1:nz+2)=w(nx  ,-1:nz+2)
    w(nx+1,-1:nz+2)=w(1   ,-1:nz+2)
    w(nx+2,-1:nz+2)=w(2   ,-1:nz+2)
  END SUBROUTINE bounds_TP
!
  SUBROUTINE bounds_VX(w)
!=====================================================================
  USE constantes, ONLY : nx,nz
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
!=====================================================================
    w(1:nx,-1  )=w(1:nx,2   )
    w(1:nx, 0  )=w(1:nx,1   )
    w(1:nx,nz+1)=w(1:nx,nz  )
    w(1:nx,nz+2)=w(1:nx,nz-1)
    w(-1  ,-1:nz+2)=w(nx-1,-1:nz+2)
    w( 0  ,-1:nz+2)=w(nx  ,-1:nz+2)
    w(nx+1,-1:nz+2)=w(1   ,-1:nz+2)
    w(nx+2,-1:nz+2)=w(2   ,-1:nz+2)
  END SUBROUTINE bounds_VX
!
!
  SUBROUTINE bounds_VZ(w)
!=====================================================================
  USE constantes, ONLY : nx,nz
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
!=====================================================================
    w(1:nx+1,-1)=-w(1:nx+1,3)
    w(1:nx+1, 0)=-w(1:nx+1,2)
    w(1:nx+1, 1)=0.
    w(1:nx+1,nz+1)=0.
    w(1:nx+1,nz+2)=-w(1:nx+1,nz  )
    w(1:nx+1,nz+3)=-w(1:nx+1,nz-1)
    w(-1  ,-1:nz+3)=w(nx-1,-1:nz+3)
    w( 0  ,-1:nz+3)=w(nx  ,-1:nz+3)
    w(nx+1,-1:nz+3)=w(1   ,-1:nz+3)
    w(nx+2,-1:nz+3)=w(2   ,-1:nz+3)
  END SUBROUTINE bounds_VZ
!
!
  SUBROUTINE bounds_cond(w,cl)
!=====================================================================
  USE constantes, ONLY : nx,nz,Pra,Ray,dz
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
  CHARACTER(LEN=1), DIMENSION(2) :: cl
!=====================================================================
!
  IF ( cl(1) == 'S') THEN
    w(  -1,-1:nz+3)=w(   3,-1:nz+3)
    w(   0,-1:nz+3)=w(   2,-1:nz+3)
    w(nx+2,-1:nz+3)=w(  nx,-1:nz+3)
    w(nx+3,-1:nz+3)=w(nx-1,-1:nz+3)
  ELSEIF ( cl(1) == 'A') THEN
    w(  -1,-1:nz+3)=-w(   3,-1:nz+3)
    w(   0,-1:nz+3)=-w(   2,-1:nz+3)
    w(   1,-1:nz+3)=0.
    w(nx+1,-1:nz+3)=0.
    w(nx+2,-1:nz+3)=-w(  nx,-1:nz+3)
    w(nx+3,-1:nz+3)=-w(nx-1,-1:nz+3)
  ELSEIF ( cl(1) == 'P') THEN
    w(  -1,-1:nz+3)=w(nx-1,-1:nz+3)
    w(   0,-1:nz+3)=w(  nx,-1:nz+3)
    w(nx+1,-1:nz+3)=w(   1,-1:nz+3)
    w(nx+2,-1:nz+3)=w(   2,-1:nz+3)
    w(nx+3,-1:nz+3)=w(   3,-1:nz+3)
  ELSE
    STOP 'Horizontal voundary conditions not implemented'
  ENDIF
!
  IF ( cl(2) == 'S') THEN
    w(-1:nx+3,  -1)=w(-1:nx+3,   3)
    w(-1:nx+3,   0)=w(-1:nx+3,   2)
    w(-1:nx+3,nz+2)=w(-1:nx+3,  nz)
    w(-1:nx+3,nz+3)=w(-1:nx+3,nz-1)
  ELSEIF ( cl(2) == 'A') THEN
    w(-1:nx+3,  -1)=-w(-1:nx+3,   3)
    w(-1:nx+3,   0)=-w(-1:nx+3,   2)
    w(-1:nx+3,   1)=0.
    w(-1:nx+3,nz+1)=0.
    w(-1:nx+3,nz+2)=-w(-1:nx+3,  nz)
    w(-1:nx+3,nz+3)=-w(-1:nx+3,nz-1)
  ELSEIF ( cl(2) == 'P') THEN
    w(-1:nx+3,  -1)=w(-1:nx+3,nz-1)
    w(-1:nx+3,   0)=w(-1:nx+3,  nz)
    w(-1:nx+3,nz+1)=w(-1:nx+3,   1)
    w(-1:nx+3,nz+2)=w(-1:nx+3,   2)
    w(-1:nx+3,nz+3)=w(-1:nx+3,   3)
  ELSEIF ( cl(2) == 'I') THEN
    w(-1:nx+3,   1)=0.
    w(-1:nx+3,   0)=2.*w(-1:nx+3,   1)-w(-1:nx+3,   2)
    w(-1:nx+3,  -1)=2.*w(-1:nx+3,   1)-w(-1:nx+3,   3)
    w(-1:nx+3,nz+1)=1.
    w(-1:nx+3,nz+2)=2.*w(-1:nx+3,nz+1)-w(-1:nx+3,  nz)
    w(-1:nx+3,nz+3)=2.*w(-1:nx+3,nz+1)-w(-1:nx+3,nz-1)
  ELSEIF ( cl(2) == 'J') THEN
    w(-1:nx+3,  -1)=w(-1:nx+3,   3)
    w(-1:nx+3,   0)=w(-1:nx+3,   2)
    w(-1:nx+3,nz+2)=-Pra*Ray*2.*dz+w(-1:nx+3,  nz)
    w(-1:nx+3,nz+3)=-Pra*Ray*4.*dz+w(-1:nx+3,nz-1)
  ELSE
    STOP 'Vertical boundary conditions not implemented'
  ENDIF
!
  END SUBROUTINE bounds_cond
!
END MODULE mod_bounds
