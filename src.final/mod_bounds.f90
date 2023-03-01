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
END MODULE mod_bounds
