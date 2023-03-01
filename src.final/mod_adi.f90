MODULE mod_ADI
CONTAINS
!

  SUBROUTINE ADI_VXd(ux,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx,dz
    USE inv_mat_band, ONLY : tridiag,cyclic,cyclic_nr,tridiag_nr
    USE mod_bounds  , ONLY : bounds_VX
    USE tableaux, ONLY : Zvx, w
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,kdttdx2,kdttdz2
    REAL(KIND=8) :: lpp0, lpp1, lpp2
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm,sol
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: ux
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    kdttdx2=diff*dt/2./dx**2
    kdttdz2=diff*dt/2./dz**2
    da=0. ; db=0. ; dc=0. ; sm=0.

!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(kdttdx2,Zvx,ux,rhs,w,dt,diff)
    DO iz=1,nz
       da(1:nx)=    -kdttdx2
       db(1:nx)=1.+2.*kdttdx2
       dc(1:nx)=da(1:nx)

        lpp0=2/((Zvx(iz-1)-Zvx(iz))*(Zvx(iz-1)-Zvx(iz+1)))
        lpp1=2/((Zvx(iz)-Zvx(iz-1))*(Zvx(iz)-Zvx(iz+1)))
        lpp2=2/((Zvx(iz+1)-Zvx(iz-1))*(Zvx(iz+1)-Zvx(iz)))

        sm(1:nx)= lpp0*ux(1:nx, iz-1)*diff*dt/2. &
                +(1.+diff*dt/2.*lpp1)*ux(1:nx, iz) &
                +lpp2*ux(1:nx, iz+1)*diff*dt/2. &
                +rhs(1:nx,iz)*dt/2.

       CALL cyclic_nr(da,db,dc,dc(nx),da(1),sm,sol,nx)
       w(1:nx,iz)=sol(1:nx)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VX(w)
!
    da=0. ; db=0. ; dc=0. ; sm=0.
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(kdttdx2,Zvx,ux,rhs,w,dt,diff)
    DO ix=1,nx
            da(1:nz)=   -diff*dt/2.*(2./((Zvx(0:nz-1)-Zvx(1:nz))*(Zvx(0:nz-1)-Zvx(2:nz+1))))
            db(1:nz)= 1.-diff*dt/2.*(2./((Zvx(1:nz)-Zvx(0:nz-1))*(Zvx(1:nz)-Zvx(2:nz+1))))
            dc(1:nz)=   -diff*dt/2.*(2./((Zvx(2:nz+1)-Zvx(0:nz-1))*(Zvx(2:nz+1)-Zvx(1:nz))))

       sm(1:nz)=w(ix,1:nz)+(w(ix+1,1:nz)-2.*w(ix,1:nz)+w(ix-1,1:nz))*kdttdx2 &
                 +rhs(ix,1:nz)*dt/2.
       db(1)=db(1)+da(1) ; da(1)=0.
       db(nz)=db(nz)+dc(nz) ; dc(nz)=0.
       CALL tridiag_nr(da,db,dc,sm,sol,nz)
       ux(ix,1:nz)=sol(1:nz)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VX(ux)
!
  END SUBROUTINE ADI_VXd
!
!
  SUBROUTINE ADI_VZd(uz,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx,dz
    USE inv_mat_band, ONLY : tridiag,cyclic,cyclic_nr,tridiag_nr
    USE mod_bounds  , ONLY : bounds_VZ
    USE tableaux,     ONLY: Zvz,w
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,kdttdx2,kdttdz2
    REAL(KIND=8) :: lpp0, lpp1, lpp2
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm,sol
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: uz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    kdttdx2=diff*dt/2./dx**2
    kdttdz2=diff*dt/2./dz**2
    da=0. ; db=0. ; dc=0. ; sm=0.
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(kdttdx2,Zvz,uz,rhs,w,dt,diff)
    DO iz=1,nz+1
        lpp0=2/((Zvz(iz-1)-Zvz(iz))*(Zvz(iz-1)-Zvz(iz+1)))
        lpp1=2/((Zvz(iz)-Zvz(iz-1))*(Zvz(iz)-Zvz(iz+1)))
        lpp2=2/((Zvz(iz+1)-Zvz(iz-1))*(Zvz(iz+1)-Zvz(iz)))

        da(1:nx)=    -kdttdx2
        db(1:nx)=1.+2.*kdttdx2
        dc(1:nx)=da(1:nx)
        
            sm(1:nx)=lpp0*uz(1:nx, iz-1)*diff*dt/2. &
                    +(1.+diff*dt/2.*lpp1)*uz(1:nx, iz) &
                    +lpp2*uz(1:nx, iz+1)*diff*dt/2. &
                    +rhs(1:nx,iz)*dt/2.

       CALL cyclic_nr(da,db,dc,dc(nx),da(1),sm,sol,nx)
       w(1:nx,iz)=sol(1:nx)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VZ(w)
!
    da=0. ; db=0. ; dc=0. ; sm=0.
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(kdttdx2,Zvz,uz,rhs,w,dt,diff)

    DO ix=1,nx+1
            da(1:nz+1)=    -diff*dt/2.*(2./((Zvz(0:nz)-Zvz(1:nz+1))*(Zvz(0:nz)-Zvz(2:nz+2))))
            db(1:nz+1)=  1.-diff*dt/2.*(2./((Zvz(1:nz+1)-Zvz(0:nz))*(Zvz(1:nz+1)-Zvz(2:nz+2))))
            dc(1:nz+1)=    -diff*dt/2.*(2./((Zvz(2:nz+2)-Zvz(0:nz))*(Zvz(2:nz+2)-Zvz(1:nz+1))))

       sm(1:nz+1)=w(ix,1:nz+1)+(w(ix+1,1:nz+1)-2.*w(ix,1:nz+1)+w(ix-1,1:nz+1))*kdttdx2 &
                 +rhs(ix,1:nz+1)*dt/2.
       db(1)=1. ; da(1)=0. ; dc(1)=0. ; sm(1)=0.
       da(2)=0.
       dc(nz)=0.
       db(nz+1)=1. ; da(nz+1)=0. ; dc(nz+1)=0. ; sm(nz+1)=0.
       CALL tridiag_nr(da,db,dc,sm,sol,nz+1)
       uz(ix,1:nz+1)=sol(1:nz+1)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VZ(uz)
!
  END SUBROUTINE ADI_VZd
!
!
END MODULE mod_ADI
