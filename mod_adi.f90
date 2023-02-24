MODULE mod_ADI
CONTAINS
!
  SUBROUTINE ADI_TP(tp,vx,vz,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx,dz,Tt,Tb
    USE inv_mat_band, ONLY : cyclic,tridiag
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,dtt,kdttdx2,kdttdz2,dttd2x,dttd2z
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: vx,vz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: tp
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    dtt=dt/2.
    kdttdx2=diff*dtt/dx**2 ; kdttdz2=diff*dtt/dz**2  
    dttd2x=dtt/dx/2. ; dttd2z=dtt/dz/2.
!
! Schema de Peaceman & Rashford, integration selon X
    w=0.
    da=0. ; db=0. ; dc=0. ; sm=0.
!
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)              &
!$OMP PRIVATE (da,db,dc,sm,iz) &
!$OMP SHARED(dtt,kdttdx2,kdttdz2,dttd2x,dttd2z,vx,vz,tp,w,rhs)
    DO iz=1,nz
       da(1:nx+1)=    -     kdttdx2               -vx(0:nx,iz) *dttd2x
       dc(1:nx+1)=    -     kdttdx2+ vx(1:nx+1,iz)             *dttd2x
       db(1:nx+1)=1.+2.*kdttdx2+(vx(1:nx+1,iz)-vx(0:nx,iz))*dttd2x
       sm(1:nx+1)=tp(1:nx+1,iz)+(tp(1:nx+1,iz+1)-2.*tp(1:nx+1,iz)+tp(1:nx+1,iz-1))*kdttdz2                             &
            -(vz(1:nx+1,iz+1)*(tp(1:nx+1,iz+1)+tp(1:nx+1,iz))-vz(1:nx+1,iz)*(tp(1:nx+1,iz)+tp(1:nx+1,iz-1)))*dttd2z &
            +rhs(1:nx+1,iz)*dtt
       CALL cyclic(da,db,dc,dc(nx),da(1),sm,nx)
       w(1:nx+1,iz)=sm(1:nx+1)
    ENDDO
!$OMP END PARALLEL DO
    w(1:nx+1,  -1)=2.*Tt-w(1:nx+1,   2)
    w(1:nx+1,   0)=2.*Tt-w(1:nx+1,   1)
    w(1:nx+1,nz+1)=2.*Tb-w(1:nx+1,nz  )
    w(1:nx+1,nz+2)=2.*Tb-w(1:nx+1,nz-1)
    w(-1  ,-1:nz+2)=w(nx-1,-1:nz+2)
    w( 0  ,-1:nz+2)=w(nx  ,-1:nz+2)
    w(nx+1,-1:nz+2)=w(1   ,-1:nz+2)
    w(nx+2,-1:nz+2)=w(2   ,-1:nz+2)
!
! Schema de Peaceman & Rashford, integration selon Z
    da=0. ; db=0. ; dc=0. ; sm=0.
!
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)              &
!$OMP PRIVATE (da,db,dc,sm,ix) &
!$OMP SHARED(Tt,Tb,dtt,kdttdx2,kdttdz2,dttd2x,dttd2z,vx,vz,tp,w,rhs)
    DO ix=1,nx+1 
       da(1:nz)=    -     kdttdz2               -vz(ix,1:nz) *dttd2z
       dc(1:nz)=    -     kdttdz2+ vz(ix,2:nz+1)             *dttd2z
       db(1:nz)=1.+2.*kdttdz2+(vz(ix,2:nz+1)-vz(ix,1:nz))*dttd2z
       sm(1:nz)=w(ix,1:nz)+(w(ix+1,1:nz)-2.*w(ix,1:nz)+w(ix-1,1:nz))*kdttdx2                         &
            -(vx(ix,1:nz)*(w(ix+1,1:nz)+w(ix,1:nz))-vx(ix-1,1:nz)*(w(ix-1,1:nz)+w(ix,1:nz)))*dttd2x &
            +rhs(ix,1:nz)*dtt
!       db( 1)=db( 1)-da( 1) ; sm( 1)=sm( 1)+da( 1) ; da( 1)=0.
!       db(nz)=db(nz)-dc(nz) ; sm(nz)=sm(nz)-dc(nz) ; dc(nz)=0.
       db( 1)=db( 1)-da( 1) ; sm( 1)=sm( 1)-2.*Tt*da( 1) ; da( 1)=0.
       db(nz)=db(nz)-dc(nz) ; sm(nz)=sm(nz)-2.*Tb*dc(nz) ; dc(nz)=0.
       CALL tridiag(da,db,dc,sm,nz)
       tp(ix,1:nz)=sm(1:nz)
    ENDDO
!$OMP END PARALLEL DO
    tp(1:nx+1,  -1)=2.*Tt-tp(1:nx+1,   2)
    tp(1:nx+1,   0)=2.*Tt-tp(1:nx+1,   1)
    tp(1:nx+1,nz+1)=2.*Tb-tp(1:nx+1,nz  )
    tp(1:nx+1,nz+2)=2.*Tb-tp(1:nx+1,nz-1)
    tp(-1  ,-1:nz+2)=tp(nx-1,-1:nz+2)
    tp( 0  ,-1:nz+2)=tp(nx  ,-1:nz+2)
    tp(nx+1,-1:nz+2)=tp(1   ,-1:nz+2)
    tp(nx+2,-1:nz+2)=tp(2   ,-1:nz+2)
!      
  END SUBROUTINE ADI_TP
!
!
  SUBROUTINE ADI_VX(ux,vx,vz,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx,dz
    USE inv_mat_band, ONLY : cyclic,tridiag
    USE mod_bounds  , ONLY : bounds_VX
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,dtt,kdttdx2,kdttdz2,dttd2x,dttd2z
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: vx,vz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: ux
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    dtt=dt/2.
    kdttdx2=diff*dtt/dx**2 ; kdttdz2=diff*dtt/dz**2  
    dttd2x=dtt/dx/2. ; dttd2z=dtt/dz/2.
!
! Schema de Peaceman & Rashford, integration selon X
    w=0.
    da=0. ; db=0. ; dc=0. ; sm=0.
!
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)              &
!$OMP PRIVATE (da,db,dc,sm,iz) &
!$OMP SHARED(dtt,kdttdx2,kdttdz2,dttd2x,dttd2z,vx,vz,ux,w,rhs)
    DO iz=1,nz
       da(1:nx)=    -     kdttdx2               -vx(1:nx,iz) *dttd2x
       dc(1:nx)=    -     kdttdx2+ vx(2:nx+1,iz)           *dttd2x
       db(1:nx)=1.+2.*kdttdx2+(vx(2:nx+1,iz)-vx(1:nx,iz))*dttd2x
       sm(1:nx)=ux(1:nx,iz)+(ux(1:nx,iz+1)-2.*ux(1:nx,iz)+ux(1:nx,iz-1))*kdttdz2                         &
            -(vz(1:nx,iz+1)*(ux(1:nx,iz+1)+ux(1:nx,iz))-vz(1:nx,iz)*(ux(1:nx,iz)+ux(1:nx,iz-1)))*dttd2z &
            +rhs(1:nx,iz)*dtt
       CALL cyclic(da,db,dc,dc(nx),da(1),sm,nx)
       w(1:nx+1,iz)=sm(1:nx+1)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VX(w)
!
! Schema de Peaceman & Rashford, integration selon Z
    da=0. ; db=0. ; dc=0. ; sm=0.
!
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)              &
!$OMP PRIVATE (da,db,dc,sm,ix) &
!$OMP SHARED(dtt,kdttdx2,kdttdz2,dttd2x,dttd2z,vx,vz,ux,w,rhs)
    DO ix=1,nx 
       da(1:nz)=    -     kdttdz2               -vz(ix,1:nz) *dttd2z
       dc(1:nz)=    -     kdttdz2+ vz(ix,2:nz+1)             *dttd2z
       db(1:nz)=1.+2.*kdttdz2+(vz(ix,2:nz+1)-vz(ix,1:nz))*dttd2z
       sm(1:nz)=w(ix,1:nz)+(w(ix+1,1:nz)-2.*w(ix,1:nz)+w(ix-1,1:nz))*kdttdx2                         &
            -(vx(ix+1,1:nz)*(w(ix+1,1:nz)+w(ix,1:nz))-vx(ix,1:nz)*(w(ix-1,1:nz)+w(ix,1:nz)))*dttd2x &
            +rhs(ix,1:nz)*dtt
       db( 1)=db( 1)+da( 1) ; da( 1)=0.
       db(nz)=db(nz)+dc(nz) ; dc(nz)=0.
       CALL tridiag(da,db,dc,sm,nz)
       ux(ix,1:nz)=sm(1:nz)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VX(ux)
!
  END SUBROUTINE ADI_VX
!
!
  SUBROUTINE ADI_VZ(uz,vx,vz,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx,dz
    USE inv_mat_band, ONLY : cyclic,tridiag
    USE mod_bounds  , ONLY : bounds_VZ
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,dtt,kdttdx2,kdttdz2,dttd2x,dttd2z
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: vx,vz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: uz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    dtt=dt/2.
    kdttdx2=diff*dtt/dx**2 ; kdttdz2=diff*dtt/dz**2  
    dttd2x=dtt/dx/2. ; dttd2z=dtt/dz/2.
!
! Schema de Peaceman & Rashford, integration selon X
    w=0.
    da=0. ; db=0. ; dc=0. ; sm=0.
!
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)              &
!$OMP PRIVATE (da,db,dc,sm,iz) &
!$OMP SHARED(dtt,kdttdx2,kdttdz2,dttd2x,dttd2z,vx,vz,uz,w,rhs)
    DO iz=2,nz
       da(1:nx)=    -     kdttdx2             -vx(0:nx-1,iz) *dttd2x
       dc(1:nx)=    -     kdttdx2+ vx(1:nx,iz)               *dttd2x
       db(1:nx)=1.+2.*kdttdx2+(vx(1:nx,iz)-vx(0:nx-1,iz))*dttd2x
       sm(1:nx)=uz(1:nx,iz)+(uz(1:nx,iz+1)-2.*uz(1:nx,iz)+uz(1:nx,iz-1))*kdttdz2                      &
            -(vz(1:nx,iz)*(uz(1:nx,iz+1)+uz(1:nx,iz))-vz(1:nx,iz-1)*(uz(1:nx,iz)+uz(1:nx,iz-1)))*dttd2z &
            +rhs(1:nx,iz)*dtt
       CALL cyclic(da,db,dc,dc(nx),da(1),sm,nx)
       w(1:nx,iz)=sm(1:nx)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VZ(w)
!
! Schema de Peaceman & Rashford, integration selon Z
    da=0. ; db=0. ; dc=0. ; sm=0.
!
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)              &
!$OMP PRIVATE (da,db,dc,sm,ix) &
!$OMP SHARED(dtt,kdttdx2,kdttdz2,dttd2x,dttd2z,vx,vz,uz,w,rhs)
    DO ix=1,nx+1
       da(2:nz)=    -     kdttdz2               -vz(ix,2:nz) *dttd2z
       dc(2:nz)=    -     kdttdz2+ vz(ix,3:nz+1)             *dttd2z
       db(2:nz)=1.+2.*kdttdz2+(vz(ix,3:nz+1)-vz(ix,2:nz))*dttd2z
       sm(2:nz)=w(ix,2:nz)+(w(ix+1,2:nz)-2.*w(ix,2:nz)+w(ix-1,2:nz))*kdttdx2                      &
            -(vx(ix,2:nz)*(w(ix+1,2:nz)+w(ix,2:nz))-vx(ix-1,2:nz)*(w(ix-1,2:nz)+w(ix,2:nz)))*dttd2x &
            +rhs(ix,2:nz)*dtt
       da( 2)=0.
       dc(nz)=0.
       CALL tridiag(da(2),db(2),dc(2),sm(2),nz-1)
       uz(ix,2:nz)=sm(2:nz)
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VZ(uz)
!
  END SUBROUTINE ADI_VZ
!
!
  SUBROUTINE ADI_VXd(ux,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx,dz
    USE inv_mat_band, ONLY : tridiag,cyclic,cyclic_nr,tridiag_nr
    USE mod_bounds  , ONLY : bounds_VX
    USE tableaux, ONLY : Zvx
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,kdttdx2,kdttdz2
    REAL(KIND=8) :: lpp0, lpp1, lpp2
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm,sol
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: ux
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    kdttdx2=diff*dt/2./dx**2
    kdttdz2=diff*dt/2./dz**2
    da=0. ; db=0. ; dc=0. ; sm=0.
    DO iz=1,nz
       da(1:nx)=    -kdttdx2
       db(1:nx)=1.+2.*kdttdx2
       dc(1:nx)=da(1:nx)
    DO ix=1, nx
        lpp0=2/((Zvx(iz-1)-Zvx(iz))*(Zvx(iz-1)-Zvx(iz+1)))
        lpp1=2/((Zvx(iz)-Zvx(iz-1))*(Zvx(iz)-Zvx(iz+1)))
        lpp2=2/((Zvx(iz+1)-Zvx(iz-1))*(Zvx(iz+1)-Zvx(iz)))

        sm(ix)= lpp0*ux(ix, iz-1)*diff*dt/2. &
                +(1.+diff*dt/2.*lpp1)*ux(ix, iz) &
                +lpp2*ux(ix, iz+1)*diff*dt/2. &
                +rhs(ix,iz)*dt/2.
    ENDDO
!       sm(1:nx)=ux(1:nx,iz)+(ux(1:nx,iz+1)-2.*ux(1:nx,iz)+ux(1:nx,iz-1))*kdttdz2 &
!               +rhs(1:nx,iz)*dt/2.
       CALL cyclic_nr(da,db,dc,dc(nx),da(1),sm,sol,nx)
       w(1:nx,iz)=sol(1:nx)
    ENDDO
    CALL bounds_VX(w)
!
    da=0. ; db=0. ; dc=0. ; sm=0.
    DO ix=1,nx
!       da(1:nz)=    -kdttdz2
!       db(1:nz)=1.+2.*kdttdz2
!       dc(1:nz)=da(1:nz)
        DO iz=1, nz
            lpp0=2/((Zvx(iz-1)-Zvx(iz))*(Zvx(iz-1)-Zvx(iz+1)))
            lpp1=2/((Zvx(iz)-Zvx(iz-1))*(Zvx(iz)-Zvx(iz+1)))
            lpp2=2/((Zvx(iz+1)-Zvx(iz-1))*(Zvx(iz+1)-Zvx(iz)))

            da(iz)=   -diff*dt/2.*lpp0
            db(iz)= 1.-diff*dt/2.*lpp1
            dc(iz)=   -diff*dt/2.*lpp2
        ENDDO
       sm(1:nz)=w(ix,1:nz)+(w(ix+1,1:nz)-2.*w(ix,1:nz)+w(ix-1,1:nz))*kdttdx2 &
                 +rhs(ix,1:nz)*dt/2.
       db(1)=db(1)+da(1) ; da(1)=0.
       db(nz)=db(nz)+dc(nz) ; dc(nz)=0.
       CALL tridiag_nr(da,db,dc,sm,sol,nz)
       ux(ix,1:nz)=sol(1:nz)
    ENDDO
    CALL bounds_VX(ux)
!
  END SUBROUTINE ADI_VXd
!
!
  SUBROUTINE ADI_VZd(uz,rhs,diff,dt)
!=====================================================================
    USE constantes  , ONLY : nx,nz,md,dx
    USE inv_mat_band, ONLY : tridiag,cyclic,cyclic_nr,tridiag_nr
    USE mod_bounds  , ONLY : bounds_VZ
    USE tableaux,     ONLY: Zvz
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: diff,dt,kdttdx2
    REAL(KIND=8) :: lpp0, lpp1, lpp2
    REAL(KIND=8), DIMENSION(1:md) :: da,db,dc,sm,sol
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(INOUT) :: uz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN)    :: rhs
!=====================================================================
!
    kdttdx2=diff*dt/2./dx**2
    da=0. ; db=0. ; dc=0. ; sm=0.
    DO iz=1,nz+1
        lpp0=2/((Zvz(iz-1)-Zvz(iz))*(Zvz(iz-1)-Zvz(iz+1)))
        lpp1=2/((Zvz(iz)-Zvz(iz-1))*(Zvz(iz)-Zvz(iz+1)))
        lpp2=2/((Zvz(iz+1)-Zvz(iz-1))*(Zvz(iz+1)-Zvz(iz)))

        da(1:nx)=    -kdttdx2
        db(1:nx)=1.+2.*kdttdx2
        dc(1:nx)=da(1:nx)
        DO ix=1, nx
            sm(ix)=lpp0*uz(ix, iz-1)*diff*dt/2. &
                    +(1.+diff*dt/2.*lpp1)*uz(ix, iz) &
                    +lpp2*uz(ix, iz+1)*diff*dt/2. &
                    +rhs(ix,iz)*dt/2.
        ENDDO
!         sm(1:nx)=uz(1:nx,iz)+(uz(1:nx,iz+1)-2.*uz(1:nx,iz)+uz(1:nx,iz-1))*kdttdz2 &
 !              +rhs(1:nx,iz)*dt/2.
       CALL cyclic_nr(da,db,dc,dc(nx),da(1),sm,sol,nx)
       w(1:nx,iz)=sol(1:nx)
    ENDDO
    CALL bounds_VZ(w)
!
    da=0. ; db=0. ; dc=0. ; sm=0.
    DO ix=1,nx+1
!       da(1:nz+1)=    -kdttdz2
!       db(1:nz+1)=1.+2.*kdttdz2
!       dc(1:nz+1)=da(1:nz+1)
        DO iz=1, nz+1
            lpp0=2/((Zvz(iz-1)-Zvz(iz))*(Zvz(iz-1)-Zvz(iz+1)))
            lpp1=2/((Zvz(iz)-Zvz(iz-1))*(Zvz(iz)-Zvz(iz+1)))
            lpp2=2/((Zvz(iz+1)-Zvz(iz-1))*(Zvz(iz+1)-Zvz(iz)))

            da(iz)=    -diff*dt/2.*lpp0
            db(iz)=  1.-diff*dt/2.*lpp1
            dc(iz)=    -diff*dt/2.*lpp2
        ENDDO
       sm(1:nz+1)=w(ix,1:nz+1)+(w(ix+1,1:nz+1)-2.*w(ix,1:nz+1)+w(ix-1,1:nz+1))*kdttdx2 &
                 +rhs(ix,1:nz+1)*dt/2.
       db(1)=1. ; da(1)=0. ; dc(1)=0. ; sm(1)=0.
       da(2)=0.
       dc(nz)=0.
       db(nz+1)=1. ; da(nz+1)=0. ; dc(nz+1)=0. ; sm(nz+1)=0.
       CALL tridiag_nr(da,db,dc,sm,sol,nz+1)
       uz(ix,1:nz+1)=sol(1:nz+1)
    ENDDO
    CALL bounds_VZ(uz)
!
  END SUBROUTINE ADI_VZd
!
!
END MODULE mod_ADI
