MODULE mod_pressure
CONTAINS
!
  SUBROUTINE Pr_FFT
!=====================================================================
! Resolution de l'equation de Poisson pourla pression
! (d2/dx2+d2/dz2)P(x,z)=rhs(x,z)
! Transformee de Fourier selon X
! Difference finies centrees du second ordre selon z
! P(k,z-1)-(2+k**2*dz**2)*P(k,z)+P(k,z+1)=rhs(k,z)*dz**2
! Conditions aux limites: 
!   Periodique en x=0,Lx
!   dp/dz=0 en z=0,Lz 
!=====================================================================
    USE constantes,   ONLY : new,nx,nz,dx,lx,lz,pi,Pra
    USE variables,    ONLY : dt
    USE tableaux,     ONLY : vx,vz,pr
    USE inv_mat_band, ONLY : tridiag
    USE fourier,      ONLY : XF,TABLE,WORK 
    USE mod_bounds,   ONLY : bounds_PR
    USE tableaux,     ONLY : Zvx, Zvz
 !!   USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ikx,iz,iz1
    REAL(KIND=8) :: kx,pmoy
    REAL(KIND=8),    DIMENSION(1:nz) :: da,db,dc,smr,smi
    COMPLEX(KIND=8), DIMENSION(1:nz) :: cz
    COMPLEX(KIND=8) :: giz1, giz1p1
!=====================================================================
    XF=0.

!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(XF,Zvz,dt,vx,vz)
    DO iz=1,nz
        XF(1:nx+1,iz)=+((vx(1:nx+1,iz)-vx(0:nx,iz))/dx     &
                   +(vz(1:nx+1,iz+1)-vz(1:nx+1,iz))/(Zvz(iz+1)-Zvz(iz)))/Pra/dt
    ENDDO
!$OMP END PARALLEL DO
! Transformee de Fourier du second membre
    CALL SCFFTM( 1, nx, nz, 1., XF, 2*(nx/2+1), XF, nx/2+1, TABLE, WORK, 0)

! Resolution de l'equation de Poisson pour la pression pour chacun des modes
! Mode fondamental
    cz=(0.,0.)
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(XF,Zvx,cz)
    DO iz=1,nz
    cz(iz)= cmplx(XF(1,1),XF(2,1))*abs(Zvx(iz)-Zvx(1))/2.*Zvx(1) &
            +cmplx(XF(1,nz),XF(2,nz))*abs(Zvx(iz)-Zvx(nz))/2.*(lz-Zvx(nz))
    DO iz1=1,nz-1
        giz1  =cmplx(XF(1,iz1  ),XF(2,iz1  ))*abs(Zvx(iz)-Zvx(iz1))/2.
        giz1p1=cmplx(XF(1,iz1+1),XF(2,iz1+1))*abs(Zvx(iz)-Zvx(iz1+1))/2.
        cz(iz)=cz(iz)+(Zvx(iz1+1)-Zvx(iz1))*(giz1+giz1p1)/2.
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    XF(1,1:nz)= real(cz(1:nz))
    XF(2,1:nz)= imag(cz(1:nz))

!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(XF,Zvx)
DO ikx=1,nx/2
   kx=2.*pi*float(ikx)/lx

   da(1:nz)= 2./((Zvx(0:nz-1)-Zvx(1:nz))*(Zvx(0:nz-1)-Zvx(2:nz+1))) !a modidier
   db(1:nz)= (2./((Zvx(1:nz)-Zvx(0:nz-1))*(Zvx(1:nz)-Zvx(2:nz+1))))-kx**2 !a modidier
   dc(1:nz)= 2./((Zvx(2:nz+1)-Zvx(0:nz-1))*(Zvx(2:nz+1)-Zvx(1:nz))) !a modidier
   smr(1:nz)=XF(2*ikx+1,1:nz)
    smi(1:nz)=XF(2*ikx+2,1:nz)
    
   db( 1)=db( 1)+da( 1) ; da( 1)=0.
   db(nz)=db(nz)+dc(nz) ; dc(nz)=0.
   CALL tridiag(da,db,dc,smr,nz)
   CALL tridiag(da,db,dc,smi,nz)
   XF(2*ikx+1,1:nz)=smr(1:nz)
   XF(2*ikx+2,1:nz)=smi(1:nz)
ENDDO
!$OMP END PARALLEL DO
!
! Transformee de Fourier inverse
    CALL CSFFTM(-1, nx, nz, 1./float(nx), XF, nx/2+1, XF, 2*(nx/2+1), TABLE, WORK, 0)
! Calcul de la moyenne
    pr(1:nx+1,1:nz)=XF(1:nx+1,1:nz)
    IF (new == 4) THEN

        pmoy=0
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(Zvz,pr,pmoy)
        DO iz=1, nz
            pmoy=pmoy+pr(1, iz)*dx/2.*(Zvz(iz+1)-Zvz(iz))
            pmoy=pmoy+SUM(pr(2:nx, iz)*dx*(Zvz(iz+1)-Zvz(iz)))
            pmoy=pmoy+pr(nx+1, iz)*dx/2.*(Zvz(iz+1)-Zvz(iz))
        ENDDO
!$OMP END PARALLEL DO
        pmoy=pmoy/lx/lz
        pr=pr-pmoy
    ENDIF
!
    CALL bounds_PR(pr)
  END SUBROUTINE Pr_FFT
!
END MODULE mod_pressure
!
!
MODULE mod_velocity
CONTAINS
!
  SUBROUTINE solKarn
!=====================================================================
! Resolution de Navier Stokes
! High Order Splitting Method
! Karniadakis et al., J. Comp. Phys., 97, 414-443, 1991.
!=====================================================================
    USE constantes,   ONLY : nx,nz,dx,Pra,eps
    USE variables,    ONLY : dt
    USE tableaux,     ONLY : vx,vz,pr, vx0, vz0, rhs
    USE tableaux,     ONLY : Zvx,Zvz
    USE mod_pressure, ONLY : Pr_FFT
    USE mod_adi,      ONLY : ADI_VXd,ADI_VZd
    USE mod_bounds,   ONLY : bounds_VX,bounds_VZ
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8)::lp0, lp1, lp2, alpha, beta, gamma
!=====================================================================
    vx0=vx
    vz0=vz
!=====================================================================
! Resolution des termes non-lineaires - Composante X
!=====================================================================
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(Zvz,Zvx,vx0,vx,dt,vz0)
    DO iz=1, nz

        vx(1:nx, iz)=vx0(1:nx, iz)-((vx0(2:nx+1,iz)+vx0(1:nx,iz))**2-(vx0(1:nx ,iz)+vx0(0:nx-1,iz))**2)*dt/dx/4. &
                   -dt*((vx0(1:nx, iz+1)*(Zvz(iz+1)-Zvx(iz))+vx0(1:nx, iz)*(Zvx(iz+1)-Zvz(iz+1)))/(Zvx(iz+1)-Zvx(iz)) &
                        *(vz0(1:nx, iz+1)+vz0(2:nx+1, iz+1))/2. &
                       -(vx0(1:nx, iz-1)*(Zvx(iz)-Zvz(iz))+vx0(1:nx, iz)*(Zvz(iz)-Zvx(iz-1)))/(Zvx(iz)-Zvx(iz-1)) &
                        *(vz0(1:nx, iz)+vz0(2:nx+1, iz))/2.)  &
                        /(Zvz(iz+1)-Zvz(iz))
ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VX(vx)
!=====================================================================
! Resolution des termes non-lineaires - Composante Z
!=====================================================================
!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(Zvz,Zvx,vz0,vz,dt,vx0)
    DO iz=1, nz+1
        lp0=(Zvz(iz)-Zvz(iz+1))/((Zvz(iz-1)-Zvz(iz))*(Zvz(iz-1)-Zvz(iz+1)))
        lp1=(Zvz(iz)-Zvz(iz-1)+Zvz(iz)-Zvz(iz+1))/((Zvz(iz)-Zvz(iz-1))*(Zvz(iz)-Zvz(iz+1)))
        lp2=(Zvz(iz)-Zvz(iz-1))/((Zvz(iz+1)-Zvz(iz-1))*(Zvz(iz+1)-Zvz(iz)))

        vz(1:nx+1, iz)= vz0(1:nx+1, iz)-dt/(2.*dx*(Zvx(iz)-Zvx(iz-1))) &
                                    *(vx0(1:nx+1,iz)*(Zvz(iz)-Zvx(iz-1)) &
                                    +vx0(1:nx+1, iz-1)*(Zvx(iz)-Zvz(iz))) &
                                    *(vz0(1:nx+1, iz)+vz0(2:nx+2, iz)) &
                              +dt/(2.*dx*(Zvx(iz)-Zvx(iz-1))) &
                                    *(vx0(0:nx,iz)*(Zvz(iz)-Zvx(iz-1)) &
                                    +vx0(0:nx, iz-1)*(Zvx(iz)-Zvz(iz))) &
                                    *(vz0(0:nx, iz)+vz0(1:nx+1, iz)) &
                              -dt*(lp0*vz0(1:nx+1, iz-1)**2+lp1*vz0(1:nx+1, iz)**2+lp2*vz0(1:nx+1, iz+1)**2) !
    ENDDO
!$OMP END PARALLEL DO

    CALL bounds_VZ(vz)
!=====================================================================
! Rajout du gradient de pression
!=====================================================================
! Calcul de la Pression
    CALL Pr_FFT
    vx(1:nx,1:nz)=vx(1:nx,1:nz)+(-(pr(2:nx+1,1:nz)-pr(1:nx,1:nz))/dx)*Pra*dt
    CALL bounds_VX(vx)

!$OMP PARALLEL DO DEFAULT(PRIVATE)  SHARED(Zvz,Zvx,pr,dt,vz)
        DO iz=1, nz+1
            alpha=((Zvx(iz)-Zvz(iz))**2-(Zvz(iz)-Zvx(iz-1))**2)/2.
            beta=((Zvx(iz+1)-Zvz(iz))**2-(Zvz(iz)-Zvx(iz-2))**2)/2.
            gamma=alpha*((Zvx(iz+1)-Zvz(iz))+(Zvz(iz)-Zvx(iz-2))) - beta*((Zvx(iz)-Zvz(iz))+(Zvz(iz)-Zvx(iz-1)))
            DO ix=1, nx+1
            IF (abs(gamma)>eps) THEN
                vz(ix, iz)=vz(ix, iz)-Pra*dt*(alpha/gamma*pr(ix, iz+1) &
                                             -beta/gamma*pr(ix,iz)     &
                                             +beta/gamma*pr(ix,iz-1)   &
                                             -alpha/gamma*pr(ix,iz-2))
            ELSE
                vz(ix, iz)= vz(ix, iz)-Pra*dt*(pr(ix,iz)-pr(ix,iz-1))/(Zvx(iz)-Zvz(iz)+Zvz(iz)-Zvx(iz-1))
            ENDIF
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    CALL bounds_VZ(vz)
!=====================================================================
! Rajout des termes lineaires - Composante X
!=====================================================================
    rhs=(vx-vx0)/dt
    vx=vx0
    CALL ADI_VXd(vx,rhs,Pra,dt)
    CALL bounds_VX(vx)
!=====================================================================
! Rajout des termes lineaires - Composante Z
!!=====================================================================
    rhs=(vz-vz0)/dt
    vz=vz0
    CALL ADI_VZd(vz,rhs,Pra,dt)
    CALL bounds_VZ(vz)
!
  END SUBROUTINE solKarn

END MODULE mod_velocity
