MODULE mod_temper
CONTAINS
!
  SUBROUTINE temper
!=====================================================================
    USE constantes, ONLY : nx,nz
    USE variables,  ONLY : dt
    USE tableaux,   ONLY : vx,vz,tp
    USE mod_adi,    ONLY : ADI_TP
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: rhs
!=====================================================================
!
    rhs=0.
!    CALL ADI_TP(tp,vx,vz,rhs,1.d0,dt)
    CALL ADI_TP(tp,vx,vz,rhs,1.,dt)
!
  END SUBROUTINE temper
!
END MODULE mod_temper
!
!
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
!$  USE OMP_LIB
    IMPLICIT NONE
    INTEGER :: ikx,ix,iz,iz1
    REAL(KIND=8) :: kx,pmoy, lpp0,lpp1,lpp2
    REAL(KIND=8),    DIMENSION(1:nz) :: da,db,dc,smr,smi
    COMPLEX(KIND=8), DIMENSION(1:nz) :: cz
    COMPLEX(KIND=8) :: giz1, giz1p1
!=====================================================================
    XF=0.
    DO ix=1,nx+1
       DO iz=1,nz
!          XF(ix,iz)=+((vx(ix,iz)-vx(ix-1,iz))/dx     &
!                     +(vz(ix,iz+1)-vz(ix,iz))/dz)/Pra/dt   &
!                   -Ray*(tp(ix,iz+1)-tp(ix,iz-1))/2./dz
            XF(ix,iz)=+((vx(ix,iz)-vx(ix-1,iz))/dx     &
                       +(vz(ix,iz+1)-vz(ix,iz))/(Zvz(iz+1)-Zvz(iz)))/Pra/dt

       ENDDO
    ENDDO
!    WRITE(6,*) 'XF avant:',minval(XF(1:nx+1,1:nz)),maxval(XF(1:nx+1,1:nz))
!
! Transformee de Fourier du second membre
    CALL SCFFTM( 1, nx, nz, 1., XF, 2*(nx/2+1), XF, nx/2+1, TABLE, WORK, 0)
!
! Resolution de l'equation de Poisson pour la pression pour chacun des modes
! Mode fondamental
    cz=(0.,0.)
    DO iz=1,nz
!! Methode des rectangles
!!       DO iz1=1,nz
!!          cz(iz)=cz(iz)+cmplx(XF(1,iz1  ),XF(2,iz1  ))*abs(iz-iz1  )*dz*dz/2.
!!       ENDDO
!! Methode des Trapezes

!       cz(iz)=(dcmplx(XF(1,1),XF(2,1))*abs(iz-1)*dz              &
!              +dcmplx(XF(1,1),XF(2,1))*abs(iz-1)*dz)*dz/8.     &
!             +(dcmplx(XF(1,nz),XF(2,nz))*abs(iz-nz  )*dz         &
!              +dcmplx(XF(1,nz),XF(2,nz))*abs(iz-nz  )*dz)*dz/8.
!       DO iz1=1,nz-1
!          cz(iz)=cz(iz)+(dcmplx(XF(1,iz1  ),XF(2,iz1  ))*abs(iz-iz1  )*dz&
!                        +dcmplx(XF(1,iz1+1),XF(2,iz1+1))*abs(iz-iz1-1)*dz)*dz/4.
!       ENDDO

    cz(iz)= cmplx(XF(1,1),XF(2,1))*abs(Zvx(iz)-Zvx(1))/2.*Zvx(1) &
            +cmplx(XF(1,nz),XF(2,nz))*abs(Zvx(iz)-Zvx(nz))/2.*(lz-Zvx(nz))
    DO iz1=1,nz-1
        giz1  =cmplx(XF(1,iz1  ),XF(2,iz1  ))*abs(Zvx(iz)-Zvx(iz1))/2.
        giz1p1=cmplx(XF(1,iz1+1),XF(2,iz1+1))*abs(Zvx(iz)-Zvx(iz1+1))/2.
        cz(iz)=cz(iz)+(Zvx(iz1+1)-Zvx(iz1))*(giz1+giz1p1)/2.
    ENDDO
!! Formule de Newton-Cotes avec n=2.
!!       DO iz1=1,nz-1,2
!!          cz(iz)=cz(iz)+(     cmplx(XF(1,iz1  ),XF(2,iz1  ))*abs(iz-iz1  )*dz/2.  &
!!                        +4.*cmplx(XF(1,iz1+1),XF(2,iz1+1))*abs(iz-iz1-1)*dz/2.  &
!!                        +     cmplx(XF(1,iz1+2),XF(2,iz1+2))*abs(iz-iz1-2)*dz/2.)*dz/3.
!!       ENDDO
!! Formule de Simpson
!!       DO iz1=1,nz-1,2
!!          cz(iz)=cz(iz)+(cmplx(XF(1,iz1  ),XF(2,iz1  ))*abs(iz-iz1  )*dz/2.     /3.  &
!!                        +cmplx(XF(1,iz1+1),XF(2,iz1+1))*abs(iz-iz1-1)*dz/2.*4./3.  &
!!                        +cmplx(XF(1,iz1+2),XF(2,iz1+2))*abs(iz-iz1-2)*dz/2.     /3.)*dz
!!       ENDDO
    ENDDO
    XF(1,1:nz)= real(cz(1:nz))
    XF(2,1:nz)= imag(cz(1:nz))
!    XF(1,1:nz)= 0.
!    XF(2,1:nz)= 0.
!
! Les autres modes
!$OMP PARALLEL DO SCHEDULE(RUNTIME) &
!$OMP DEFAULT(NONE)                 &
!$OMP PRIVATE (ikx,kx,da,db,dc,smr,smi) &
!$OMP SHARED(XF)
!    DO ikx=1,nx/2
!       kx=2.*pi*float(ikx)/lx
!       da(1:nz)= 1.
!       db(1:nz)=-2.-(kx*dz)**2
!       dc(1:nz)= 1.
!       smr(1:nz)=XF(2*ikx+1,1:nz)*dz**2
!       smi(1:nz)=XF(2*ikx+2,1:nz)*dz**2
!       db( 1)=db( 1)+da( 1) ; da( 1)=0.
!       db(nz)=db(nz)+dc(nz) ; dc(nz)=0.
!       CALL tridiag(da,db,dc,smr,nz)
!       CALL tridiag(da,db,dc,smi,nz)
!       XF(2*ikx+1,1:nz)=smr(1:nz)
!       XF(2*ikx+2,1:nz)=smi(1:nz)
!    ENDDO

DO ikx=1,nx/2
   kx=2.*pi*float(ikx)/lx
    DO iz=1, nz
    !rajouter les derivées seconde Lagrange p19
    lpp0=2./((Zvx(iz-1)-Zvx(iz))*(Zvx(iz-1)-Zvx(iz+1)))
    lpp1=2./((Zvx(iz)-Zvx(iz-1))*(Zvx(iz)-Zvx(iz+1)))
    lpp2=2./((Zvx(iz+1)-Zvx(iz-1))*(Zvx(iz+1)-Zvx(iz)))
   da(iz)=   lpp0            !a modidier
   db(iz)=lpp1-kx**2      !a modidier
   dc(iz)= lpp2              !a modidier
   smr(iz)=XF(2*ikx+1,iz)
    smi(iz)=XF(2*ikx+2,iz)
    !faire un write pour comparer da/db/dc nouvaeu et ancien
!    WRITE(17,*) 'da(',iz,')=',da(iz),'daold(',iz,')=',1./dz**2
!    WRITE(18,*) 'db(',iz,')= ',db(iz),'dbold(',iz,')=',-2./dz**2-kx**2
!    WRITE(19,*) 'dc(',iz,')=',dc(iz),'dcold(',iz,')=',1./dz**2
    ENDDO
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
!      pmoy=SUM(pr(1,1:nz))*dx*dz/2.+SUM(pr(nx+1,1:nz))*dx*dz/2.+SUM(pr(2:nx,1:nz))*dx*dz
!      pmoy=pmoy/lx/lz
!      pr=pr-pmoy
        pmoy=0
        DO iz=1, nz
            pmoy=pmoy+pr(1, iz)*dx/2.*(Zvz(iz+1)-Zvz(iz))
            DO ix=2, nx
                pmoy=pmoy+pr(ix, iz)*dx*(Zvz(iz+1)-Zvz(iz))
            ENDDO
            pmoy=pmoy+pr(nx+1, iz)*dx/2.*(Zvz(iz+1)-Zvz(iz))
        ENDDO
        pmoy=pmoy/lx/lz
        !pmoy1=0.
        !pmoy1=(SUM(pr(1,1:nz))*dx*dz/2.+SUM(pr(nx+1,1:nz))*dx*dz/2.+SUM(pr(2:nx,1:nz))*dx*dz)/lx/lz
        !WRITE(13,*) 'pmoy=',pmoy,' pmoy1=', pmoy1

        pr=pr-pmoy
    ENDIF
!
    CALL bounds_PR(pr)
!    WRITE(6,*) 'Pr:',minval(pr(1:nx+1,1:nz)),maxval(pr(1:nx+1,1:nz))
!    pr(-1:nx+3,  -1)=-0.5*Ray*2.*dz+pr(-1:nx+3,   2)
!    pr(-1:nx+3,   0)=-0.5*Ray*1.*dz+pr(-1:nx+3,   1)
!    pr(-1:nx+3,nz+1)=-0.5*Ray*1.*dz+pr(-1:nx+3,  nz)
!    pr(-1:nx+3,nz+2)=-0.5*Ray*2.*dz+pr(-1:nx+3,nz-1)
!    pr(-1,-1:nz+2)=pr(nx-1,-1:nz+2)
!    pr( 0,-1:nz+2)=pr(nx  ,-1:nz+2)
!    pr(nx+1,-1:nz+2)=pr(1,-1:nz+2)
!    pr(nx+2,-1:nz+2)=pr(2,-1:nz+2)
!    pr(nx+3,-1:nz+2)=pr(3,-1:nz+2)
!
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
    USE tableaux,     ONLY : vx,vz,pr
    USE tableaux,     ONLY : Zvx,Zvz
    USE mod_pressure, ONLY : Pr_FFT
    USE mod_adi,      ONLY : ADI_VXd,ADI_VZd
    USE mod_bounds,   ONLY : bounds_VX,bounds_VZ
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: rhs
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: vx0,vz0
    REAL(KIND=8)::lp0, lp1, lp2, alpha, beta, gamma
!=====================================================================
    vx0=vx
    vz0=vz
!=====================================================================
! Resolution des termes non-lineaires - Composante X
!=====================================================================
DO ix=1,nx
    DO iz=1, nz
        lp0=(Zvx(iz)-Zvx(iz+1))/((Zvx(iz-1)-Zvx(iz))*(Zvx(iz-1)-Zvx(iz)+1))
        lp1=(Zvx(iz)-Zvx(iz-1)+Zvx(iz)-Zvx(iz+1))/((Zvx(iz)-Zvx(iz-1))*(Zvx(iz)-Zvx(iz+1)))
        lp2=(Zvx(iz)-Zvx(iz-1))/((Zvx(iz+1)-Zvx(iz-1))*(Zvx(iz+1)-Zvx(iz)))
        vx(ix, iz)=vx0(ix, iz)-((vx0(ix+1,iz)+vx0(ix,iz))**2-(vx0(ix ,iz)+vx0(ix-1,iz))**2)*dt/dx/4. &
                   -dt*((vx0(ix, iz+1)*(Zvz(iz+1)-Zvx(iz))+vx0(ix, iz)*(Zvx(iz+1)-Zvz(iz+1)))/(Zvx(iz+1)-Zvx(iz)) &
                        *(vz0(ix, iz+1)+vz0(ix+1, iz+1))/2. &
                       -(vx0(ix, iz-1)*(Zvx(iz)-Zvz(iz))+vx0(ix, iz)*(Zvz(iz)-Zvx(iz-1)))/(Zvx(iz)-Zvx(iz-1)) &
                        *(vz0(ix, iz)+vz0(ix+1, iz))/2.)  &
                        /(Zvz(iz+1)-Zvz(iz))

    ENDDO
ENDDO
    CALL bounds_VX(vx)
!=====================================================================
! Resolution des termes non-lineaires - Composante Z
!=====================================================================

DO ix=1,nx+1
    DO iz=1, nz+1
        lp0=(Zvz(iz)-Zvz(iz+1))/((Zvz(iz-1)-Zvz(iz))*(Zvz(iz-1)-Zvz(iz+1)))
        lp1=(Zvz(iz)-Zvz(iz-1)+Zvz(iz)-Zvz(iz+1))/((Zvz(iz)-Zvz(iz-1))*(Zvz(iz)-Zvz(iz+1)))
        lp2=(Zvz(iz)-Zvz(iz-1))/((Zvz(iz+1)-Zvz(iz-1))*(Zvz(iz+1)-Zvz(iz)))

        vz(ix, iz)= vz0(ix, iz)-dt/(2.*dx*(Zvx(iz)-Zvx(iz-1))) &
                                    *(vx0(ix,iz)*(Zvz(iz)-Zvx(iz-1)) &
                                    +vx0(ix, iz-1)*(Zvx(iz)-Zvz(iz))) &
                                    *(vz0(ix, iz)+vz0(ix+1, iz)) & !fin premeiere ligne
                              +dt/(2.*dx*(Zvx(iz)-Zvx(iz-1))) &
                                    *(vx0(ix-1,iz)*(Zvz(iz)-Zvx(iz-1)) &
                                    +vx0(ix-1, iz-1)*(Zvx(iz)-Zvz(iz))) &
                                    *(vz0(ix-1, iz)+vz0(ix, iz)) & !fin seconde ligne
                              -dt*(lp0*vz0(ix, iz-1)**2+lp1*vz0(ix, iz)**2+lp2*vz0(ix, iz+1)**2) ! - remplacer par +
    ENDDO
ENDDO

    CALL bounds_VZ(vz)
!=====================================================================
! Rajout du gradient de pression
!=====================================================================
! Calcul de la Pression
    CALL Pr_FFT
!    pr=pr_a
    vx(1:nx,1:nz)=vx(1:nx,1:nz)+(-(pr(2:nx+1,1:nz)-pr(1:nx,1:nz))/dx)*Pra*dt
    CALL bounds_VX(vx)
    DO ix=1, nx+1
        DO iz=1, nz+1
            alpha=((Zvx(iz)-Zvz(iz))**2-(Zvz(iz)-Zvx(iz-1))**2)/2. !alpha = (d+^2 - d-^2)/2
            beta=((Zvx(iz+1)-Zvz(iz))**2-(Zvz(iz)-Zvx(iz-2))**2)/2. !beta = (d++^2 - d--^2)/2
            gamma=alpha*((Zvx(iz+1)-Zvz(iz))+(Zvz(iz)-Zvx(iz-2))) - beta*((Zvx(iz)-Zvz(iz))+(Zvz(iz)-Zvx(iz-1)))
            IF (abs(gamma)>eps) THEN
                vz(ix, iz)=vz(ix, iz)-Pra*dt*(alpha/gamma*pr(ix, iz+1) &
                                             -beta/gamma*pr(ix,iz)     &
                                             +beta/gamma*pr(ix,iz-1)   &
                                             -alpha/gamma*pr(ix,iz-2))
            ELSE
                vz(ix, iz)= vz(ix, iz)-Pra*dt*(pr(ix,iz)-pr(ix,iz-1))/(Zvx(iz)-Zvz(iz)+Zvz(iz)-Zvx(iz-1))
!                   vz(ix,iz)=vz(ix,iz)-(pr(ix,iz)-pr(ix,iz-1))/dz*Pra*dt
            ENDIF
        ENDDO
    ENDDO
!    vz(1:nx+1,1:nz+1)=vz(1:nx+1,1:nz+1)+(-(pr(1:nx+1,1:nz+1)-pr(1:nx+1,0:nz))/dz &
!                                        -Ray*((tp(1:nx+1,1:nz+1)+tp(1:nx+1,0:nz))/2.-T0))*Pra*dt
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
!
!
  SUBROUTINE solKA1_vel
!=====================================================================
! Resolution de Navier Stokes
! High Order Splitting Method
! Karniadakis et al., J. Comp. Phys., 97, 414-443, 1991.
!=====================================================================
    USE constantes,   ONLY : nx,nz,dx,dz,Ray,Pra
    USE variables,    ONLY : dt
    USE tableaux,     ONLY : vx,vx0,vz,vz0,pr,tp
    USE mod_pressure, ONLY : Pr_FFT
    USE mod_adi,      ONLY : ADI_VX,ADI_VZ
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), PARAMETER :: zero=0.
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: rhs
!=====================================================================
    vx0=vx
    vz0=vz
!=====================================================================
! Resolution des termes non-lineaires - Composante X
!=====================================================================
    vx(1:nx,1:nz)=vx0(1:nx,1:nz)                                                                &
         -((vx0(2:nx+1,1:nz)+vx0(1:nx  ,1:nz))**2                                               &
         -(vx0(1:nx  ,1:nz)+vx0(0:nx-1,1:nz))**2)*dt/dx/4.                                    &
         -((vz0(2:nx+1,2:nz+1)+vz0(1:nx,2:nz+1))*(vx0(1:nx,2:nz+1)+vx0(1:nx,1:nz  ))            &
         -(vz0(2:nx+1,1:nz  )+vz0(1:nx,1:nz  ))*(vx0(1:nx,1:nz  )+vx0(1:nx,0:nz-1)))*dt/dz/4.
! Conditions aux limites
    vx(1:nx+1,  -1)=vx(1:nx+1,   2)
    vx(1:nx+1,   0)=vx(1:nx+1,   1)
    vx(1:nx+1,nz+1)=vx(1:nx+1,nz  )
    vx(1:nx+1,nz+2)=vx(1:nx+1,nz-1)
    vx(-1  ,-1:nz+2)=vx(nx-1,-1:nz+2)
    vx( 0  ,-1:nz+2)=vx(nx  ,-1:nz+2)
    vx(nx+1,-1:nz+2)=vx(1   ,-1:nz+2)
    vx(nx+2,-1:nz+2)=vx(2   ,-1:nz+2)
    vx(nx+3,-1:nz+2)=vx(3   ,-1:nz+2)
!=====================================================================
! Resolution des termes non-lineaires - Composante Z
!=====================================================================
    vz(1:nx+1,1:nz+1)=vz0(1:nx+1,1:nz+1)                                                            &
         -((vx0(1:nx+1,1:nz+1)+vx0(1:nx+1,0:nz))*(vz0(2:nx+2,1:nz+1)+vz0(1:nx+1,1:nz+1))            &
         -(vx0(0:nx  ,1:nz+1)+vx0(0:nx  ,0:nz))*(vz0(1:nx+1,1:nz+1)+vz0(0:nx  ,1:nz+1)))*dt/dx/4. &
         -((vz0(1:nx+1,2:nz+2)+vz0(1:nx+1,1:nz+1))**2                                               &
         -(vz0(1:nx+1,1:nz+1)+vz0(1:nx+1,0:nz  ))**2)*dt/dz/4.
! Conditions aux limites
    vz(1:nx,-1)=-vz(1:nx,3)
    vz(1:nx, 0)=-vz(1:nx,2)
    vz(1:nx, 1)=0. 
    vz(1:nx,nz+1)=0.
    vz(1:nx,nz+2)=-vz(1:nx,nz  )
    vz(1:nx,nz+3)=-vz(1:nx,nz-1)
    vz(-1  ,-1:nz+3)=vz(nx-1,-1:nz+3)
    vz( 0  ,-1:nz+3)=vz(nx  ,-1:nz+3)
    vz(nx+1,-1:nz+3)=vz(   1,-1:nz+3)
    vz(nx+2,-1:nz+3)=vz(   2,-1:nz+3)
!=====================================================================
! Rajout du gradient de pression
!=====================================================================
    CALL Pr_FFT
    vx(1:nx,1:nz)=vx(1:nx,1:nz)+(-(pr(2:nx+1,1:nz)-pr(1:nx,1:nz))/dx)*Pra*dt
    vx(1:nx,-1  )=vx(1:nx,2   )
    vx(1:nx, 0  )=vx(1:nx,1   )
    vx(1:nx,nz+1)=vx(1:nx,nz  )
    vx(1:nx,nz+2)=vx(1:nx,nz-1)
    vx(-1  ,-1:nz+2)=vx(nx-1,-1:nz+2)
    vx( 0  ,-1:nz+2)=vx(nx  ,-1:nz+2)
    vx(nx+1,-1:nz+2)=vx(1   ,-1:nz+2)
    vx(nx+2,-1:nz+2)=vx(2   ,-1:nz+2)
    vz(1:nx+1,1:nz+1)=vz(1:nx+1,1:nz+1)+(-(pr(1:nx+1,1:nz+1)-pr(1:nx+1,0:nz))/dz&
         -Ray*(tp(1:nx+1,1:nz+1)+tp(1:nx+1,0:nz))/2.)*Pra*dt
    vz(1:nx+1,-1)=-vz(1:nx+1,3)
    vz(1:nx+1, 0)=-vz(1:nx+1,2)
    vz(1:nx+1, 1)=0.
    vz(1:nx+1,nz+1)=0.
    vz(1:nx+1,nz+2)=-vz(1:nx+1,nz  )
    vz(1:nx+1,nz+3)=-vz(1:nx+1,nz-1)
    vz(-1  ,-1:nz+3)=vz(nx-1,-1:nz+3)
    vz( 0  ,-1:nz+3)=vz(nx  ,-1:nz+3)
    vz(nx+1,-1:nz+3)=vz(1   ,-1:nz+3)
    vz(nx+2,-1:nz+3)=vz(2   ,-1:nz+3)
!=====================================================================
! Resolution des termes lineaires - Composante X
!=====================================================================
    rhs=vx-vx0
    vx=vx0
    CALL ADI_VX(vx,zero,zero,rhs,Pra,dt)
    vx(1:nx+1,  -1)=vx(1:nx+1,   2)
    vx(1:nx+1,   0)=vx(1:nx+1,   1)
    vx(1:nx+1,nz+1)=vx(1:nx+1,nz  )
    vx(1:nx+1,nz+2)=vx(1:nx+1,nz-1)
    vx(-1  ,-1:nz+2)=vx(nx-1,-1:nz+2)
    vx( 0  ,-1:nz+2)=vx(nx  ,-1:nz+2)
    vx(nx+1,-1:nz+2)=vx(1   ,-1:nz+2)
    vx(nx+2,-1:nz+2)=vx(2   ,-1:nz+2)
    vx(nx+3,-1:nz+2)=vx(3   ,-1:nz+2)
!=====================================================================
! Resolution des termes lineaires - Composante Z
!=====================================================================
    rhs=vz-vz0
    vz=vz0
    CALL ADI_VZ(vz,zero,zero,rhs,Pra,dt)
    vz(1:nx,-1)=-vz(1:nx,3)
    vz(1:nx, 0)=-vz(1:nx,2)
    vz(1:nx, 1)=0. 
    vz(1:nx,nz+1)=0.
    vz(1:nx,nz+2)=-vz(1:nx,nz  )
    vz(1:nx,nz+3)=-vz(1:nx,nz-1)
    vz(-1  ,-1:nz+3)=vz(nx-1,-1:nz+3)
    vz( 0  ,-1:nz+3)=vz(nx  ,-1:nz+3)
    vz(nx+1,-1:nz+3)=vz(   1,-1:nz+3)
    vz(nx+2,-1:nz+3)=vz(   2,-1:nz+3)
!
!
  END SUBROUTINE solKA1_vel
!
END MODULE mod_velocity
