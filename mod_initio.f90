MODULE mod_initio
CONTAINS
!
!
  SUBROUTINE init
!=====================================================================
    USE constantes,   ONLY : nx,nz,lx,lz,dx,dz,ampV,ampT,eps, &
                             Pra,Ray,new,kx0,kz0,alf0
    USE variables,    ONLY : it,itmax,itold,dt,temps,cour
    USE tableaux,     ONLY : vx,vx_a,vz,vz_a,pr,pr_a,tp
    USE mod_figure,   ONLY : dessin
    USE mod_pressure, ONLY : Pr_FFT
    USE mod_bounds,   ONLY : bounds_PR,bounds_VX,bounds_VZ,bounds_TP
    USE fourier,      ONLY : XF,TABLE,WORK
    IMPLICIT NONE
    INTEGER :: ix,iz
    REAL(KIND=8) :: x,z,a,b,c,delta,tau1,tau2,ampPsi
    REAL(KIND=8), DIMENSION(-1:nx/2+3,-1:nz/2+3) :: vxu,vzu,pru,tpu
!=====================================================================
!
    WRITE(6,'(2(a,e16.8))') 'lx=',lx,' lz=',lz
    WRITE(6,'(2(a,i4))') 'nx=',nx,' nz=',nz
    WRITE(6,'(2(a,e16.8))') 'dx=',dx,' dz=',dz
    WRITE(6,'(2(a,e16.8))') 'Prandtl=',Pra,' Rayleigh=',Ray
!
! Initialisation des FFT
    XF=0.
    CALL SCFFTM(0, nx, nz, 1., XF, 2*(nx/2+1), XF, nx/2+1, TABLE, WORK, 0)
!
    IF (new == 0 ) THEN
       WRITE(6,'(a)') 'Condition initiale: nouveau cas'
       itold=0
       temps=0.
       a=1.
       b=(Pra+1.)*(kx0**2+kz0**2)
       c=((kx0**2+kz0**2)**3-Ray*kx0**2)*Pra/(kx0**2+kz0**2)
       delta=b**2-4.*a*c
       tau1=(-b-sqrt(delta))/2./a
       tau2=(-b+sqrt(delta))/2./a
       WRITE(6,'(a,e16.8)') 'Taux de croissance:',tau2
       ampPsi=(tau2+kx0**2+kz0**2)*ampV/kx0
       DO iz=-1,nz+2
          DO ix=-1,nx+2
             x=(float(ix)-0.5)*dx
             z=(float(iz)-0.5)*dz
             vx(ix,iz)=-sin(kx0*x)*cos(kz0*z)*ampV/kx0
          ENDDO
       ENDDO
       DO iz=-1,nz+3
          DO ix=-1,nx+3
             x=float(ix-1)*dx
             z=float(iz-1)*dz
             vz(ix,iz)= cos(kx0*x)*sin(kz0*z)*ampV/kz0
          ENDDO
       ENDDO
       DO iz=-1,nz+2
          DO ix=-1,nx+3
             x=float(ix-1)*dx
             z=(float(iz)-0.5)*dz
             pr(ix,iz)= (cos(2.*kx0*x)*(ampV*kz0/2.)**2+cos(2.*kz0*z)*(ampV*kx0/2.)**2)
          ENDDO
       ENDDO
       DO iz=-1,nz+2
          DO ix=-1,nx+3
             x=float(ix-1)*dx
             z=(float(iz)-0.5)*dz
             tp(ix,iz)=(z-0.5)+cos(kx0*x)*sin(kz0*z)*ampT
          ENDDO
       ENDDO
       dt=cour*MIN(dx/(MAXVAL(ABS(vx(1:nx,1:nz)))+eps),dz/(MAXVAL(ABS(vz(1:nx+1,1:nz+1)))+eps),&
                   1./(1./dx**2+1./dz**2))
       CALL Pr_FFT
    ELSEIF (new == 11) THEN
       WRITE(6,'(a)') 'Condition initiale: poursuite'
       OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
       READ(20) itold,temps,vx,vz,pr,tp
       CLOSE(20)
       dt=cour*MIN(dx/(MAXVAL(ABS(vx(1:nx,1:nz)))+eps),dz/(MAXVAL(ABS(vz(1:nx+1,1:nz+1)))+eps),&
                   1./(1./dx**2+1./dz**2))
    ELSEIF (new == 2) THEN
       WRITE(6,'(a)') 'Condition initiale: poursuite avec remise a zero du temps'
       OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
       READ(20) itold,temps,cour,vx,vz,pr,tp
       CLOSE(20)
       dt=cour*MIN(dx/(MAXVAL(ABS(vx(1:nx,1:nz)))+eps),dz/(MAXVAL(ABS(vz(1:nx+1,1:nz+1)))+eps),&
                   1./(1./dx**2+1./dz**2))
       itold=0
       temps=0.
    ELSEIF (new == 3) THEN
       WRITE(6,'(a)') 'Condition initiale: poursuite avec interpolation'
       OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
       READ(20) itold,temps,cour,vxu,vzu,pru,tpu
       CLOSE(20)
       CALL interp(vxu,vx)
       CALL bounds_VX(vx)
       CALL interp(vzu,vz)
       CALL bounds_VZ(vz)
       CALL interp(pru,pr)
       CALL bounds_PR(pr)
       CALL interp(tpu,tp)
       CALL bounds_TP(tp) 
       dt=cour*MIN(dx/(MAXVAL(ABS(vx(1:nx,1:nz)))+eps),dz/(MAXVAL(ABS(vz(1:nx+1,1:nz+1)))+eps),&
                   1./(1./dx**2+1./dz**2))
       CALL Pr_FFT
    ELSEIF (new == 31) THEN
       WRITE(6,'(a)') 'Condition initiale: poursuite avec interpolation'
       OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
       READ(20) itold,temps,cour,vxu,vzu,pru,tpu
       CLOSE(20)
       CALL interp(vxu,vx)
       CALL bounds_VX(vx)
       CALL interp(vzu,vz)
       CALL bounds_VZ(vz)
       CALL interp(pru,pr)
       CALL bounds_PR(pr)
       CALL interp(tpu,tp)
       CALL bounds_TP(tp) 
       dt=cour*MIN(dx/(MAXVAL(ABS(vx(1:nx,1:nz)))+eps),dz/(MAXVAL(ABS(vz(1:nx+1,1:nz+1)))+eps),&
                   1./(1./dx**2+1./dz**2))
       itold=0
       temps=0.
       CALL Pr_FFT
    ELSEIF (new == 4) THEN
       WRITE(6,'(a)') 'Test Navier-Stokes: Vortex'
       itold=0
       temps=0.
       dt=cour*MIN(dx**2,dz**2)
       itmax=NINT(1./alf0/dt)
       dt=1./alf0/DBLE(itmax)
       CALL sol_analytique(temps)
       vx=vx_a
       vz=vz_a
       pr=pr_a
       tp=0.
    ELSE
       WRITE(6,'(a)') 'Condition initiale: poursuite'
       OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
       READ(20) itold,temps,cour,vx,vz,pr,tp
       CLOSE(20)
       dt=cour*MIN(dx/(MAXVAL(ABS(vx(1:nx,1:nz)))+eps),dz/(MAXVAL(ABS(vz(1:nx+1,1:nz+1)))+eps),&
                   1./(1./dx**2+1./dz**2))
    ENDIF
    WRITE(6,'(a,i5,2(a,e16.8))') 'itmax=',itmax,' dt=',dt,' tmax=',dt*float(itmax)
    it=itold
!
    CALL dessin
    
!
  END SUBROUTINE init
!
!
  SUBROUTINE interp(wu,w)
!=====================================================================
    USE constantes, ONLY : nx,nz
    USE mod_bounds, ONLY : bounds_cond
    IMPLICIT NONE
    INTEGER, PARAMETER :: nxu=nx/2,nzu=nz/2
    REAL(KIND=8), DIMENSION(-1:nxu+3,-1:nzu+3), INTENT(IN) :: wu
    REAL(KIND=8), DIMENSION(-1:nx +3,-1:nz +3), INTENT(OUT):: w
!    CHARACTER(LEN=1), DIMENSION(2), INTENT(IN) :: cl
    INTEGER :: ix,iz
!=====================================================================
    DO iz=1,nzu+1
       DO ix=1,nxu+1
         w(2*ix-1,2*iz-1)=wu(ix,iz)
       ENDDO
       DO ix=1,nxu
         w(2*ix,2*iz-1)=(wu(ix,iz)+wu(ix+1,iz))/2. 
       ENDDO
    ENDDO
    DO ix=1,nx+1
       DO iz=2,nz,2
          w(ix,iz)=(w(ix,iz-1)+w(ix,iz+1))/2.
       ENDDO
    ENDDO
!    CALL bounds_cond(w,cl)
!
  END SUBROUTINE interp
!
!
  SUBROUTINE sol_analytique(temps)
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,dz,ampV,kx0,kz0,alf0
    USE tableaux,   ONLY : vx_a,vz_a,pr_a
    IMPLICIT NONE
    INTEGER                  :: ix,iz
    REAL(KIND=8)             :: x,z
    REAL(KIND=8), INTENT(IN) :: temps
!=====================================================================
! conditions initiale sur la vitesse
! On suppose  la solution a t=0 
! psi=   A*sin(kx*x)*sin(kz*z)
! vx=-kz*A*sin(kx*x)*cos(kz*z) ; dvx/dx=-kx*kz*A*cos(kx*x)*cos(kz*z) ; dvx/dz= kz**2*A*sin(kx*x)*sin(kz*z)
! vz=+kx*A*cos(kx*x)*sin(kz*z) ; dvz/dx=-kx**2*A*sin(kx*x)*sin(kz*z) ; dvz/dz= kx*kz*A*cos(kx*x)*cos(kz*z)
! pr=cos(2*kx*x)*(A*kz)**2/4+cos(2*kz*z)*(A*kx)**2/4
! (v.grad)vx=A**2*kx*kz**2*sin(2*kx*x)/2
! (v.grad)vz=A**2*kx**2*kz*sin(2*kz*z)/2
! div((v.grad)v)=(A*kx*kz)**2*(cos(2*kx*x)+cos(2*kz*z))
    DO ix=-1,nx+3
       DO iz=-1,nz+3
          x=(float(ix)-0.5)*dx
          z=(float(iz)-0.5)*dz
          vx_a(ix,iz)=-kz0*ampV*sin(kx0*x)*cos(kz0*z)*exp(-alf0*temps)
       ENDDO
    ENDDO
    DO ix=-1,nx+3
       DO iz=-1,nz+3
          x=(float(ix)-1.)*dx
          z=(float(iz)-1.)*dz
          vz_a(ix,iz)= kx0*ampV*cos(kx0*x)*sin(kz0*z)*exp(-alf0*temps)
       ENDDO
    ENDDO
    DO ix=-1,nx+3
       DO iz=-1,nz+3
          x=(float(ix)-1.)*dx
          z=(float(iz)-0.5)*dz
          pr_a(ix,iz)= (cos(2.*kx0*x)*(ampV*kz0/2.)**2+cos(2.*kz0*z)*(ampV*kx0/2.)**2)*exp(-2.*alf0*temps)
       ENDDO
    ENDDO
!
  END SUBROUTINE sol_analytique
!
!
  FUNCTION ms_error(sol1,sol2,cas)
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,dz,lx,lz
    IMPLICIT NONE
    INTEGER, INTENT(IN)                      :: cas
    REAL(KIND=8)                             :: ms_error,moy1,moy2
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: sol1,sol2
!=====================================================================
!
    IF (cas == 1) THEN
! Vx rms
       moy1=moyenne(sol1,cas)
       moy2=moyenne(sol2,cas)
       ms_error=SUM((sol1(1:nx,1:nz)-sol2(1:nx,1:nz)+moy2-moy1)**2)*dx*dz
    ELSEIF (cas == 2) THEN
! Vz rms
       moy1=moyenne(sol1,cas)
       moy2=moyenne(sol2,cas)
       ms_error=(sol1(   1,   1)-sol2(   1,   1)+moy2-moy1)**2*dx*dz/4. &
               +(sol1(nx+1,   1)-sol2(nx+1,   1)+moy2-moy1)**2*dx*dz/4. &
               +(sol1(   1,nz+1)-sol2(   1,nz+1)+moy2-moy1)**2*dx*dz/4. &
               +(sol1(nx+1,nz+1)-sol2(nx+1,nz+1)+moy2-moy1)**2*dx*dz/4. &
               +SUM((sol1(2:nx,   1)-sol2(2:nx,   1)+moy2-moy1)**2)*dx*dz/2. &
               +SUM((sol1(2:nx,nz+1)-sol2(2:nx,nz+1)+moy2-moy1)**2)*dx*dz/2. &
               +SUM((sol1(   1,2:nz)-sol2(   1,2:nz)+moy2-moy1)**2)*dx*dz/2. &
               +SUM((sol1(nx+1,2:nz)-sol2(nx+1,2:nz)+moy2-moy1)**2)*dx*dz/2.
       ms_error=ms_error+SUM((sol1(2:nx,2:nz)-sol2(2:nx,2:nz)+moy2-moy1)**2)*dx*dz
    ELSEIF (cas == 3) THEN
! Pr rms
       moy1=moyenne(sol1,cas)
       moy2=moyenne(sol2,cas)
       ms_error=SUM((sol1(   1,1:nz)-sol2(   1,1:nz)+moy2-moy1)**2)*dx*dz/2. &
               +SUM((sol1(nx+1,1:nz)-sol2(nx+1,1:nz)+moy2-moy1)**2)*dx*dz/2.
       ms_error=ms_error+SUM((sol1(2:nx,1:nz)-sol2(2:nx,1:nz)+moy2-moy1)**2)*dx*dz 
    ELSE
       STOP 'RMS non definie'
    ENDIF
       ms_error=ms_error/lx/lz
    RETURN
!
  END FUNCTION ms_error
!
!
  SUBROUTINE bilan 
!=====================================================================
    USE constantes, ONLY : new,nx,nz,dx,dz,Pra,eps
    USE variables,  ONLY : it,dt,temps,cour
    USE tableaux,   ONLY : vx,vx_a,vz,vz_a,pr,pr_a,tp
    IMPLICIT NONE
    REAL(KIND=8) :: ecm,dv,tm,pegx,pegz,flt,flb
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
!=====================================================================
!
! Moyenne de l'energie cinetique, RMS de la divergence
    w(1:nx+1,1:nz)=((vx(1:nx+1,1:nz  )+vx(0:nx  ,1:nz))/2.)**2  &
                  +((vz(1:nx+1,2:nz+1)+vz(1:nx+1,1:nz))/2.)**2
    ecm=moyenne(w,1)
    w(1:nx+1,1:nz)=((vx(1:nx+1,1:nz  )-vx(0:nx  ,1:nz))/dx    &
                   +(vz(1:nx+1,2:nz+1)-vz(1:nx+1,1:nz))/dz)**2
    dv=sqrt(moyenne(w,1))/(sqrt(ecm)+eps)
    IF (dv > 1.) THEN
      STOP 'ERROR DV > 1'
    ENDIF
    IF((dv > 1.e-02) .AND. (new /= 4)) THEN
      cour=cour/2.
    ELSEIF ((dv < 1.e-04) .AND. (new /= 4)) THEN
      cour=min(cour*10.,1.e-02)
    ENDIF
    ecm=ecm/2.
!
    IF (new == 4) THEN
       WRITE(10,'(i6,7e16.8)') it,dt,temps,ecm,dv,sqrt(ms_error(vx,vx_a,1)),sqrt(ms_error(vz,vz_a,2)),sqrt(ms_error(pr,pr_a,3))
    ELSE
! Temperature moyenne
       tm=moyenne(tp,1)
! Flux de chaleur
       flt=flux(tp,1)
       flb=flux(tp,nz+1)
! Calcul du nouveau pas de temps
       pegx=MAXVAL(ABS(vx(1:nx  ,1:nz  )))
       pegz=MAXVAL(ABS(vz(1:nx+1,1:nz+1)))
       dt=MIN(dx/(pegx+eps),dz/(pegz+eps),1./(1./dx**2+1./dz**2))
       dt=MIN(dt,1./(1./dx**2+1./dz**2)/Pra)
       dt=cour*dt
       pegx=pegx*dx
       pegz=pegz*dz
! Ecriture
       WRITE(10,'(i10,10e20.12)') it,dt,temps,flt,flb,tm,ecm,dv,pegx,pegz,cour
    ENDIF
!
  END SUBROUTINE bilan
!
! 
  FUNCTION flux(w,iz)
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,dz,lx
    IMPLICIT NONE
    INTEGER :: iz
    REAL(KIND=8) :: flux
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: w
!=====================================================================
    flux=(w(1,iz)-w(1,iz-1))*dx/2./dz+(w(nx+1,iz)-w(nx+1,iz-1))*dx/2./dz
    flux=flux+SUM(w(2:nx,iz)-w(2:nx,iz-1))*dx/dz
    flux=flux/lx
    RETURN
  END FUNCTION flux
!
!
  FUNCTION moyenne(w,cas)
!=====================================================================
! Moyenne sur les points scalaires
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,dz,lx,lz
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: cas
    REAL(KIND=8) :: moyenne
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN) :: w
!=====================================================================
    IF (cas == 1) THEN
! Moyenne sur les points de Vx
       moyenne=SUM(w(1:nx,1:nz))*dx*dz
    ELSEIF (cas == 2) THEN
! Moyenne sur les points de Vz
       moyenne=w(   1,   1)*dx*dz/4. &
              +w(nx+1,   1)*dx*dz/4. &
              +w(   1,nz+1)*dx*dz/4. &
              +w(nx+1,nz+1)*dx*dz/4. &
              +SUM(w(2:nx,   1))*dx*dz/2. &
              +SUM(w(2:nx,nz+1))*dx*dz/2. &
              +SUM(w(   1,2:nz))*dx*dz/2. &
              +SUM(w(nx+1,2:nz))*dx*dz/2.
       moyenne=moyenne+SUM(w(2:nx,2:nz))*dx*dz
    ELSEIF (cas == 3) THEN
! Moyenne sur les points de Pr
       moyenne=SUM(w(   1,1:nz))*dx*dz/2. &
              +SUM(w(nx+1,1:nz))*dx*dz/2.
       moyenne=moyenne+SUM(w(2:nx,1:nz))*dx*dz 
    ELSE
       STOP 'RMS non definie'
    ENDIF
    moyenne=moyenne/lx/lz
    RETURN
  END FUNCTION moyenne
!
  SUBROUTINE ecriture
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,dz
    USE variables,  ONLY : it,temps,cour
    USE tableaux,   ONLY : vx,vz,pr,tp
    IMPLICIT NONE
    INTEGER :: i,j
    CHARACTER(LEN=8) :: num
!=====================================================================
    WRITE(num,'(i8.8)') it
    OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
    WRITE(20) it,temps,cour,vx,vz,pr,tp
    CLOSE(20)
    OPEN(UNIT=60,FILE='tp'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
    DO i=1,nx+1
      DO j=1,nz+1
        WRITE(60,'(3f16.8)') (float(i)-0.5)*dx,1.-(float(j)-1.)*dz,(tp(i,j)+tp(i-1,j))/2.
      ENDDO
      WRITE(60,*) 
    ENDDO
    CLOSE(60)
  END SUBROUTINE ecriture
!

END MODULE mod_initio
