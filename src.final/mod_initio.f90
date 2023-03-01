MODULE mod_initio
CONTAINS
!
!
  SUBROUTINE init
!=====================================================================
    USE constantes,   ONLY : nx,nz,lx,lz,dx,eps, &
                             Pra,Ray,new,alf0,pi,dx,dz
    USE variables,    ONLY : it,itmax,itold,dt,temps,cour
    USE tableaux,     ONLY : vx,vx_a,vz,vz_a,pr,pr_a,tp,Zvx,Zvz
    USE mod_pressure, ONLY : Pr_FFT
    USE mod_bounds,   ONLY : bounds_PR,bounds_VX,bounds_VZ
    USE fourier,      ONLY : XF,TABLE,WORK
    IMPLICIT NONE
    INTEGER :: iz
!=====================================================================
!
    WRITE(6,'(2(a,e16.8))') 'lx=',lx,' lz=',lz
    WRITE(6,'(2(a,i4))') 'nx=',nx,' nz=',nz
    WRITE(6,'(2(a,e16.8))') 'dx=',dx,' dz=',dz
    WRITE(6,'(2(a,e16.8))') 'Prandtl=',Pra,' Rayleigh=',Ray

!initialisation Zvz
!    DO iz=0, nz+2
!      Zvz(iz)=float(iz-1)*1./float(nz)

DO iz=0, nz+2
!      Zvz(iz)=float(iz-1)*1./float(nz)
  Zvz(iz+1)=0.5-cos(float(2*iz+1)/2./float(nz+1)*pi)/(cos(float(1)/2./float(nz+1)*pi)-cos(float(2*nz+1)/2./float(nz+1)*pi))
ENDDO
    Zvz(1)=0.
    Zvz(nz+1)=1.
    Zvz(0)=-Zvz(2)
    Zvz(-1)=-Zvz(3)
    Zvz(nz+2)=lz+lz-Zvz(nz)
    Zvz(nz+3)=lz+lz-Zvz(nz-1)

!initialistation Zvx
    DO iz=1, nz
      Zvx(iz)=Zvz(iz)+(Zvz(iz+1)-Zvz(iz))/2.
    ENDDO
    Zvx(0)=-Zvx(1)
    Zvx(-1)=-Zvx(2)
    Zvx(nz+1)= lz+lz-Zvx(nz)
    Zvx(nz+2)= lz+lz-Zvx(nz-1)
    Zvx(nz+3)= lz+lz-Zvx(nz-2)


! Initialisation des FFT
    XF=0.
    CALL SCFFTM(0, nx, nz, 1., XF, 2*(nx/2+1), XF, nx/2+1, TABLE, WORK, 0)
!
    
    IF (new == 4) THEN
       WRITE(6,'(a)') 'Test Navier-Stokes: Vortex'
       itold=0
       temps=0.
       dt=cour*MIN(dx**2,MINVAL(Zvz(2:nz+1)-Zvz(1:nz)))
       itmax=NINT(1./alf0/dt)
       dt=1./alf0/float(itmax)
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
    WRITE(6,'(a,i7,2(a,e16.8))') 'itmax=',itmax,' dt=',dt,' tmax=',dt*float(itmax)
    it=itold
!
    CALL ecriture
!
  END SUBROUTINE init
!
  SUBROUTINE sol_analytique(temps)
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,ampV,kx0,kz0,alf0
    USE tableaux,   ONLY : vx_a,vz_a,pr_a,Zvz,Zvx
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
          z=Zvx(iz)
          vx_a(ix,iz)=-kz0*ampV*sin(kx0*x)*cos(kz0*z)*exp(-alf0*temps)
       ENDDO
    ENDDO
    DO ix=-1,nx+3
       DO iz=-1,nz+3
          x=(float(ix)-1.)*dx
          z=Zvz(iz)
          vz_a(ix,iz)= kx0*ampV*cos(kx0*x)*sin(kz0*z)*exp(-alf0*temps)
       ENDDO
    ENDDO
    DO ix=-1,nx+3
       DO iz=-1,nz+3
          x=(float(ix)-1.)*dx
          z=Zvx(iz)
          pr_a(ix,iz)= (cos(2.*kx0*x)*(ampV*kz0/2.)**2+cos(2.*kz0*z)*(ampV*kx0/2.)**2)*exp(-2.*alf0*temps)
       ENDDO
    ENDDO
!
  END SUBROUTINE sol_analytique
!
!
  FUNCTION ms_error(sol1,sol2,cas)
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,lx,lz
    USE tableaux,   ONLY : Zvx, Zvz
    IMPLICIT NONE
    INTEGER :: iz, ix
    INTEGER, INTENT(IN)                      :: cas
    REAL(KIND=8)                             :: ms_error
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: sol1,sol2
!=====================================================================
!
    IF (cas == 1) THEN
! Vx rms
        ms_error=0.
        DO ix=1, nx
            DO iz=1, nz
                ms_error=ms_error+(sol1(ix, iz)-sol2(ix, iz))**2*dx*(Zvz(iz+1)-Zvz(iz))
            ENDDO
        ENDDO
        
    ELSEIF (cas == 2) THEN
! Vz rms
        ms_error=0.
        DO iz=2, nz
            ms_error=ms_error+(sol1(1, iz)-sol2(1, iz))**2*dx*(Zvx(iz)-Zvx(iz-1))/2.
            DO ix=2, nx
               ms_error=ms_error+(sol1(ix, iz)-sol2(ix, iz))**2*dx*(Zvx(iz)-Zvx(iz-1))
            ENDDO
            ms_error=ms_error+(sol1(nx+1, iz)-sol2(nx+1, iz))**2*dx*(Zvx(iz)-Zvx(iz-1))/2.
        ENDDO
        ms_error=ms_error+SUM((sol1(2:nx,1)-sol2(2:nx, 1))**2)*dx*(Zvx(1)-Zvx(0))/2.
        ms_error=ms_error+SUM((sol1(2:nx, nz+1)-sol2(2:nx, nz+1))**2)*dx*(Zvx(nz+1)-Zvx(nz))/2.
        ms_error=ms_error+(sol1(1,1)-sol2(1,1))**2*dx*(Zvx(1)-Zvx(0))/4.
        ms_error=ms_error+(sol1(nx+1,1)-sol2(nx+1,1))**2*dx*(Zvx(1)-Zvx(0))/4.
        ms_error=ms_error+(sol1(1,nz+1)-sol2(1,nz+1))**2*dx*(Zvx(nz+1)-Zvx(nz))/4.
        ms_error=ms_error+(sol1(nx+1,nz+1)-sol2(nx+1,nz+1))**2*dx*(Zvx(nz+1)-Zvx(nz))/4.
            ELSEIF (cas == 3) THEN
        ms_error=0
        DO iz=1, nz
            ms_error=ms_error+(sol1(1, iz)-sol2(1, iz))**2*dx*(Zvz(iz+1)-Zvz(iz))/2.
            DO ix=2, nx
                ms_error=ms_error+(sol1(ix, iz)-sol2(ix, iz))**2*dx*(Zvz(iz+1)-Zvz(iz))
            ENDDO
            ms_error=ms_error+(sol1(nx+1, iz)-sol2(nx+1, iz))**2*dx*(Zvz(iz+1)-Zvz(iz))/2.
        ENDDO
        
    ELSE
       STOP 'RMS non definie'
    ENDIF
       ms_error=ms_error/lx/lz
    RETURN
!
  END FUNCTION ms_error
!
  FUNCTION moyenne(w,cas)
!=====================================================================
! Moyenne sur les points scalaires
!=====================================================================
    USE constantes, ONLY : nx,nz,dx,lx,lz
    USE tableaux,   ONLY : Zvz, Zvx
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: cas
    INTEGER ::ix, iz
    REAL(KIND=8) :: moyenne
    REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3), INTENT(IN) :: w
!=====================================================================
    IF (cas == 1) THEN
! Moyenne sur les points de Vx
    moyenne=0.
    DO ix=1, nx
        DO iz=1, nz
            moyenne=moyenne+(w(ix, iz))*dx*(Zvz(iz+1)-Zvz(iz))
        ENDDO
    ENDDO
    ELSEIF (cas == 2) THEN
! Moyenne sur les points de Vz
moyenne=0
DO iz=2, nz
    moyenne=moyenne+w(1, iz)*dx*(Zvx(iz)-Zvx(iz-1))/2.
    DO ix=2, nx
       moyenne=moyenne+w(ix, iz)*dx*(Zvx(iz)-Zvx(iz-1))
    ENDDO
    moyenne=moyenne+w(nx+1, iz)*dx*(Zvx(iz)-Zvx(iz-1))/2.
ENDDO
moyenne=moyenne+SUM(w(2:nx,1))*dx*(Zvx(1)-Zvx(0))/2.
moyenne=moyenne+SUM(w(2:nx, nz+1))*dx*(Zvx(nz+1)-Zvx(nz))/2.
moyenne=moyenne+w(1,1)*dx*(Zvx(1)-Zvx(0))/4.
moyenne=moyenne+w(nx+1,1)*dx*(Zvx(1)-Zvx(0))/4.
moyenne=moyenne+w(1,nz+1)*dx*(Zvx(nz+1)-Zvx(nz))/4.
moyenne=moyenne+w(nx+1,nz+1)*dx*(Zvx(nz+1)-Zvx(nz))/4.
    ELSEIF (cas == 3) THEN
! Moyenne sur les points de Pr
moyenne=0
DO iz=1, nz
    moyenne=moyenne+w(1, iz)*dx*(Zvz(iz+1)-Zvz(iz))/2.
    DO ix=2, nx
        moyenne=moyenne+w(ix, iz)*dx*(Zvz(iz+1)-Zvz(iz))
    ENDDO
    moyenne=moyenne+w(nx+1, iz)*dx*(Zvz(iz+1)-Zvz(iz))/2.
ENDDO
    ELSE
       STOP 'RMS non definie'
    ENDIF
    moyenne=moyenne/lx/lz
    RETURN
  END FUNCTION moyenne
!
  SUBROUTINE ecriture
!=====================================================================
    USE constantes, ONLY : nx,nz,dx
    USE variables,  ONLY : it,temps,cour
    USE tableaux,   ONLY : vx,vz,pr,tp,Zvx, Zvz, vx_a, vz_a, pr_a
    IMPLICIT NONE
    INTEGER :: i,j
    CHARACTER(LEN=8) :: num
    REAL (KIND=8) :: vxmin, vxmax, vzmin, vzmax, prmin, prmax
!=====================================================================
    WRITE(num,'(i8.8)') it
    OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
    WRITE(20) it,temps,cour,vx,vz,pr,tp
    CLOSE(20)

OPEN(UNIT=60,FILE='vx'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
vxmin=MINVAL( vx(1:nx,1:nz) )
vxmax=MAXVAL( vx(1:nx,1:nz) )
WRITE(6,*) 'it=',it,' Vxmin=',vxmin,' Vxmax=',vxmax
DO i=1,nx
  DO j=1,nz
    WRITE(60,'(3f16.8)') (float(i)-0.5)*dx,Zvx(j),vx(i,j)
  ENDDO
  WRITE(60,*)
ENDDO
CLOSE(60)

OPEN(UNIT=60,FILE='vx_a'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
vxmin=MINVAL( vx_a(1:nx,1:nz) )
vxmax=MAXVAL( vx_a(1:nx,1:nz) )
WRITE(6,*) 'it=',it,' Vx_a min=',vxmin,' Vx_a max=',vxmax
DO i=1,nx
  DO j=1,nz
    WRITE(60,'(3f16.8)') (float(i)-0.5)*dx,Zvx(j),vx_a(i,j)
  ENDDO
  WRITE(60,*)
ENDDO
CLOSE(60)

OPEN(UNIT=60,FILE='vx-vx_a'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
prmin=MINVAL( vx(1:nx,1:nz)-vx_a(1:nx,1:nz) )
prmax=MAXVAL( vx(1:nx,1:nz)-vx_a(1:nx,1:nz) )
WRITE(6,*) 'it=',it,' Vx-Vx_a min=',prmin,' Vx-Vx_a max=',prmax
DO i=1,nx
DO j=1,nz
WRITE(60,'(3f16.8)')    (float(i)-0.5)*dx,Zvx(j),vx(i,j)-vx_a(i,j)
ENDDO
WRITE(60,*)
ENDDO
CLOSE(60)

OPEN(UNIT=60,FILE='vz'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
vzmin=MINVAL( vz(1:nx+1,1:nz+1) )
vzmax=MAXVAL( vz(1:nx+1,1:nz+1) )
WRITE(6,*) 'it=',it,' Vzmin=',vzmin,' Vzmax=',vzmax
DO i=1,nx+1
  DO j=1,nz+1
    WRITE(60,'(3f16.8)') (float(i)-1.)*dx,Zvz(j),vz(i,j)
  ENDDO
  WRITE(60,*)
ENDDO
CLOSE(60)

OPEN(UNIT=60,FILE='vz_a'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
vzmin=MINVAL( vz_a(1:nx+1,1:nz+1) )
vzmax=MAXVAL( vz_a(1:nx+1,1:nz+1) )
WRITE(6,*) 'it=',it,' Vz_a min=',vzmin,' Vz_a max=',vzmax
DO i=1,nx+1
DO j=1,nz+1
WRITE(60,'(3f16.8)')    (float(i)-1.)*dx,Zvz(j),vz_a(i,j)
ENDDO
WRITE(60,*)
ENDDO
CLOSE(60)

OPEN(UNIT=60,FILE='vz-vz_a'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
prmin=MINVAL( vz(1:nx+1,1:nz+1)-vz_a(1:nx+1,1:nz+1) )
prmax=MAXVAL( vz(1:nx+1,1:nz+1)-vz_a(1:nx+1,1:nz+1) )
WRITE(6,*) 'it=',it,' Vz-Vz_a min=',prmin,' Vz-Vz_a max=',prmax
DO i=1,nx+1
DO j=1,nz+1
WRITE(60,'(3f16.8)')    (float(i)-1.)*dx,Zvz(j),vz(i,j)-vz_a(i,j)
ENDDO
WRITE(60,*)
ENDDO
CLOSE(60)

OPEN(UNIT=60,FILE='pr'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
prmin=MINVAL( pr(1:nx+1,1:nz) )
prmax=MAXVAL( pr(1:nx+1,1:nz) )
WRITE(6,*) 'it=',it,' Prmin=',prmin,' Prmax=',prmax
DO i=1,nx+1
  DO j=1,nz
    WRITE(60,'(3f16.8)') (float(i)-1.)*dx,Zvx(j),pr(i,j)
  ENDDO
  WRITE(60,*)
ENDDO
CLOSE(60)

    OPEN(UNIT=60,FILE='pr_a'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
    prmin=MINVAL( pr_a(1:nx+1,1:nz) )
    prmax=MAXVAL( pr_a(1:nx+1,1:nz) )
    WRITE(6,*) 'it=',it,' Pr_a min=',prmin,' Pr_a max=',prmax
    DO i=1,nx+1
    DO j=1,nz
    WRITE(60,'(3f16.8)')    (float(i)-1.)*dx,Zvx(j),pr_a(i,j)
    ENDDO
    WRITE(60,*)
    ENDDO
    CLOSE(60)

OPEN(UNIT=60,FILE='pr-pr_a'//num//'.plt',FORM='FORMATTED',STATUS='REPLACE')
prmin=MINVAL( pr(1:nx+1,1:nz)-pr_a(1:nx+1,1:nz) )
prmax=MAXVAL( pr(1:nx+1,1:nz)-pr_a(1:nx+1,1:nz) )
WRITE(6,*) 'it=',it,' Pr-Pr_a min=',prmin,' Pr-Pr_a max=',prmax
DO i=1,nx+1
DO j=1,nz
WRITE(60,'(3f16.8)')    (float(i)-1.)*dx,Zvx(j),pr(i,j)-pr_a(i,j)      
ENDDO
WRITE(60,*)
ENDDO
CLOSE(60)

  END SUBROUTINE ecriture
!

END MODULE mod_initio
