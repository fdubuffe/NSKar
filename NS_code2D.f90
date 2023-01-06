PROGRAM Navier_Stokes
!=====================================================================
  USE constantes,   ONLY : nx,nz,dx,dz,new,eps
  USE variables,    ONLY : it,itmax,itold,dt,temps
  USE tableaux,     ONLY : vx,vx_a,vz,vz_a,pr,pr_a,tp
  USE mod_velocity
  USE mod_initio,   ONLY : init,ecriture,bilan,ms_error,sol_analytique
  USE mod_temper,   ONLY : temper
  USE mod_figure,   ONLY : dessin
  IMPLICIT NONE
  INTEGER :: i,j
  REAL(KIND=8) :: err_vx=0.,err_vz=0.,err_pr=0.
!=====================================================================
!
  OPEN(UNIT=95,FILE='gmt_exe',FORM='FORMATTED')
!
  CALL init
!
  IF     (new == 0) THEN
     OPEN(UNIT=10,FILE='Results.dat',FORM='FORMATTED',STATUS='REPLACE')
  ELSEIF ((new == 2) .OR. (new == 31) .OR. (new == 4)) THEN
     itold=0
     it=0
     OPEN(UNIT=10,FILE='Results.dat',FORM='FORMATTED',STATUS='REPLACE')
  ELSE
     OPEN(UNIT=10,FILE='Results.dat',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
  ENDIF
!
  DO WHILE (it < itold+itmax)
     it=it+1
     temps=temps+dt
     IF (new == 4) THEN
        CALL sol_analytique(temps)
     ELSE
        CALL temper
     ENDIF
     CALL solKarn
     IF (new == 4) THEN
        err_vx=err_vx+ms_error(vx,vx_a,1)*dt
        err_vz=err_vz+ms_error(vz,vz_a,2)*dt
        err_pr=err_pr+ms_error(pr,pr_a,3)*dt
     ENDIF
     CALL bilan
     IF (mod(float(it),float(itmax)/10.) == 0.) THEN
!     IF (mod(float(it),100.) == 0.) THEN
       CALL dessin
       CALL ecriture
!       OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
!       WRITE(20) it,temps,vx,vz,pr,tp
!       CLOSE(20)
     ENDIF
  ENDDO
  CALL dessin
  CALL ecriture
  IF (new == 4) THEN
     OPEN(UNIT=11,FILE='RMS_error.dat',FORM='FORMATTED',STATUS='REPLACE')
     WRITE(11,'(4e16.8)') dx,sqrt(err_vx/(temps+eps)),sqrt(err_vz/(temps+eps)),sqrt(err_pr/(temps+eps))
     CLOSE(11)
  ENDIF
!  OPEN(UNIT=20,FILE='DATA_UNF',FORM='UNFORMATTED')
!  WRITE(20) it,temps,vx,vz,pr,tp
!  CLOSE(20)
!
  CLOSE(10)
  CLOSE(95)
!
    DO i=1,nx+1
  DO j=1,nz+1
       WRITE(99,'(3f16.8)') (float(i)-0.5)*dx,1.-(float(j)-1.)*dz,(tp(i,j)+tp(i-1,j))/2.
    ENDDO
    WRITE(99,*) 
  ENDDO

!
END PROGRAM Navier_Stokes
