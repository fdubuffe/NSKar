PROGRAM Navier_Stokes
!=====================================================================
  USE constantes,   ONLY : dx,new,eps
  USE variables,    ONLY : it,itmax,itold,dt,temps
  USE tableaux,     ONLY : vx,vx_a,vz,vz_a,pr,pr_a
  USE mod_velocity
  USE mod_initio,   ONLY : init,ecriture,ms_error,sol_analytique
  IMPLICIT NONE
  REAL(KIND=8) :: err_vx=0.,err_vz=0.,err_pr=0.
!=====================================================================
!
  OPEN(UNIT=95,FILE='gmt_exe',FORM='FORMATTED')
!
  CALL init
!
  OPEN(UNIT=10,FILE='Results.dat',FORM='FORMATTED',STATUS='OLD',POSITION='APPEND')
!
  DO WHILE (it < itold+itmax)
     it=it+1
     temps=temps+dt
     IF (new == 4) THEN
        CALL sol_analytique(temps)
     ENDIF
     CALL solKarn
     IF (new == 4) THEN
        err_vx=err_vx+ms_error(vx,vx_a,1)*dt
        err_vz=err_vz+ms_error(vz,vz_a,2)*dt
        err_pr=err_pr+ms_error(pr,pr_a,3)*dt
     ENDIF
     IF (abs(mod(float(it),float(itmax)/10.)) < 1.e-14) THEN
       CALL ecriture
     ENDIF
  ENDDO
  CALL ecriture
  IF (new == 4) THEN
     OPEN(UNIT=11,FILE='RMS_error.dat',FORM='FORMATTED',STATUS='REPLACE')
     WRITE(11,'(4e16.8)') dx,sqrt(err_vx/(temps+eps)),sqrt(err_vz/(temps+eps)),sqrt(err_pr/(temps+eps))
     CLOSE(11)
  ENDIF

  CLOSE(10)
  CLOSE(95)

!
END PROGRAM Navier_Stokes
