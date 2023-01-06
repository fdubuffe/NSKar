MODULE mod_figure
CONTAINS
!
!
  SUBROUTINE dessin
!===============================================================================
!    USE constantes, ONLY : nx,nz,lx,dx,dz
    USE constantes, ONLY : nx,nz,lx
    USE variables,  ONLY : it
    USE tableaux,   ONLY : vx,vx_a,vz,vz_a,pr,pr_a,tp
    IMPLICIT NONE
    INTEGER :: i,j,power
    REAL(KIND=8), DIMENSION(1:nx+1,1:nz+1) :: ww
    REAL(KIND=8)  :: xpos,zpos,dimd,wmin,wmax,deltax,deltaz
    CHARACTER(LEN=2)   :: name  
    CHARACTER(LEN=10)   :: iti
    CHARACTER(LEN=12)   :: name1
!===============================================================================
!
!  -K     More PostScript code will be appended later
!  -O     Selects Overlap plot mode
    WRITE(iti,'(i10.10)') it
!
    dimd=0.2*lx
    deltax=0.
    deltaz=5.
!
    name='vx'
    ww(1:nx+1,1:nz+1)=(vx(1:nx+1,1:nz+1)+vx(0:nx,1:nz+1)+vx(1:nx+1,0:nz)+vx(0:nx,0:nz))/4.
!    ww(1:nx+1,1:nz+1)=vz(1:nx+1,1:nz+1)
    DO i=1,nx+1
      DO j=1,nz+1
        IF (abs(ww(i,j)) < 1.e-12) ww(i,j)=0.
      ENDDO
    ENDDO
    xpos=2.
    zpos=22.
    name1=name//iti
    wmin=minval(ww(1:nx+1,1:nz+1))
    wmax=maxval(ww(1:nx+1,1:nz+1))
    CALL GMT_format(ww,nx+1,nz+1,power,100,name1,wmin,wmax)
!    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-P -K >','-O -K>>')
    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-P -K >',' -O >> ')

    name='vz'
    ww(1:nx+1,1:nz+1)=vz(1:nx+1,1:nz+1)
!    ww(1:nx+1,1:nz+1)=vz(1:nx+1,1:nz+1)-vz_a(1:nx+1,1:nz+1)
    DO i=1,nx+1
      DO j=1,nz+1
        IF (abs(ww(i,j)) < 1.e-12) ww(i,j)=0.
      ENDDO
    ENDDO
    xpos=xpos
    zpos=zpos-deltaz
    name1=name//iti
    wmin=minval(ww(1:nx+1,1:nz+1))
    wmax=maxval(ww(1:nx+1,1:nz+1))
    CALL GMT_format(ww,nx+1,nz+1,power,100,name1,wmin,wmax)
!    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-O -K>>','-O -K>>')
    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-P -K> ',' -O >> ')
!
    name='pr'
    ww(1:nx+1,1:nz+1)=(pr(1:nx+1,1:nz+1)+pr(1:nx+1,0:nz))/2.
!    name='Wy'
!    ww(1:nx+1,1:nz+1)=(vx(1:nx+1,2:nz+2)-vx(1:nx+1,0:nz  ))/dz/2.&
!                     -(vz(2:nx+2,1:nz+1)-vz(0:nx  ,1:nz+1))/dx/2.
    DO i=1,nx+1
      DO j=1,nz+1
        IF (abs(ww(i,j)) < 1.e-12) ww(i,j)=0.
      ENDDO
    ENDDO
    xpos=2.
    zpos=zpos-deltaz
    name1=name//iti
    wmin=minval(ww(1:nx+1,1:nz+1))
    wmax=maxval(ww(1:nx+1,1:nz+1))
    CALL GMT_format(ww,nx+1,nz+1,power,100,name1,wmin,wmax)
!    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-O -K>>','-O -K>>')
    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-P -K> ',' -O >> ')
!
!    name='pr'
    name='tp'
!    ww(1:nx+1,1:nz+1)=(pr(1:nx+1,1:nz+1)+pr(1:nx+1,0:nz))/2.
!    ww(1:nx+1,1:nz+1)=(vx(2:nx+2,1:nz+1)-vx(0:nx,1:nz+1))/2./dx &
!                     +(vz(1:nx+1,2:nz+2)-vz(1:nx+1,0:nz))/2./dz
    ww(1:nx+1,1:nz+1)=(tp(1:nx+1,1:nz+1)+tp(1:nx+1,0:nz))/2.
!    ww(1:nx+1,1:nz+1)=(pr(1:nx+1,1:nz+1)+pr(1:nx+1,0:nz))/2. &
!                     -(pr_a(1:nx+1,1:nz+1)+pr_a(1:nx+1,0:nz))/2.
!    ww(1:nx+1,1:nz+1)=(pr_a(1:nx+1,1:nz+1)+pr_a(1:nx+1,0:nz))/2.
    DO i=1,nx+1
      DO j=1,nz+1
        IF (abs(ww(i,j)) < 1.e-12) ww(i,j)=0.
      ENDDO
    ENDDO
    xpos=xpos
    zpos=zpos-deltaz
    name1=name//iti
    wmin=minval(ww(1:nx+1,1:nz+1))
    wmax=maxval(ww(1:nx+1,1:nz+1))
    CALL GMT_format(ww,nx+1,nz+1,power,100,name1,wmin,wmax)
!    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-O -K>>',' -O >> ')
    CALL GMT_dessin(name,power,iti,xpos,zpos,dimd,wmin,wmax,'-P -K> ',' -O >> ')
!
  END SUBROUTINE dessin
!
!
  SUBROUTINE GMT_format(ww,nx,nz,power,ncl,name1,wmin,wmax)
!===============================================================================
    USE constantes, ONLY : lx,lz
    IMPLICIT NONE
    INTEGER,      INTENT(IN)    :: nx,nz,ncl
    REAL(KIND=8), INTENT(INOUT) :: wmin,wmax
    REAL(KIND=8), INTENT(INOUT) :: ww(1:nx,1:nz)
    INTEGER :: i,j,k,itergmt,R,G,B,ncol,power
    CHARACTER(LEN=12),  INTENT(IN)    :: name1
    CHARACTER(LEN=16) ::  name2
    REAL(KIND=8) :: x,z,gap,zcl(1:3)
    REAL(KIND=8) :: val(0:250)
    REAL(KIND=8) :: wmoy(1:nz)
!===============================================================================
!
    name2=name1//'_gmt'
    power=0
    IF (abs(wmax) < 1.e-14) wmax=0.
    IF (abs(wmin) < 1.e-14) wmin=0.
    IF (abs(wmax-wmin) < 1.e-14) THEN
       IF (wmax < 0.)  THEN
         wmin=wmin*1.1
         wmax=wmax*0.9
       ELSEIF (wmin > 0.) THEN
         wmin=wmin*0.9 
         wmax=wmax*1.1
       ELSEIF (wmin == 0.) THEN
         wmin=-1.
         wmax=+1.
       ENDIF
    ENDIF
    power=NINT(log10(abs(wmax-wmin)))
    ww(1:nx,1:nz)=ww(1:nx,1:nz)/(10.**power)
    wmin=MINVAL(ww(1:nx,1:nz))
    wmax=MAXVAL(ww(1:nx,1:nz))
!
    IF (wmin*wmax < 0.) THEN
       IF ((wmin .lt. 0.) .AND. (abs(wmin) .gt. abs(wmax))) wmax=-wmin
       IF ((wmin .lt. 0.) .AND. (abs(wmax) .gt. abs(wmin))) wmin=-wmax
    ENDIF
    IF (wmin == wmax) wmax=wmin+1.
    IF ((abs(wmin)<1) .AND. (abs(wmax) <1)) THEN
      power=power-1
      ww(1:nx,1:nz)=ww(1:nx,1:nz)*10. 
      wmin=MINVAL(ww(1:nx,1:nz))
      wmax=MAXVAL(ww(1:nx,1:nz))
    ENDIF
!
    OPEN(1,FORM='formatted',FILE=name2)
    DO i=1,nx
       DO j=1,nz
          x=lx*(float(i-1))/float(nx-1)
          z=lz-lz*(float(j-1))/float(nz-1)
          WRITE(1,'(3e18.10)') x,z,ww(i,j)
       ENDDO
    ENDDO
    CLOSE(1)
!
    gap=(wmax-wmin)/float(ncl)
    ncol=255
    IF (gap /= 0.) THEN 
       itergmt=NINT((wmax-wmin)/gap)
       zcl(1)=wmin+(wmax-wmin)/4.
       zcl(2)=wmin+(wmax-wmin)/2.
       zcl(3)=wmin+3.*(wmax-wmin)/4.
       DO i=0,itergmt
          val(i)=wmin+float(i)*(wmax-wmin)/float(itergmt)
       ENDDO
    ENDIF
    OPEN(1,FORM='formatted',FILE=name2//'.cpt')
    WRITE(1,'(a)') '#   cpt file'
    WRITE(1,'(a)') '#COLOR_MODEL = RGB'
    WRITE(1,'(a,i3)') '#',power
    IF (gap /= 0.) THEN
       DO i=0,itergmt-1
          IF (val(i) .le. zcl(1)) THEN
             R=0
             G=ncol-NINT(float(ncol)*(val(i)-wmin)/(zcl(1)-wmin))
             B=ncol
          ELSEIF (val(i) .le. zcl(2)) THEN
             R=0
             G=0
             B=ncol-NINT(float(ncol)*(val(i)-zcl(1))/(zcl(2)-zcl(1)))
          ELSEIF (val(i) .le. zcl(3)) THEN
             R=NINT(float(ncol)*(val(i)-zcl(2))/(zcl(3)-zcl(2)))
             G=0
             B=0
          ELSE
             R=ncol
             G=NINT(float(ncol)*(val(i)-zcl(3))/(wmax-zcl(3)))
             B=0
          ENDIF
          WRITE(1,'(2(e16.8,3(i7)))') val(i),R,G,B,val(i+1),R,G,B
       ENDDO
    ELSE
       R=0
       G=ncol
       B=ncol
       WRITE(1,'(2(e16.8,3(i7)))') wmin,R,G,B,wmax+1.,R,G,B
    ENDIF
    WRITE(1,'(a,3(i4))') 'B',0,255,255
    WRITE(1,'(a,3(i4))') 'F',R,G,B
    WRITE(1,'(a,3(i4))') 'N',128,128,128
    CLOSE(1)
!
    DO k=1,nz
       wmoy(k)=SUM(ww(1:nx,k))/float(nx)
    ENDDO
    OPEN(1,FORM='formatted',FILE=name2//'.tmoy')
    DO k=1,nz
       WRITE(1,'(3e18.10)') wmoy(k),float(k-1)/float(nz-1)
    ENDDO
    CLOSE(1)
!     
  END SUBROUTINE GMT_format
!
!
  SUBROUTINE GMT_dessin(name2,power,iti,xpos,ypos,dimd,wmin,wmax,ph1,ph2)
!===============================================================================
    USE constantes, ONLY : dx,dz,lx,lz
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: xpos,ypos,dimd,wmin,wmax
    CHARACTER(LEN=10),  INTENT(IN) :: iti
    CHARACTER(LEN=7),  INTENT(IN) :: ph1,ph2
    CHARACTER(LEN=2),  INTENT(IN) :: name2
    REAL(KIND=8) :: xech,dime,gap
    INTEGER :: power
    CHARACTER(LEN=16) :: name1
!===============================================================================
    name1=name2//iti//'_gmt'
    xech=xpos+dimd*lx/2.
    dime=4.*dimd*lx/5.
    gap=(wmax-wmin)/2.
    gap=float(INT(10.*gap))/10.
!
    WRITE(95,'(4(a,e14.8),a)')'gmt xyz2grd -R0./0./',lx,'/',lz,&
    'r -I',dx,'/',dz,' '//name1//' -G'//name1//'.grd'
! 
    WRITE(95,'(6(a,e14.8),a)')&
    'gmt grdimage -B1.::/1.::WSne -Jx',dimd,&
    '/',dimd,' -R0./',lx,'/0./',&
    lz,' -E300 -Xa',xpos,' -Ya',ypos,' -C'//name1//'.cpt '//&
    name1//'.grd '//ph1//name2//iti//'.ps'
!
    IF (power == 0) THEN
       WRITE(95,'(a,e9.3,a,e9.3,a,i1.1,a,e9.3,a,e9.3,a)')&
       'gmt psscale -P -D0./0./',dime,'/0.5ch -B',gap,&
       ':@%29%'//name2//':/:10@+',power,'@+:WS -Xa',xech,&
       ' -Ya',ypos-1.,' -C'//name1//'.cpt '//ph2//name2//iti//'.ps'
    ELSEIF (power < 0) THEN
       WRITE(95,'(a,e9.3,a,e9.3,a,i3.2,a,e9.3,a,e9.3,a)')&
       'gmt psscale -P -D0./0./',dime,'/0.5ch -B',gap,&
       ':@%29%'//name2//':/:10@+',power,'@+:WS -Xa',xech,&
       ' -Ya',ypos-1.,' -C'//name1//'.cpt '//ph2//name2//iti//'.ps'
    ELSEIF (power > 0) THEN
       WRITE(95,'(a,e9.3,a,e9.3,a,i2.2,a,e9.3,a,e9.3,a)')&
       'gmt psscale -P -D0./0./',dime,'/0.5ch -B',gap,&
       ':@%29%'//name2//':/:10@+',power,'@+:WS -Xa',xech,&
       ' -Ya',ypos-1.,' -C'//name1//'.cpt '//ph2//name2//iti//'.ps'
    ENDIF
    WRITE(95,'(a)') 'ps2epsi '//name2//iti//'.ps'
!
  END SUBROUTINE GMT_dessin
!
END MODULE mod_figure
