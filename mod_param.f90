MODULE constantes
  REAL(KIND=8), PARAMETER :: lc=2.**(3./2.),lx=2.*lc,lz=1.,pi=acos(-1.),&
                             Tt=-0.5,Tb=+0.5,T0=Tt,&
                             Ray=10.**(3./8.)*1.e+03,Pra=1.e+00,&
                             ampV=1.e+00,ampT=1.e-02,eps=1.e-20
  INTEGER,      PARAMETER :: nz=64,nx=NINT(lx/lz)*nz,new=1,md=MAX(nx+1,nz+1)
  REAL(KIND=8), PARAMETER :: dx=lx/float(nx),dz=lz/float(nz), &
                             kx0=2.*pi,kz0=1.*pi,alf0=(kx0**2+kz0**2)*Pra
!                             kx0=2.*pi,kz0=2.*pi,alf0=(kx0**2+kz0**2)*Pra
END MODULE constantes
!
MODULE variables
  INTEGER      :: it,itmax=1000000,itold,it_pr
  REAL(KIND=8) :: dt,alpha,temps,cour=1.e-05
END MODULE variables
!
MODULE tableaux
  USE constantes, ONLY : nx,nz
  REAL(KIND=8), DIMENSION(-1:nx+3,-1:nz+3) :: vx,vx0,vx_a,vz,vz0,vz_a,pr,pr0,pr_a,tp,tp0
END MODULE tableaux
!
MODULE fourier
  USE constantes, ONLY : nx,nz
  REAL(KIND=8), DIMENSION(2*(nx/2+1),nz) :: XF
  REAL(KIND=8), DIMENSION(100+2*nx)      :: TABLE
  REAL(KIND=8), DIMENSION((2*nx+4)*nz)   :: WORK
END MODULE fourier
