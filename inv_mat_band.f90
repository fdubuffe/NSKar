MODULE inv_mat_band
CONTAINS
!
!
  SUBROUTINE tridiag(a,b,c,r,n)
!===============================================================================
    IMPLICIT NONE
    INTEGER                      :: i,n
    REAL(KIND=8)                 :: bet
    REAL(KIND=8), DIMENSION(1:n) :: a,b,c,u,r,gam
!===============================================================================
!
    bet=b(1)
    u(1)=r(1)/bet
    DO i=2,n
       gam(i)=c(i-1)/bet
       bet=b(i)-a(i)*gam(i)
       u(i)=(r(i)-a(i)*u(i-1))/bet
    END DO
    DO i=n-1,1,-1
       u(i)=u(i)-gam(i+1)*u(i+1)
    END DO
    r(1:n)=u(1:n)
!
  END SUBROUTINE tridiag
!
!
  SUBROUTINE tridiag_nr(a,b,c,r,u,n)
!=================================================================================================
! Numerical Recipes
! Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
! serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N,
! while a and c (oﬀ-diagonal elements) are size N − 1.
!=================================================================================================
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a,b,c,r
  REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: u
  REAL(KIND=8), DIMENSION(size(b)) :: gam
  INTEGER :: n,j
  REAL(KIND=8) :: bet
!=================================================================================================
    bet=b(1)
    u(1)=r(1)/bet
    do j=2,n
       gam(j)=c(j-1)/bet
!       bet=b(j)-a(j-1)*gam(j)
       bet=b(j)-a(j)*gam(j)
!       u(j)=(r(j)-a(j-1)*u(j-1))/bet
       u(j)=(r(j)-a(j)*u(j-1))/bet
    end do
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridiag_nr
!
!
  SUBROUTINE solve_tridiag(a,b,c,v,n)
!===============================================================================
    IMPLICIT NONE 
    INTEGER,INTENT(IN) :: n
    REAL(KIND=8),DIMENSION(n),INTENT(IN)    :: a,b,c
    REAL(KIND=8),DIMENSION(n),INTENT(INOUT) :: v
    REAL(KIND=8),DIMENSION(n) :: x,bp,vp
    REAL(KIND=8) :: m
    INTEGER :: i
!===============================================================================
    bp(1)=b(1)
    vp(1)=v(1)
! 
    DO i=2,n
       m=a(i)/bp(i-1)
       bp(i)=b(i)-m*c(i-1)
       vp(i)=v(i)-m*vp(i-1)
    ENDDO
    x(n)=vp(n)/bp(n)
    DO i=n-1,1,-1
       x(i)=(vp(i)-c(i)*x(i+1))/bp(i)
    ENDDO
    v=x
! 
  END SUBROUTINE solve_tridiag
!
!
  SUBROUTINE cyclic(a,b,c,alpha,beta,r,n)
!===============================================================================
    IMPLICIT NONE
    INTEGER :: n
    REAL(KIND=8) :: alpha,beta,fact,delta
    REAL(KIND=8), DIMENSION(1:n) :: a,b,bb,c,r,u
!===============================================================================
    delta=-b(1)
    bb(1)=b(1)-delta
    bb(n)=b(n)-alpha*beta/delta
    bb(2:n-1)=b(2:n-1)
    CALL solve_tridiag(a,bb,c,r,n)
    u(1)=delta
    u(n)=alpha
    u(2:n-1)=0.
    call solve_tridiag(a,bb,c,u,n)
    fact=(r(1)+beta*r(n)/delta)/(1.+u(1)+beta*u(n)/delta)
    r(1:n)=r(1:n)-fact*u(1:n)
  END SUBROUTINE cyclic
!
!
  SUBROUTINE cyclic_nr(a,b,c,alpha,beta,r,x,n)
!=================================================================================================
! Numerical Recipes
! Solves for a vector x(1:n) the “cyclic” set of linear equations given by equation (2.7.9).
! a, b, c, and r are input vectors, while alpha and beta are the corner entries in the matrix.
! The input is not modiﬁed.
!=================================================================================================
    REAL(KIND=8), DIMENSION(:), INTENT(IN):: a,b,c,r
    REAL(KIND=8), INTENT(IN) :: alpha,beta
    REAL(KIND=8), DIMENSION(:), INTENT(OUT):: x
    INTEGER :: n
    REAL(KIND=8) :: fact,gamma
    REAL(KIND=8), DIMENSION(size(x)) :: bb,u,z
!=================================================================================================
    gamma=-b(1)
    bb(1)=b(1)-gamma
    bb(n)=b(n)-alpha*beta/gamma
    bb(2:n-1)=b(2:n-1)
    call tridiag_nr(a,bb,c,r,x,n)
!    call tridag_nr(a(2:n),bb,c(1:n-1),r,x,n)
    u(1)=gamma
    u(n)=alpha
    u(2:n-1)=0.
    call tridiag_nr(a,bb,c,u,z,n)
!    call tridag(a(2:n),bb,c(1:n-1),u,z,n)
    fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
    x=x-fact*z
  END SUBROUTINE cyclic_nr
!
END MODULE inv_mat_band
