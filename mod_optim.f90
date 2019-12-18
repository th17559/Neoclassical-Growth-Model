MODULE optim 
  IMPLICIT NONE
  private
  public :: golden
CONTAINS
  !Numerical recipes simple bracketing minimizer.
  FUNCTION golden(ax,bx,cx,func,tol,xmin) 
    IMPLICIT NONE
    real(8), INTENT(IN) :: ax,bx,cx,tol
    real(8), INTENT(OUT) :: xmin
    real(8) :: golden
    INTERFACE
       FUNCTION func(x)
         IMPLICIT NONE
         real(8), INTENT(IN) :: x
         real(8) :: func
       END FUNCTION func
    END INTERFACE
    real(8), PARAMETER :: R=0.61803399d0,C=1.0d0-R
    real(8) :: f1,f2,x0,x1,x2,x3
    x0=ax
    x3=cx
    if (abs(cx-bx) > abs(bx-ax)) then
       x1=bx
       x2=bx+C*(cx-bx)
    else
       x2=bx
       x1=bx-C*(bx-ax)
    end if
    f1=func(x1)
    f2=func(x2)
    do
       if (abs(x3-x0) <= tol*(abs(x1)+abs(x2))) exit
       if (f2 < f1) then
          call shft3(x0,x1,x2,R*x2+C*x3)
          call shft2(f1,f2,func(x2))
       else
          call shft3(x3,x2,x1,R*x1+C*x0)
          call shft2(f2,f1,func(x1))
       end if
    end do
    if (f1 < f2) then
       golden=f1
       xmin=x1
    else
       golden=f2
       xmin=x2
    end if
  END FUNCTION golden

  SUBROUTINE shft2(a,b,c)
    real(8), INTENT(OUT) :: a
    real(8), INTENT(INOUT) :: b
    real(8), INTENT(IN) :: c
    a=b
    b=c
  END SUBROUTINE shft2

  SUBROUTINE shft3(a,b,c,d)
    real(8), INTENT(OUT) :: a
    real(8), INTENT(INOUT) :: b,c
    real(8), INTENT(IN) :: d
    a=b
    b=c
    c=d
  END SUBROUTINE shft3
END MODULE optim
