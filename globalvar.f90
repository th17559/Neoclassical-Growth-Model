
!----DECLARATION SECTION-----!
! - Procedure type (module) and name (globalvar). 
module globalvar
! - Always declare "implicit none"
implicit none
! - Then, for each variable/array, declare 
!		type, dimention, size, and other attributes. 

! model parameters (b=beta, d=delta, a=alpha)
real(8)  	:: b=0.99d0, d=0.025d0 , a=0.36d0
! bounds and increments of Kgrid
real(8),parameter  	:: klb=1.0d-5 , kub=55.0d0
! length of Kgrid
INTEGER,parameter 	:: length_grid_k = 500
! tolerance level
real(8)  	:: valtoler=1.0d-7, poltoler=1.0d-7
! declare Kgrid
real(8),DIMENSION(length_grid_k) :: Kgrid
! declare value and g_k0
real(8),DIMENSION(length_grid_k, 2) :: value, g_k0
! declare Zgrid
real(8),DIMENSION(2) :: Zgrid = (/1.05, 0.84/)
! declare transition matrix
real(8),DIMENSION(2,2) :: P = &
							RESHAPE((/0.977d0, 0.074d0, 0.023d0, 0.926d0/),(/2,2/))

!----EXECUTION SECTION-----!
! In a module, the execution section can only contain internal procedures.

CONTAINS
	
	!----DECLARATION SECTION OF INTERNAL PROCEDURE "initial"-----!
	subroutine initial
	implicit none
	INTEGER :: i,j
	
	!----EXECUTION SECTION OF INTERNAL PROCEDURE "initial"-----!	
	! Initialize capital grid
	do i=1,length_grid_k
		Kgrid(i) = log(klb) + (i-2)*(log(kub)-log(klb))/dble(length_grid_k-1)
		Kgrid(i) = exp(Kgrid(i))
	end do
	
	! Initial Guess for Value and Policy
	do  i=1,length_grid_k
		do	j=1,2
			value(i,j) = log(max(klb,Zgrid(j)*Kgrid(i)**a-d*Kgrid(i)))/(1.0d0-b)
			g_k0(i,j)=Kgrid(i)
		end do
	end do

	!----TERMINATION SECTION OF INITIAL PROCEDURE "initial"-----!
	end subroutine initial
	
	!----SECOND INTERNAL PROCEDURE-----!
	!metric for evaluation of distance between vectors.
	double precision function dist(inp1,inp2) 
	implicit none
	real(8), dimension(length_grid_k,2), intent(in) :: inp1, inp2
	real(8) :: avg
	!integer :: dim
	! Find the length of the input vector
	!dim=size(inp1)
	! Dist: sup-norm normalize by "avg"
	avg=(maxval(abs(inp1))+maxval(abs(inp2)))/dble(2)
	dist=maxval(abs(inp1-inp2))/avg
	!dist=maxval(abs(inp1-inp2))
	! If the sup norm is not well defined, set the distance
	!  to a large number (relative to your tolerance)
	if(isnan(dist)) dist=1.0d3
	! End of the function
	end function dist
	

!----TERMINATION SECTION-----!
end module globalvar