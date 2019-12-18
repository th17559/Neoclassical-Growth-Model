
!=============================================================================
! Fortran 90 code for example in Lecture 2 (deterministic growth model)
! Topics in Advanced Macroeconomics, University of Konstanz
! Created by Haomin Wang
!=============================================================================


!----DECLARATION SECTION-----!
! - Procedure type (program) and name (main). 
! Only one "program" per set of code. 
program  main
! - "use" statement gives the current procedure access other modules
use globalvar
use solver
! - Always declare "implicit none"
implicit none
! - Then, for each variable/array, declare 
!		type, dimention, size, and other attributes. 
real(8) :: tic,toc

!----EXECUTION SECTION-----!

! "tic" start time of execution
call cpu_time(tic)

! initial is a subroutine in module "globalvar"
call initial
! solution is a subroutine in module "solver"
call solution

! Print execution time 
call cpu_time(toc)
PRINT*,'--------------------------------------------------'
PRINT*,'execution time =', toc-tic, 'sec.'
PRINT*,'--------------------------------------------------'


!----TERMINATION SECTION-----!
end program main








