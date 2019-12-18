
! WITH CUBIC SPLINE PROJECTION 

module solver
use globalvar
implicit none

!Derived data type for cubic spline (x,y=f(x),y2=f''(x))
type csvar 
	real(8) :: x 
	real(8) :: y
	real(8) :: y2
end type csvar

! Current states for rhs function. 
!	Cannot be directly passed through optimizer. 
!	Hence, pass it through the global variable rhs_k.
real(8) :: k_current, z_current
! Current value spline parameters
type(csvar), dimension(length_grid_k) :: value_current_1, value_current_2, value_h_1, value_h_2


CONTAINS
!----------------------------------------------
	!cubic spline evaluation.
	function cseval(inp,eval) 
	use spline
    implicit none
    type(csvar) :: inp(:)
    real(8) :: eval
    real(8) :: cseval
    call splint(inp%x,inp%y,inp%y2,eval,cseval)
	end function cseval
	
!----------------------------------------------	
	!cubic spline update.
	subroutine spline_update(inp) 
	use spline
	implicit none
	type(csvar), intent(inout) :: inp(:)
	! nkcspline: "not-a-knot" cubic spline is a subroutine
	!	in the module "spline"
	call nkcspline(inp%x,inp%y,inp%y2)
	end subroutine spline_update

!----------------------------------------------	  
  	! Evaluates the right hand side of the mapping. 
  	! Takes "..._current" as module-wide global input
	function rhs(kp) 
	implicit none
	real(8), intent(in) :: kp
	real(8) :: rhs,c
	integer :: z_state
	
	! Consumption
	c=z_current*k_current**a+(1.0d0-d)*k_current-kp
	c=max(c,tiny(1.0d0)) ! c must be strictly positive
	
	if (z_current==Zgrid(1)) then
		z_state = 1
	else
		z_state = 2
	end if
	! Right-hand side of value function mapping
	rhs=log(c)+b*(P(z_state,1)*cseval(value_current_1,kp) + P(z_state,2)*cseval(value_current_2,kp))
	
	! reverse the sign for minimizer
	rhs=-rhs

	end function rhs

!----------------------------------------------  
	! Bracket minimum such that on ouput 
	!	rhs(b)<rhs(a) and rhs(b)<rhs(c) 
	! fail = .true. if bracket is not found.
	subroutine bracket(a,b,c,fail) 
	implicit none
	real(8), intent(inout) :: a,b,c
	integer, parameter :: maxit=100
	integer :: it
	logical, intent(out) :: fail
	fail=.FALSE.
	it=0
	if(rhs(a)<rhs(c)) then
		do while (rhs(a)<=rhs(b) .and. it<maxit)
		  c=b
		  b=(a+c)/2.0d0
		  it=it+1
		end do
	else
		do while (rhs(c)<=rhs(b) .and. it<maxit)
		  a=b
		  b=(a+c)/2.0d0
		  it=it+1
		end do
	end if
	if(it == maxit) fail=.TRUE.
	end subroutine bracket

!----------------------------------------------  
	! The main solver procedure
	subroutine solution
	! The "use" statement gives the current subroutine 
	!	access to module "globalvar".
	USE optim
	INTEGER :: iter, index_k, index_z, index_kp, i, index_h, index_aux, z_s
	real(8) :: valdist, poldist,distvec(2), k, kp, c, lb,ub,gs,gk
	! place holders used in value function interations
	real(8),dimension(length_grid_k,2) ::  value_new, g_k, value_new_new
	LOGICAL :: fail
	
	! Initialize value function iteration trackers
	iter = 1 ! # of iterateration
	valdist = 1000.0d0 ! value function distance
	poldist = 1000.0d0 ! policy function distance
	
	! Initial guess value spline parameters
	! "value_current" is a module-wide global object
	value_current_1%x=Kgrid
	value_current_1%y=value(:,1)
	value_current_1%y2=0.0d0
	
	value_current_2%x=Kgrid
	value_current_2%y=value(:,2)
	value_current_2%y2=0.0d0
	

	! Start value function iteration loop
	do while (valdist >= valtoler &
		.and. iter<=1000)
		
		! Update "value_current" (value spline parameters)
		call spline_update(value_current_1)
		call spline_update(value_current_2)
		
		! Looping over the state space
		do index_k = 1, length_grid_k			! capital grid
				!***Update module-wide global state variable
				k_current = Kgrid(index_k)
				
				do index_z = 1, 2
					z_current = Zgrid(index_z)
				
					!***Find a bracket for the "golden" minimizer
					! guess lb: the lowest capital grid point
					lb = klb 
					! guess ub: maximum savings 
					ub = z_current*k_current**a+(1.0d0-d)*k_current - 1.0d-14
					gs = (lb+ub)/2.0d0
					! call procedure bracket such that
					!	rhs(gs)< rhs(lb) and rhs(gs)<rhs(ub)
					call bracket(lb,gs,ub,fail)
				
					!*** Minimize "rhs"
					if (fail) then ! solution is at a bound	
						lb = klb
						ub = z_current*k_current**a+(1.0d0-d)*k_current - 1.0d-14
					if (rhs(lb)<rhs(ub)) then
						g_k(index_k,index_z)=lb
						value_new(index_k,index_z)=-rhs(lb)
					else
						g_k(index_k,index_z)=ub
						value_new(index_z,1)=-rhs(ub)
					end if
					
					else ! if bracket is found, minimize using "golden" 
					value_new(index_k,index_z)=-golden(lb,gs,ub,rhs,1.0d-8,kp)
					g_k(index_k,index_z)=kp
					end if				


				end do ! end of productivity grid loop (state space)

				
		end do ! end of capital grid loop (state space)
	
		!***Distance (check convergence)
		valdist=dist(value_new,value)
		poldist=dist(g_k,g_k0)

		write(*,*) 'Iteration =',iter,'valdist =',valdist,'poldist =',poldist
		
		value_h_1%x=Kgrid
		value_h_1%y=value_new(:,1)
		value_h_1%y2=0.0d0
			
		value_h_2%x=Kgrid
		value_h_2%y=value_new(:,2)
		value_h_2%y2=0.0d0
		
		do index_h=1,(iter**2)
		
			call spline_update(value_h_1)
			call spline_update(value_h_2)
			
			! Looping over the state space
			do index_k = 1, length_grid_k

				k_current = Kgrid(index_k)
				
				do index_z = 1, 2
					z_current = Zgrid(index_z)
					
					gk = g_k(index_k,index_z)
					
					c = z_current*k_current**a+(1.0d0-d)*k_current - gk
					
					if (z_current==Zgrid(1)) then
					z_s = 1
					else
					z_s = 2
					end if
					
					value_new_new(index_k,index_z)=log(c)+b*(P(z_s,1)*cseval(value_h_1,gk) + P(z_s,2)*cseval(value_h_2,gk))
					
				end do
				
			end do

			value_h_1%y = value_new_new(:,1)
			value_h_2%y = value_new_new(:,2)
			
					
		end do
		
		value_new = value_new_new
			
	
		!***Update value, policy, and "iter"
		value_current_1%y = value_new(:,1)
		value_current_2%y = value_new(:,2)
		value=value_new
		g_k0=g_k
		iter = iter+1
		
		!***Write results to external files.
		open (UNIT=1,FILE='valiter_cspline',STATUS='replace')
			do index_k = 1, length_grid_k
				WRITE(UNIT=1,FMT=*) Kgrid(index_k),value(index_k,1),value(index_k,2),g_k0(index_k,1), g_k0(index_k,2)
			end do
		close (UNIT=1)
		
	end do ! end of Value function iteration loop
	

	end subroutine solution
!----------------------------------------------  

end module solver
