##########################
# makefile for valiter in Lecture 2
# Topics in Advanced Macroeconomics, University of Konstanz
# Created by Haomin Wang
#########################


# Select Compiler
COMPILER =  gfortran

# Select files with source code (w.o projection)
SRCS= globalvar.f90\
	 solver.f90\
	main.f90

# Select files with source code (w. projection)
   SRCS_cspline= mod_optim.f90\
   spline.f90\
	globalvar.f90\
	solver_cspline.f90\
   main.f90


# command line
# W/o projection, use SRCS
# W. projection, use SRCS_cspline
gfortran:
#	$(COMPILER)  $(SRCS)   -g -fbacktrace  -fbounds-check  -o run
#	$(COMPILER)  $(SRCS_cspline)   -g -fbacktrace  -fbounds-check  -o run_cspline
	$(COMPILER)  $(SRCS_cspline)   -O1 -o run

