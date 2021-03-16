#ifndef LINSOLVE_H
#define LINSOLVE_H

/*
  Simple wrappers for LAPACK linear equation
  solvers. 
 */
  
/*
  Solve the linear system AX=rhs, using LAPACK. A is an NxN matrix
  stored in column-packed format, with a distance LDA between
  elements in each row. rhs is NxNRHS matrix with a distance LDrhs.
  
  If the solution was found linsolve() returns 0 and writes the
  solution into rhs. The A matrix is overwritten in the process. If
  you need the factorization of A you should call DGESV directly
  instead of using this wrapper.
*/
int linsolve(int N, int NRHS, double *A, int LDA,
	     double *rhs,  int LDrhs);

/*
  Like linsolve(), but solves the least square problem suitable for
  rank-deficit A. Uses DGELSD internally.  A is an MxN matrix in
  this case. Note that rhs is NxNRHS on input and MxNRHS on output.
  
  Singular values of A smaller than rcond will be treated as
  zero. Using rcond < 0 will set rcond to machine precision.
  
  This function may allocate quite a bit of memory for the SVD.
*/
int linsolve_least_squares(int M, int N, int NRHS, double *A, int LDA,
			   double *rhs,  int LDrhs, double rcond);


#endif
