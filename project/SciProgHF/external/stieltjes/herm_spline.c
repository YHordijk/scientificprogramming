#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int test_sign(double,double);
static double max(double,double);
static double min(double,double);
static void herm_spline_nj(double X1, double X2, double F1, double F2, double D1,
    double D2, int NE, const double *XE, double *FE, int *NEXT, int *IERR);

/************************************************
 * Parameters:
 * N -- (input) number of data points, >= 2
 *
 * X -- (input) array (double) of N increasing 
 *    independent variable values
 *
 * F -- (input) array (double) of N dependent
 *    variable values
 *
 * D -- (input) array (double) of derivative
 *    values at the data points
 *
 * M -- (input) number of evaluation points, >= 1
 *
 * PX -- (input) array (double) of points at which
 *     the function is to be evaluated
 *
 * PF -- (output) array (double) of values of the
 *     cubic Hermite function
 * ************************************************/
void herm_spline_(const int *N, const double *X, const double *F,
    const double *D, const int *M, const double *PX, double *PF,
    int *IERR) {

  int ierc,ir,j,jfirst,nj;
  int next[2];

  *IERR = 0;
  jfirst = 0;
  for (ir=1; ir<*N; ir++) {
    /* locate all the PX points located in the interval (X[ir-1],X[ir]) */
    for (j=jfirst; j<*M; j++)
      if (PX[j] >= X[ir])
        break;
    if (ir == *N-1)
      j = *M;

    nj = j - jfirst;
    if (nj != 0) {
      herm_spline_nj(X[ir-1],X[ir],F[ir-1],F[ir],D[ir-1],D[ir],
          nj,PX+jfirst,PF+jfirst,next,&ierc);

      if (next[1] > 0) {
        /* There are next[1] points to the right of X[ir] in the current set of
         * PX points */
        if (ir == *N-1)
          *IERR = *IERR - next[1];
        else {
          *IERR = 6;
          return;
        }
      }

      if (next[0] > 0) {
        /* There are next[0] points to the left of X[ir-1] in the current set of
         * PX points */
        if (ir == 1)
          *IERR = *IERR - next[0];
        else {
          *IERR = 6;
          return;
        }
      }
    }
    jfirst = j;
    if (jfirst >= *M)
      break;
  }
  return;
}


/************************************************ 
 * Parameters:
 * N -- (input) number of data points; >= 2
 *
 * X -- (input) array (double) of N increasing 
 *    independent variable values
 *
 * F -- (input) array (double) of N dependent
 *    variable values
 *
 * D -- (output) array (double) of derivative
 *    values at the data points
 *
 * INFCD -- (input) increment between succesive 
 *    values in F and D (provided primarily for
 *    2D applications); INFCD.ge.1
 * ************************************************/
void herm_spline_init_(const int *N, const double *X, const double *F,
    double *D) {

  int ierr,i,nless1,tmp;
  double del1,del2,dmax,dmin,drat1,drat2,dsave,h1,h2,hsum,hsumt3,
         w1,w2;

  ierr = 0;
  nless1 = *N-1;
  h1 = X[1]-X[0];
  del1 = (F[1]-F[0])/h1;
  dsave = del1;
  /* N=2 -- use linear interpolation */
  if (nless1 <= 1) {
    D[0] = del1;
    D[nless1] = del1;
    return;
  }
  /* N>2 -- normal case */
  h2 = X[2]-X[1];
  del2 = (F[2]-F[1])/h2;
  /* Set D[0] via non-centered three point formula, adjusted to be
   * shape-preserving */
  hsum = h1+h2;
  w1 = (h1+hsum)/hsum;
  w2 = -h1/hsum;
  D[0] = w1*del1+w2*del2;
  if (test_sign(D[0],del1) < 0)
    D[0] = 0.0;
  else if (test_sign(del1,del2) < 0) {
    /* necessary if monotonicity switches */
    dmax = 3.0*del1;
    if (fabs(D[0]) > fabs(dmax))
      D[0] = dmax;
  }
  /* Loop through interior points */
  for (i=1; i<nless1; i++) {
    if (i != 1) {
      h1 = h2;
      h2 = X[i+1]-X[i];
      hsum = h1+h2;
      del1 = del2;
      del2 = (F[i+1]-F[i])/h2;
    }
    /* Set D[i]=0 unless data are strictly monotonic */
    D[i] = 0.0;
    tmp = test_sign(del1,del2);
    if (tmp < 0) {
      ierr--;
      dsave = del2;
    }
    else if (tmp == 0) {
      /* Count number of changes in direction of monotonicity */
      if (del2 != 0.0) {
        if (test_sign(dsave,del2)<0)
          ierr--;
        dsave = del2;
      }
    }
    else {
      /* Use Brodlie modification of Butland formula */
      hsumt3 = hsum+hsum+hsum;
      w1 = (hsum+h1)/hsumt3;
      w2 = (hsum+h2)/hsumt3;
      dmax = max(fabs(del1),fabs(del2));
      dmin = min(fabs(del1),fabs(del2));
      drat1 = del1/dmax;
      drat2 = del2/dmax;
      D[i] = dmin/(w1*drat1+w2*drat2);
    }
  }

  /* Set D(N-1) via non-centered three point formula, adjusted to be
   * shape-preserving */
  w1 = -h2/hsum;
  w2 = (h2+hsum)/hsum;
  D[*N-1] = w1*del1+w2*del2;
  if (test_sign(D[*N-1],del2) <= 0)
    D[*N-1] = 0.0;
  else if (test_sign(del1,del2) < 0) {
    dmax = 3.0*del2;
    if (fabs(D[*N-1]) > fabs(dmax))
      D[*N-1] = dmax;
  }
}


int test_sign(double x1, double x2) {
  double sign;
  if (x1 == 0.0 || x2 == 0.0)
    return(0);
  else {
    sign = x1/fabs(x1)*x2/fabs(x2);
    return((int) sign);
  }
}

double max(double x1, double x2) {
  if (x1 >= x2)
    return(x1);
  else
    return(x2);
}

double min(double x1, double x2) {
  if (x1 <= x2)
    return(x1);
  else
    return(x2);
}

/************************************************ 
 * Parameters:
 *
 * X1,X2 --  (input) endpoints of the interval of definition of cubic
 *
 * F1,F2 -- (input) values of the function at X1 and X2, respectively
 *
 * D1,D2 -- (input) values of derivative at X1 and X2, respectively
 *
 * NE -- (input) number of evaluation points
 *
 * XE -- (input) array (double) of points at which the function are to be
 * evaluated
 *
 * FE -- (output) array (double) of values of the cubic function at the XE
 * points
 *
 * NEXT -- (output) array indicating number of evaluation points:
 *    NEXT[0] = number of evaluation points to the left of the interval
 *    NEXT[1] = number of evaluation points to the right of the interval
 *
 * IERR -- (output) error flag
 * ************************************************/ 
void herm_spline_nj(double X1, double X2, double F1, double F2, double D1,
    double D2, int NE, const double *XE, double *FE, int *NEXT, int *IERR) {

  double c2,c3,del1,del2,delta,h,x,xma,xmi;
  int i;

  if (NE < 1) {
    *IERR = 1;
    return;
  }

  h = X2 - X1;
  if (h == 0.0) {
    *IERR = 2;
    return;
  }

  *IERR = 0;
  NEXT[0] = 0;
  NEXT[1] = 0;
  xmi = min(0.0,h);
  xma = max(0.0,h);
  delta = (F2-F1)/h;
  del1 = (D1-delta)/h;
  del2 = (D2-delta)/h;
  c2 = -(del1+del1+del2);
  c3 = (del1+del2)/h;
  for (i=0; i<NE; i++) {
    x = XE[i]-X1;
    FE[i] = F1+x*(D1+x*(c2+x*c3));
    if (x < xmi)
      NEXT[0]++;
    else if (x > xma)
      NEXT[1]++;
  }
  return;
}
