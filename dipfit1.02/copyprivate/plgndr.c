/*
 *  (c) 2002 Robert Oostenveld
 *
 * This is an attempt to make it execute faster, but I do not 
 * notice any difference. Probably the Matlab overhead is larger 
 * than the actual time spent in this function.
 */

#include <math.h>
#include "mex.h"
#include "matrix.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int l, m;
  double x, y;
  double *pd;

  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;

  if (nrhs != 3)
    mexErrMsgTxt ("invalid number of arguments for PLGNDR");

  l = mxGetScalar (prhs[0]);
  m = mxGetScalar (prhs[1]);
  x = mxGetScalar (prhs[2]);

  if (m < 0 || m > l || x > 1.0 || x < -1.0)
    mexErrMsgTxt ("Bad arguments in routine plgndr");

  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  pd = mxGetData (plhs[0]);

  /* this is from Numerical Recipes 2.0 */

  pmm=1.0;
  if (m > 0) {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
   for (i=1;i<=m;i++) {
     pmm *= -fact*somx2;
     fact += 2.0;
    }
  }

  if (l == m)
    pd[0] = pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      pd[0] = pmmp1;
    else {
      for (ll=m+2;ll<=l;ll++) {
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;
      }
      pd[0] = pll;
    }
  }

  return;
}

