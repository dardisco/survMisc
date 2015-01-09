#include "R.h"
#include "Rinternals.h"
#include <R_ext/Utils.h>  
#ifdef USING_R
   /* typedef int Sint; */
#define S_EVALUATOR    /* Turn this into a "blank line" in R */
#else
/*
** Splus definitions, to use R type calls
*/
typedef long Sint;
/* 
**   I am using the 8.1 R*.h files courtesy of Bill Dunlap
*/
#ifdef  defineVar
#undef  defineVar
#endif
#define defineVar(a,b,c) ASSIGN_IN_FRAME(a,b, INTEGER_VALUE(c))
//
#ifdef  eval
#undef  eval
#endif
#define eval(a, b)  EVAL_IN_FRAME(a, INTEGER_VALUE(b))
//
/*
** These two refer to undefined functions, so use the 8.0.1 defs
*/
#ifdef asInteger
#undef asInteger
#endif
#define asInteger(a) INTEGER_VALUE(a)
#ifdef asReal
#undef asReal
#endif
#define asReal(a) NUMERIC_VALUE(a)
//
#endif
//
/*
** Memory defined with ALLOC is removed automatically by S.
**  That with "Calloc" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls.
*/
#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif
/*
** Prototype for callback function
**
*/
#ifdef USING_R
void cox_callback(int which, double *coef, double *first, double *second,
		  double *penalty, int *flag, int p, SEXP fexpr, SEXP rho);
#endif
//
/*
** Cox regression fit, replacement for coxfit2 in order
**    to be more frugal about memory: specificly that we 
**    don't make copies of the input data.
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       time(n)      :time of event or censoring for person i
**       status(n)    :status for the ith person    1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :offset for the linear predictor
**       weights(n)   :case weights
**       init         :initial estimate for the coefficients
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       chol_tol     : tolerance for the Cholesky decompostion
**       method       : 0=Breslow, 1=Efron
**       doscale      : 0=don't scale the X matrix, 1=scale the X matrix
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u(nv)        :score vector
**       imat(nv,nv)  :the variance matrix at beta=final
**                      (returned as a vector)
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
**       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       iter         :actual number of iterations used
**
**  work arrays
**       mark(n)
**       wtave(n)
**       a(nvar), a2(nvar)
**       cmat(nvar,nvar)       ragged array
**       cmat2(nvar,nvar)
**       newbeta(nvar)         always contains the "next iteration"
**       maxbeta(nvar)         limits on beta
**
**  calls functions:  cholesky2, chsolve2, chinv2
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>
//#include "survS.h"
//#include "survproto.h"
//
double **dmatrix(double *array, int ncol, int nrow){
  int i;
  double **pointer;
  pointer = (double **) ALLOC(nrow, sizeof(double *));
  for (i=0; i<nrow; i++) {
    pointer[i] = array;
    array += ncol;
  }
  return(pointer);
}
//
int cholesky2(double **matrix, int n, double toler){
/* $Id: cholesky2.c 11357 2009-09-04 15:22:46Z therneau $
**
** subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**     toler     the threshold value for detecting "singularity"
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.
**    The lower triangle need not be filled in at the start.
**
**  Return value:  the rank of the matrix (non-negative definite), or -rank
**     it not SPD or NND
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau
*/
  double temp;
  int  i,j,k;
  double eps, pivot;
  int rank;
  int nonneg;
  //printf (" Cholesky decomposition ");
  nonneg=1;
  eps =0;
  for (i=0; i<n; i++) {
    if (matrix[i][i] > eps)  eps = matrix[i][i];
    for (j=(i+1); j<n; j++)  matrix[j][i] = matrix[i][j];
  }
  eps *= toler;
  //
  rank =0;
  for (i=0; i<n; i++) {
    pivot = matrix[i][i];
    if (pivot < eps) {
      matrix[i][i] =0;
      if (pivot < -8*eps) nonneg= -1;
    }
    else  {
      rank++;
      for (j=(i+1); j<n; j++) {
	temp = matrix[j][i]/pivot;
	matrix[j][i] = temp;
	matrix[j][j] -= temp*temp*pivot;
	for (k=(j+1); k<n; k++) matrix[k][j] -= temp*matrix[k][i];
      }
    }
  }
  return(rank * nonneg);
}
//
//
void chsolve2(double **matrix, int n, double *y){
/*  $Id: chsolve2.c 11376 2009-12-14 22:53:57Z therneau $
**
** Solve the equation Ab = y, where the cholesky decomposition of A and y
**   are the inputs.
**
** Input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**        y[n] contains the right hand side
**
**  y is overwriten with b
**
**  Terry Therneau
*/
  register int i,j;
  register double temp;
  //printf (" Cholesky solve Ab = y ");
  /*
  ** solve Fb =y
  */
  for (i=0; i<n; i++) {
    temp = y[i] ;
    for (j=0; j<i; j++)
      temp -= y[j] * matrix[i][j] ;
    y[i] = temp ;
  }
     /*
     ** solve DF'z =b
     */
  for (i=(n-1); i>=0; i--) {
    if (matrix[i][i]==0)  y[i] =0;
    else {
      temp = y[i]/matrix[i][i];
      for (j= i+1; j<n; j++)
	temp -= y[j]*matrix[j][i];
      y[i] = temp;
    }
  }
}
//
void chinv2(double **matrix , int n){
/* $Id: chinv2.c 11357 2009-09-04 15:22:46Z therneau $
**
** matrix inversion, given the FDF' cholesky decomposition
**
** input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**
** returned: the upper triangle + diagonal contain (FDF')^{-1}
**            below the diagonal will be F inverse
**
**  Terry Therneau
*/
  register double temp;
  register int i,j,k;
  //printf (" Cholesky FDF' invert  ");
  /*
  ** invert the cholesky in the lower triangle
  **   take full advantage of the cholesky's diagonal of 1's
  */
  for (i=0; i<n; i++){
    if (matrix[i][i] >0) {
      matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
      for (j= (i+1); j<n; j++) {
	matrix[j][i] = -matrix[j][i];
	for (k=0; k<i; k++)     /*sweep operator */
	  matrix[j][k] += matrix[j][i]*matrix[i][k];
      }
    }
  }
  /*
  ** lower triangle now contains inverse of cholesky
  ** calculate F'DF (inverse of cholesky decomp process) to get inverse
  **   of original matrix
  */
  for (i=0; i<n; i++) {
    if (matrix[i][i]==0) {  /* singular row */
      for (j=0; j<i; j++) matrix[j][i]=0;
      for (j=i; j<n; j++) matrix[i][j]=0;
    }
    else {
      for (j=(i+1); j<n; j++) {
	temp = matrix[j][i]*matrix[j][j];
	if (j!=i) matrix[i][j] = temp;
	for (k=i; k<j; k++)
	  matrix[i][k] += temp*matrix[j][k];
      }
    }
  }
}
//
SEXP cf6(SEXP maxiter2,  SEXP time2,   SEXP status2, 
	 SEXP covar2,    SEXP offset2, SEXP weights2,
	 SEXP strata2,   SEXP method2, SEXP eps2, 
	 SEXP toler2,    SEXP ibeta,    SEXP doscale2) {
  int i, j, k, person;
  //printf("cf6 \n");
  double **covar, **cmat, **imat;  /*ragged arrays */
  double  wtave;
  double *a, *newbeta;
  double *a2, **cmat2;
  double *scale;
  double  denom=0, zbeta, risk;
  double  temp, temp2;
  int     ndead;  /* number of death obs at a time point */
  double  tdeath=0;  /* ndead= total at a given time point, tdeath= all */
  double  newlk=0;
  double  dtime, d2;
  double  deadwt;  /*sum of case weights for the deaths*/
  double  efronwt; /* sum of weighted risk scores for the deaths*/
  int     halving;    /*are we doing step halving at the moment? */
  int     nrisk;   /* number of subjects in the current risk set */
  double  *maxbeta;
  /* copies of scalar input arguments */
  int     nused, nvar, maxiter;
  int     method;
  double  eps, toler;
  int doscale;
  /* vector inputs */
  double *time, *weights, *offset;
  int *status, *strata;
  /* returned objects */
  SEXP imat2, means2, beta2, u2, loglik2;
  double *beta, *u, *loglik, *means;
  SEXP sctest2, flag2, iter2;
  double *sctest;
  int *flag, *iter;
  SEXP rlist, rlistnames;
  int nprotect;  /* number of protect calls I have issued */
  /* get local copies of some input args */
  nused = LENGTH(status2);
  nvar  = ncols(covar2);
  //printf ("nvar = %d \n", nvar);
  method = asInteger(method2);
  maxiter = asInteger(maxiter2);
  eps  = asReal(eps2);     /* convergence criteria */
  toler = asReal(toler2);  /* tolerance for cholesky */
  doscale = asInteger(doscale2);
  time = REAL(time2);
  weights = REAL(weights2);
  offset= REAL(offset2);
  status = INTEGER(status2);
  strata = INTEGER(strata2);
  /*
  **  Set up the ragged arrays and scratch space
  **  Normally covar2 does not need to be duplicated, even though
  **  we are going to modify it, due to the way this routine was
  **  was called.  In this case NAMED(covar2) will =0
  */
  nprotect = 0;
  if (NAMED(covar2) > 0){
    PROTECT(covar2 = duplicate(covar2)); 
    nprotect++;
  }
  covar = dmatrix(REAL(covar2), nused, nvar);
  //printf ("covar = %g ", **covar);
  for (i=0; i<nvar; i++){
    for (j=0; j<nused; j++){
      //printf ("covar[%i][%i] = %g ", i, j, covar[i][j]);
    }}; //printf("\n");
  PROTECT(imat2 = allocVector(REALSXP, nvar*nvar)); 
  nprotect++;
  imat = dmatrix(REAL(imat2),  nvar, nvar);
  a = (double *) R_alloc(2 * nvar * nvar + 5 * nvar, sizeof(double));
  for (i=0; i<nvar; i++) //printf ("a[%i] = %g", i, a[i]);
  //printf("\n");
  newbeta = a + nvar;
  for (i=0; i<nvar; i++) //printf ("newbeta[%i] = %g ", i, newbeta[i]);
  a2 = newbeta + nvar;
  for (i=0; i<nvar; i++) //printf ("a2[%i] = %g ", i, a2[i]);
  //printf("\n");
  maxbeta = a2 + nvar;
  for (i=0; i<nvar; i++){
    //printf ("i, maxbeta = %i, %g", i, maxbeta[i]);
  }; //printf("\n");
  scale = maxbeta + nvar;
  //printf(" scale = %g ", *scale);
  cmat = dmatrix(scale + nvar, nvar, nvar);
  //printf (" cmat = %g ", **cmat);
  int n = sizeof(**cmat) / sizeof(*cmat);
  //printf (" cmatsize = %d ", n);
  cmat2 = dmatrix(scale + nvar +nvar*nvar, nvar, nvar);
  //printf (" cmat2 = %g \n", **cmat2);
  /* 
  ** create output variables
  */ 
  PROTECT(beta2 = duplicate(ibeta));
  beta = REAL(beta2);
  PROTECT(means2 = allocVector(REALSXP, nvar));
  means = REAL(means2);
  PROTECT(u2 = allocVector(REALSXP, nvar));
  u = REAL(u2);
  PROTECT(loglik2 = allocVector(REALSXP, 2)); 
  loglik = REAL(loglik2);
  PROTECT(sctest2 = allocVector(REALSXP, 1));
  sctest = REAL(sctest2);
  PROTECT(flag2 = allocVector(INTSXP, 1));
  flag = INTEGER(flag2);
  PROTECT(iter2 = allocVector(INTSXP, 1));
  iter = INTEGER(iter2);
  nprotect += 7;
  /*
  ** Subtract the mean from each covar, as this makes the regression
  **  much more stable.
  */
  tdeath=0; temp2=0;
  for (i=0; i<nused; i++) {
    temp2 += weights[i];
    tdeath += weights[i] * status[i];
  }
  //printf ("tdeath = %g \n", tdeath);	
  for (i=0; i<nvar; i++) scale[i] = 1.0;
  //printf ("\nInitial iteration step\n");	
  /*
  ** do the initial iteration step
  */
  strata[nused-1] =1;
  for (i=0; i<nused; i++) {
    //printf ("strata[%i], %d", i, strata[i]);	
  };
  //printf ("\n");
  loglik[1] =0;
  for (i=0; i<nvar; i++) {
    u[i] = 0;
    a2[i] = 0;
    for (j=0; j<nvar; j++) {
      imat[i][j] = 0;
      cmat2[i][j] = 0;
    }
  }
  //printf ("u = %g, a2 = %g \n", *u, *a2);	
  for (person=nused-1; person>=0; ) {
    if (strata[person] == 1) {
      nrisk = 0 ;  
      denom = 0;
      for (i=0; i<nvar; i++) {
	a[i] = 0;
	for (j=0; j<nvar; j++) cmat[i][j] = 0;
      }
    }
    dtime = time[person];
    //printf ("\n\ndtime = %g \n", dtime);
    ndead = 0; /*number of deaths at this time point */
    deadwt = 0;  /* sum of weights for the deaths */
    efronwt= 0;  /* sum of weighted risks for the deaths */
    while(person >=0 && time[person]==dtime) {
      //printf ("person = %d \n", person);	
      /* walk through the this set of tied times */
      nrisk++;
      zbeta = offset[person];    /* form the term beta*z (vector mult) */
      for (i=0; i<nvar; i++) //printf (" beta[%i] = %g ", i, *beta);
      //printf("\n");	
      for (i=0; i<nvar; i++){
	//printf (" cov[%i, person] = %g ", i, covar[i][person]);	
	zbeta += beta[i]*covar[i][person];
      }; //printf("\n");
      //printf ("zbeta = %g ", zbeta);	
      risk = exp(zbeta) * weights[person];
      //printf ("risk = %g ", risk);	
      denom += risk;
      //printf ("denom = |%g|\n", denom);	
      /* a is the vector of weighted sums of x, cmat sums of squares */
      for (i=0; i<nvar; i++){ 
	a[i] += risk*covar[i][person];
	//printf ("a[%i] = %g ", i, a[i]);	
	for (j=0; j<=i; j++)
	  cmat[i][j] += risk*covar[i][person]*covar[j][person];
      }
      //printf("\n");
      for (i=0; i<nvar; i++) {
	for (j=0; j<=i; j++){
	  //printf ("cmat[%i][%i] = %g ", i, j, cmat[i][j]);
	}}; //printf ("\n");	
      //
      if (status[person]==1) {
	ndead++;
	deadwt += weights[person];
	efronwt += risk;
	//printf ("efronwt = %g ", efronwt);	
	loglik[1] += weights[person]*zbeta;
	//printf ("loglik[1] = %g ", loglik[1]);
	for (i=0; i<nvar; i++){ 
	  u[i] += weights[person]*covar[i][person];
	  //printf("u[%i] = %g ", i, u[i]);
	}			    
	if (method==1) { /* Efron */
	  for (i=0; i<nvar; i++) {
	    a2[i] +=  risk*covar[i][person];
       ////printf ("a2[i] (i in nvar) = |%g|\n", a2[i]);	
	    for (j=0; j<=i; j++){
	      cmat2[i][j] += risk*covar[i][person]*covar[j][person];
	    }
	    ////printf ("cmat2 = |%g|\n", cmat2[i][j]);	
	  }
	}
      };
      //	    
      // //printf ("While loop: person before = %d ", person);
      person--;
      // //printf ("While loop: person after = %d \n", person);
      if (strata[person]==1) break;  /*ties don't cross strata */
    }
    //printf ("ndead = %d \n", ndead);	
    if (ndead >0) {  /* we need to add to the main terms */
      if (method==0) { /* Breslow */
	loglik[1] -= deadwt*log(denom);
	//printf ("loglik[1] = %g ", loglik[1]);	
	for (i=0; i<nvar; i++) {
	  temp2= a[i]/ denom;  /* mean */
	  //printf ("temp2[%i] = %g ", i, temp2);	
	}; //printf("\n");
	for (i=0; i<nvar; i++) {
	  u[i] -=  deadwt* temp2;
	  //printf ("u[%i] = %g ", i, u[i]);	
	}; //printf("\n");
	for (i=0; i<nvar; i++) {
	  for (j=0; j<=i; j++){
	    imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
	    //printf (" imat[%i][%i] = %g ", i, j, imat[j][i]);	
	  }}
      }
      else { /* Efron */
	/*
	** If there are 3 deaths we have 3 terms: in the first the
	** three deaths are all in, in the second they are 2/3
	** in the sums, and in the last 1/3 in the sum.  Let k go
	** from 0 to (ndead -1), then we will sequentially use
	** denom - (k/ndead)*efronwt as the denominator
	** a - (k/ndead)*a2 as the "a" term
	** cmat - (k/ndead)*cmat2 as the "cmat" term
	** and reprise the equations just above.
	*/
	for (k=0; k<ndead; k++) {
	  temp = (double)k/ ndead;
	  wtave = deadwt/ndead;
	  //printf ("wtave (deadwt/ndead) = %g \n", wtave);	
	  d2 = denom - temp*efronwt;
	  //printf ("d2 = denom-k/ndead = %g \n", d2);	
	  loglik[1] -= wtave* log(d2);
	  for (i=0; i<nvar; i++) {
	    temp2 = (a[i] - temp*a2[i])/ d2;
	    u[i] -= wtave *temp2;
	    for (j=0; j<=i; j++)
	      imat[j][i] +=  (wtave/d2) *
		((cmat[i][j] - temp*cmat2[i][j]) -
		 temp2*(a[j]-temp*a2[j]));
	    // //printf ("imat[j][i] = |%g|\n", imat[j][i]);	
	  }
	}
	for (i=0; i<nvar; i++) {
	  a2[i]=0;
	  for (j=0; j<nvar; j++) cmat2[i][j]=0;
	}
      }
    }
  } /* end  of accumulation loop */
  //printf ("\n\n End accumulation loop \n");	
  //
  loglik[0] = loglik[1]; /* save the loglik for iter 0 */
  //printf ("loglik = %g \n", loglik[1]);	
    /*
    ** Use the initial variance matrix to set a maximum coefficient
    **  (The matrix contains the variance of X * weighted number of deaths)
    */
 for (i=0; i<nvar; i++){ 
   maxbeta[i] = 20* sqrt(imat[i][i]/tdeath);
   //printf ("imat[%i][%i] = %g ", i, i, imat[i][i]);	
   //printf ("maxbeta[%i] = %g \n", i, maxbeta[i]);	
 }
 /* am I done?
 **   update the betas and test for convergence
 */
 for (i=0; i<nvar; i++){ /*use 'a' as a temp to save u0, for the score test*/
   a[i] = u[i];
   //printf ("a[%i] = %g ", i, a[i]);	
 }; //printf("\n");
 //
 for (i=0; i<nvar; i++){
   for (j=0; j<=i; j++){
     //printf ("imat[%i][%i] = %g ", i, j, imat[i][j]);
   }}; 
 //
    *flag = cholesky2(imat, nvar, toler);
    //printf ("flag = %i \n", *flag);
    //
for (i=0; i<nvar; i++){
for (j=0; j<=i; j++){
  //printf ("imat[%i][%i] = %g ", i, j, imat[i][j]);
 }}

    chsolve2(imat,nvar,a);        
for (i=0; i<nvar; i++){
for (j=0; j<=i; j++){
  //printf ("imat[%i][%i] = %g", i, j, imat[i][j]);
 }}

/* a replaced by  a *inverse(i) */
    for (i=0; i<nvar; i++){
      //printf ("a[%i] = |%g| ", i, a[i]);
    }	
    temp=0;
    for (i=0; i<nvar; i++){
      //printf ("u[%i] = %g ", i, u[i]);
	temp +=  u[i]*a[i];
	//printf ("st/temp = %g ", temp);
    }
    *sctest = temp;  /* score test */
    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
      newbeta[i] = beta[i] + a[i];
      //printf ("newbeta[%i] = %g \n", i, newbeta[i]);
    }
    //
    //
    if (maxiter==0) {
	chinv2(imat, nvar);
	for (i=0; i<nvar; i++) {
	    beta[i] *= scale[i];  /*return to original scale */
	    u[i] /= scale[i];
	    imat[i][i] *= scale[i]*scale[i];
	    for (j=0; j<i; j++) {
		imat[j][i] *= scale[i]*scale[j];
		imat[i][j] = imat[j][i];
		}
	    }
	goto finish;
    }
    //
    /*
    ** here is the main loop
    */
    //
    //
    //
    halving =0 ;             /* =1 when in the midst of "step halving" */
    //printf ("\n\nMAIN LOOP");
    for (*iter=1; *iter<= maxiter; (*iter)++) {
      //printf ("\n\n ITER = |%d | \n", *iter);		
      newlk =0;
      for (i=0; i<nvar; i++) {
	u[i] =0;
	for (j=0; j<nvar; j++)
	  imat[i][j] =0;
      }
      /*
      ** The data is sorted from smallest time to largest
      ** Start at the largest time, accumulating the risk set 1 by 1
      */
      for (person=nused-1; person>=0; ) {
	if (strata[person] == 1) { /* rezero temps for each strata */
	  denom = 0;
	  nrisk =0;
	  for (i=0; i<nvar; i++) {
	    a[i] = 0;
	    for (j=0; j<nvar; j++) cmat[i][j] = 0;
	  }
	}
	dtime = time[person];
	//printf ("\n\ndtime = %g \n", dtime);		
	deadwt =0;
	ndead =0;
	efronwt =0;
	while(person>=0 && time[person]==dtime) {
	  nrisk++;
	  zbeta = offset[person];
	  for (i=0; i<nvar; i++) //printf("newbeta[%i] = %g ", i, *newbeta);	
	  for (i=0; i<nvar; i++) //printf ("cov[%i] = %g ", i, covar[i][person]);	
	  for (i=0; i<nvar; i++) {
	    zbeta += newbeta[i]*covar[i][person];
	    //printf ("zbeta = %g ", zbeta);	
	  };
	  risk = exp(zbeta) * weights[person];
	  //printf ("risk = %g ", risk);	
	  denom += risk;
	  //printf ("denom = %g \n", denom);	
	  for (i=0; i<nvar; i++) {
	    a[i] += risk*covar[i][person];
	    //printf ("a[%i] = %g ", i, a[i]);	
	    for (j=0; j<=i; j++)
	      cmat[i][j] += risk*covar[i][person]*covar[j][person];
	  }
	  for (i=0; i<nvar; i++) {
	    for (j=0; j<=i; j++){
	      //printf ("cmat[%i][%i] = %g ", i, j, cmat[i][j]);
	    }}; //printf ("\n");	
	  //
	  if (status[person]==1) {
	    ndead++;
	    deadwt += weights[person];
	    newlk += weights[person] * zbeta;
	    //printf ("newlk = %g ", newlk);     
	    for (i=0; i<nvar; i++){ 
	      u[i] += weights[person] *covar[i][person];
	      //printf("u[%i] = %g ", i, u[i]);
	    }
	    if (method==1) { /* Efron */
	      efronwt += risk;
	      for (i=0; i<nvar; i++) {
		a2[i] +=  risk*covar[i][person];
		for (j=0; j<=i; j++)
				cmat2[i][j] += risk*covar[i][person]*covar[j][person];
	      }   
	    }
	  }
	  
		person--;
		if (strata[person]==1) break; /*tied times don't cross strata*/
  }
	//printf ("ndead = %d \n", ndead);	
	if (ndead >0) {  /* add up terms*/
	  if (method==0) { /* Breslow */
		    newlk -= deadwt* log(denom);
		    //printf ("newlk[1] = %g ", newlk);	
		    for (i=0; i<nvar; i++) {
		      temp2 = a[i] / denom; /* mean */
		      //printf ("temp2[%i] = %g ", i, temp2);	
			u[i] -= deadwt * temp2;
			//printf ("u[%i] = %g ", i, u[i]);
		    }; //printf("\n");
		    for (i=0; i<nvar; i++) {	
		      for (j=0; j<=i; j++){
			imat[j][i] +=  (deadwt/denom)*
			  (cmat[i][j] - temp2*a[j]);
			//printf("imat[%i][%i] = %g ", j, i, imat[j][i]);
		      }}; //printf("\n");
	  }
	  else  { /* Efron */
	    for (k=0; k<ndead; k++) {
	      temp = (double)k / ndead;
	      wtave= deadwt/ ndead;
	      d2= denom - temp* efronwt;
	      newlk -= wtave* log(d2);
	      for (i=0; i<nvar; i++) {
		temp2 = (a[i] - temp*a2[i])/ d2;
		u[i] -= wtave*temp2;
		for (j=0; j<=i; j++)
		  imat[j][i] +=  (wtave/d2)*
		    ((cmat[i][j] - temp*cmat2[i][j]) -
		     temp2*(a[j]-temp*a2[j]));
	      }
    		        }
	    //
	    for (i=0; i<nvar; i++) { /*in anticipation */
	      a2[i] =0;
		      for (j=0; j<nvar; j++) cmat2[i][j] =0;
	    }
	  }
	}
      }   /* end  of accumulation loop  */
      /* am I done?
      **   update the betas and test for convergence
      */
      *flag = cholesky2(imat, nvar, toler);
      //printf ("loglik[1] = %g \n", loglik[1]);	
      //printf ("newlk = %g \n", newlk);	
      if (fabs(1-(loglik[1]/newlk))<= eps && halving==0) { /* all done */
	loglik[1] = newlk;
	//
	for (i=0; i<nvar; i++){
	  for (j=0; j<nvar; j++){
	    //printf ("imat[%i][%i] = |%g| ", i, j, imat[i][j]);	
	  }}
	//    
	    chinv2(imat, nvar);     /* invert the information matrix */
	    //
	    for (i=0; i<nvar; i++){
	      for (j=0; j<nvar; j++){
		//printf ("imat[%i][%i] = |%g| ", i, j, imat[i][j]);	
	      }}
	    //
	    for (i=0; i<nvar; i++) {
	      beta[i] = newbeta[i]*scale[i];
	      u[i] /= scale[i];
	      imat[i][i] *= scale[i]*scale[i];
	      for (j=0; j<i; j++) {
		imat[j][i] *= scale[i]*scale[j];
		    imat[i][j] = imat[j][i];
	      }
	    }
	    goto finish;
      }
      //
      if (*iter== maxiter) break;  /*skip the step halving calc*/
      //
      if (newlk < loglik[1])   {    /*it is not converging ! */
	//printf ("NOT converging");
	halving =1;
	for (i=0; i<nvar; i++)
	  newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
      }
      else {
	halving=0;
	loglik[1] = newlk;
	for (i=0; i<nvar; i++) //printf("u[%i] = %g ", i, u[i]);
	chsolve2(imat,nvar,u);
        for (i=0; i<nvar; i++) //printf("u[%i] = %g ", i, u[i]);
	j=0;
	for (i=0; i<nvar; i++) {
	  beta[i] = newbeta[i];
	  newbeta[i] = newbeta[i] +  u[i];
	  if (newbeta[i] > maxbeta[i]) newbeta[i] = maxbeta[i];
	  else if (newbeta[i] < -maxbeta[i]) newbeta[i] = -maxbeta[i];
	}
      }
    }   /* return for another iteration */
    /*
    ** We end up here only if we ran out of iterations 
    */
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
	beta[i] = newbeta[i]*scale[i];
	u[i] /= scale[i];
	imat[i][i] *= scale[i]*scale[i];
	for (j=0; j<i; j++) {
	    imat[j][i] *= scale[i]*scale[j];
	    imat[i][j] = imat[j][i];
	    }
	}
    *flag = 1000;
    //
finish:
    /*
    ** create the output list
    */
    PROTECT(rlist= allocVector(VECSXP, 8));
    SET_VECTOR_ELT(rlist, 0, beta2);
    SET_VECTOR_ELT(rlist, 1, means2);
    SET_VECTOR_ELT(rlist, 2, u2);
    SET_VECTOR_ELT(rlist, 3, imat2);
    SET_VECTOR_ELT(rlist, 4, loglik2);
    SET_VECTOR_ELT(rlist, 5, sctest2);
    SET_VECTOR_ELT(rlist, 6, iter2);
    SET_VECTOR_ELT(rlist, 7, flag2);
    /* add names to the objects */
    PROTECT(rlistnames = allocVector(STRSXP, 8));
    SET_STRING_ELT(rlistnames, 0, mkChar("coef"));
    SET_STRING_ELT(rlistnames, 1, mkChar("means"));
    SET_STRING_ELT(rlistnames, 2, mkChar("u"));
    SET_STRING_ELT(rlistnames, 3, mkChar("imat"));
    SET_STRING_ELT(rlistnames, 4, mkChar("loglik"));
    SET_STRING_ELT(rlistnames, 5, mkChar("sctest"));
    SET_STRING_ELT(rlistnames, 6, mkChar("iter"));
    SET_STRING_ELT(rlistnames, 7, mkChar("flag"));
    setAttrib(rlist, R_NamesSymbol, rlistnames);
    //
    unprotect(nprotect+2);
    //printf("\n");
    return(rlist);
}
