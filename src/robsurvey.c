/* =========================================================================
 * PROJECT: RobSurv
 * SUBJECT: C functions
 * AUTHORS: Tobias Schoch, March 23, 2019
 * LICENSE: MIT
 * COMMENT:
 * ========================================================================= */
#include "robsurvey.h"

/*some macros*/
#define _WGT_HUBER(_x, _k) ((fabs(_x) >= _k) ? _k / fabs(_x) : 1.)
#define _WGT_HUBERasym(_x, _k) ((fabs(_x) >= _k) ? _k / _x : 1.)
#define _POWER2(_x) ((_x) * (_x))

/*declaration of 'local' functions (inline imperative is GCC specific)*/
static inline double euclid_metric(const double*, const double*, int)
    __attribute__((always_inline));

/*declaration of 'exported' functions*/
void fitwls(double*, double*, double*, double*, double*, int*, int*, int*);
void rwlslm(double*, double*, double*, double*, double*, double*, int*, int*,
      double*, double*, double*, int*, double*, int*);
void wmad(double*, double*, double*, int*);
void wquantile(double*, double*, double*, double*, int*);
void wmeantrimmed(double*, double*, double*, double*, double*, int*);
void wmeanwinsorized(double*, double*, double*, double*, double*, int*);


void part(double*, double, int, int*, int*);
void wpart(double*, double*, double, int, int*, int*);
void med(double*, int*, double*);
double select(double*, int, int);
double heap_select(double*, int, int);




/* rwls: iteratively reweighted least squares
 *---------------------------------------------------------------------------
 * PARAMETERS
 *   X         vectorized design matrix, array[n * p]
 *   y         response vectirls.exeor, array[n]
 *   w	      weight, array[n]
 *   resid     residuals vector, array[n]
 *   infl      influence function values, array[n]
 *   robwgt    robustness weights, array[n]
 *   ptr_n     (array dimensions)
 *   ptr_p     (array dimensions)
 *   ptr_k     r obustness tuning constant
 *   beta0     coefficient vector, array[p]
 *   ptr_irls.exescale scale estimate
 *   ptr_maxit max iterations (on input); iterations (on return)
 *   ptr_tol   numerical tolerance criterion (stoping rule in rwls updating
 *             rule)
 *   ptr_psi   psi-function (0 = Huber, 1 = asymmetric Huber)
 */
void rwlslm(double *X, double *y, double *w, double *resid, double *infl,
	    double *robwgt, int *ptr_n, int *ptr_p, double *ptr_k,
	    double *beta0, double *ptr_scale, int *ptr_maxit, double *ptr_tol,
	    int *ptr_psi)
{
   int iterations = 0, converged = 0;
   double *beta_new, *modweight;
   /*STEP 1: initialize beta by weighted least squares*/
   beta_new = (double* ) Calloc(*ptr_p, double);
   modweight = (double *) Calloc(*ptr_n, double);
   /*determine optimal size of array 'work' in 'fitwls'*/
   int lwork = -1;
   fitwls (X, y, w, resid, beta0, ptr_n, ptr_p, &lwork);
   /*compute wls fit (using the optimal size of 'work')*/
   fitwls (X, y, w, resid, beta0, ptr_n, ptr_p, &lwork);
   /*STEP 2: initialize scale estimate*/
   wmad(resid, w, ptr_scale, ptr_n);
   /*STEP 3: irls updating*/
   while(!converged && ++iterations < *ptr_maxit){
      /*STEP 3.1: update beta*/
      if(*ptr_psi == 0){ // Huber psi
	 for(int i = 0; i < *ptr_n; i++){
	    modweight[i] = w[i] * _WGT_HUBER( resid[i] / *ptr_scale, *ptr_k);
	 }
      }else{ // asymmetric Huber psi
	 for(int i = 0; i < *ptr_n; i++){
	    modweight[i] = w[i] * _WGT_HUBERasym( resid[i] / *ptr_scale, *ptr_k);
	 }
      }
      fitwls (X, y, modweight, resid, beta_new, ptr_n, ptr_p, &lwork);
      /*STEP 3.2: update A*/
      wmad(resid, w, ptr_scale, ptr_n);
      /*check for convergence*/
      converged = (euclid_metric(beta0, beta_new, *ptr_p) < *ptr_tol) ? 1: 0;
      /*prepare the next while 'run'*/
      Memcpy(beta0, beta_new, *ptr_p);
   }
   *ptr_maxit = (converged) ? iterations : 0;
   /*influence function values and robweights*/
   if(*ptr_psi == 0){ // Huber psi
      for(int j = 0; j < *ptr_n; j++){
	 robwgt[j] = _WGT_HUBER(resid[j] / *ptr_scale, *ptr_k);
	 infl[j] = robwgt[j] * resid[j];
      }
   }else{ // asymmetric Huber psi
      for(int j = 0; j < *ptr_n; j++){
	 robwgt[j] = _WGT_HUBERasym(resid[j] / *ptr_scale, *ptr_k);
	 infl[j] = robwgt[j] * resid[j];
      }
   }
   /*house keeping*/
   Free(beta_new); Free(modweight);
}
/* fitwls: weighted least squares estimate
 *--------------------------------------------------------------------------
 * PARAMETERS
 *   X           vectorized design matrix, array[n * p]
 *   y           response vector, array[n]
 *   w           weights vector, array[n]
 *   resid       residuals vector, array[n]
 *   beta0       coefficient vector, array[p]
 *   n, p        (array dimensions)
 *   ptr_lwork   size of array 'work' (if < 0, then 'dgels' determines
 *	        optimal size)
 *
 * DEPENDENCIES
 *     dgels    <R_ext/LAPLACK.h>
 *     dgemv    <R_ext/BLAS.h>
 */
void fitwls(double *X, double *y, double *w, double *resid, double *beta0,
	    int *ptr_n, int *ptr_p, int *ptr_lwork)
{
    /*define constants for the call of 'dgels'*/
    const int int_1 = 1;
    int info = 1;
    double *work, *wy, *wX, tmp;
    /*STEP 0: compute optimal size of array work*/
    if (*ptr_lwork < 0) {
        work = (double *) Calloc(*ptr_n, double);
    } else {
        work = (double *) Calloc(*ptr_lwork, double);
    }
    wy = (double *) Calloc(*ptr_n, double);
    wX   = (double *) Calloc(*ptr_n * *ptr_p, double);
    /*pre-multiply the design matrix and the response vector by sqrt(w)*/
    for (int i = 0; i < *ptr_n; i++) {
        tmp = sqrt(w[i]);
        wy[i] = y[i] * tmp;
        for (int j = 0; j < *ptr_p; j++) {
            wX[*ptr_n * j + i] = X[*ptr_n * j + i] * tmp;
        }
    }
    /*compute the (weighted) least squares estimate (LAPACK::dgels)*/
    F77_CALL(dgels)("N", ptr_n, ptr_p, &int_1, wX, ptr_n, wy, ptr_n, work,
        ptr_lwork, &info);
    /*STEP 1: compute weighted ls*/
    if( *ptr_lwork < 0 ) {
        *ptr_lwork = (int) work[0]; // optimal value of 'lwork'
    } else {
        const double double_minus1 = -1., double_1 = 1.;
        /*STEP 2: obtain 'betacoefficients'*/
        Memcpy(beta0, wy, *ptr_p);
        /*STEP 3: compute the residuals (BLAS::dgemv)*/
        Memcpy(resid, y, *ptr_n);
        F77_CALL(dgemv)("N", ptr_n, ptr_p, &double_minus1, X, ptr_n, beta0,
            &int_1, &double_1, resid, &int_1);
    }
    /*housekeeping*/
    Free(work); Free(wy); Free(wX);
}
/* euclid_metric: Euclidean metric
 * -------------------------------------------------------------------------
 * PARAMETERS
 *    x	    array
 *    y	    array
 *    n	    array dimension
 */
static inline double euclid_metric(const double *x, const double *y,
				   const int n)
{
    double s = 0.;
    for (int i = 0; i < n; i++) {
	    s += _POWER2(x[i] - y[i]);
    }
    return(sqrt(s));
}
/* wquantile: weighted lower quantile
 * -------------------------------------------------------------------------
 * PARAMETERS
 *    x		  data
 *    w		  weight
 *    ptr_prob	  probability of quantile
 *    ptr_quant	  weighted lower quantile (on return)
 *    ptr_n	  array dimension
 *
 * DEPENDENCIES
 *    rsort_with_index	<R_ext/Utils.h>
 */
void wquantile(double *x, double *w, double *ptr_prob, double *ptr_quant,
	       int *ptr_n)
{
   int i = 1;
   int *iarray;
   double s = 0.0, total = 0.0;
   double *x_cpy, *w_cpy;
   /*allocate array 'iarray'*/
   iarray = (int*) Calloc(*ptr_n, int);
   /*populate 'iarray'*/
   for(int j = 0; j < *ptr_n; j++){
      iarray[j] = j;
   }
   /*allocate copies of 'x' and 'w' (the order of the copies is modified)*/
   x_cpy = (double*) Calloc(*ptr_n, double);
   Memcpy(x_cpy, x, *ptr_n);
   w_cpy = (double*) Calloc(*ptr_n, double);
   /*sort array 'x_cpy' and return index numbers in 'iarray'*/
   rsort_with_index(x_cpy, iarray, *ptr_n);
   /*order 'w' according to 'iarray' and compute weight total*/
   for(int j = 0; j < *ptr_n; j++){
      w_cpy[iarray[j]] = w[j];
      total += w[j];
   }
   /*compute the weighted quantile*/
   total *= *ptr_prob;
   s = w_cpy[0];
   while(s < total){
      s += w_cpy[i];
      i++;
   }
   *ptr_quant = x_cpy[i - 1];
   /*housekeeping*/
   Free(w_cpy); Free(x_cpy); Free(iarray);
}
/* wmad: weighted median absolute deviation from the median
 * -------------------------------------------------------------------------
 * PARAMETERS
 *    x	       data
 *    w	       weight
 *    ptr_mad  weighted mad (on return)
 *    ptr_n    array dimension
 */
void wmad(double *x, double *w, double *ptr_mad, int *ptr_n)
{
 //  double *absdev;
   double med = 0.0, prob = 0.5;
   double *absdev;
   /*compute weighted median of x*/
   wquantile(x, w, &prob, &med, ptr_n);
   /*compute absolute deviation from the weighted median*/
   absdev = (double *) Calloc(*ptr_n, double);
   for(int i = 0; i < *ptr_n; i++){
      absdev[i] = fabs(x[i] - med);
   }
   /*weighted median of absdev*/
   wquantile(absdev, w, &prob, ptr_mad, ptr_n);
   /*housekeeping*/
   Free(absdev);
}
/* wmeantrimmed: weighted trimmed mean
 * -------------------------------------------------------------------------
 * PARAMETERS
 *    x	       vector of data
 *    w	       vector of weights
 *    LB_ptr   lower bound
 *    UB_ptr   upper bound
 *    ptr_mean estimated mean
 *    ptr_n    array dimension
 *
 * DEPENDENCIES
 *    rsort_with_index	<R_ext/Utils.h>
 */
void wmeantrimmed(double *x, double *w, double *LB_ptr, double *UB_ptr,
		  double *ptr_mean, int *ptr_n)
{
   int i = 0;
   int *iarray;
   double lb = 0.0, ub = 0.0, cumsum_x = 0.0, cumsum_w = 0.0, total = 0.0,
	  wsumtrim = 0.0;
   double *w_cpy, *x_cpy;
   if (*ptr_n == 1){
      *ptr_mean = x[0];
   }else{
      /*allocate array 'iarray'*/
      iarray = (int*) Calloc(*ptr_n, int);
      /*populate 'iarray'*/
      for(int j = 0; j < *ptr_n; j++){
	 iarray[j] = j;
      }
      /*allocate copies of 'x' and 'w' (the order of the copies is modified)*/
      x_cpy = (double*) Calloc(*ptr_n, double);
      Memcpy(x_cpy, x, *ptr_n);
      w_cpy = (double*) Calloc(*ptr_n, double);
      /*sort array 'x_cpy' and return index numbers in 'iarray'*/
      rsort_with_index(x_cpy, iarray, *ptr_n);
      /*order 'w' according to 'iarray' and compute weight total*/
      for(int j = 0; j < *ptr_n; j++){
	 w_cpy[iarray[j]] = w[j];
	 total += w[j];
      }
      /*thresholds, lower and upper bounds*/
      lb = *LB_ptr * total;
      ub = *UB_ptr * total;
      /*compute weighted trimmed mean*/
      do {
	 cumsum_w += w_cpy[i];
	 if(cumsum_w > lb){
	    cumsum_x += x_cpy[i] * w_cpy[i];
	    wsumtrim += w_cpy[i];
	 }
	 i++;
      }while(cumsum_w < ub);
      if(wsumtrim > 0.0){
	 *ptr_mean = cumsum_x / wsumtrim;
      }else{
	 *ptr_mean = 0.0;
      }
   /*housekeeping*/
   Free(w_cpy); Free(x_cpy); Free(iarray);
   }
}
/* wmeanwinsorized: weighted winsorized mean
 * -------------------------------------------------------------------------
 * PARAMETERS
 *    x	       vector of data
 *    w	       vector of weights
 *    LB_ptr   lower bound
 *    UB_ptr   upper bound
 *    ptr_mean estimated mean
 *    ptr_n    array dimension
 *
 * DEPENDENCIES
 *    rsort_with_index	<R_ext/Utils.h>
 */
void wmeanwinsorized(double *x, double *w, double *LB_ptr, double *UB_ptr,
		     double *ptr_mean, int *ptr_n)
{
   int i = 0, lastlow = 0;
   int *iarray;
   double lb = 0.0, ub = 0.0, cumsum_x = 0.0, cumsum_w = 0.0, total = 0.0,
	  wsumwin = 0.0, sumlow = 0.0;
   double *w_cpy, *x_cpy;
   if (*ptr_n == 1){
      *ptr_mean = x[0];
   }else{
      /*allocate array 'iarray'*/
      iarray = (int*) Calloc(*ptr_n, int);
      /*populate 'iarray'*/
      for(int j = 0; j < *ptr_n; j++){
	 iarray[j] = j;
      }
      /*allocate copies of 'x' and 'w' (the order of the copies is modified)*/
      x_cpy = (double*) Calloc(*ptr_n, double);
      Memcpy(x_cpy, x, *ptr_n);
      w_cpy = (double*) Calloc(*ptr_n, double);
      /*sort array 'x_cpy' and return index numbers in 'iarray'*/
      rsort_with_index(x_cpy, iarray, *ptr_n);
      /*order 'w' according to 'iarray' and compute weight total*/
      for(int j = 0; j < *ptr_n; j++){
	 w_cpy[iarray[j]] = w[j];
	 total += w[j];
      }
      /*thresholds, lower and upper bounds*/
      lb = *LB_ptr * total;
      ub = *UB_ptr * total;
      /*compute weighted winsorized mean*/
      while(cumsum_w + w_cpy[i] < ub){
	 cumsum_w += w_cpy[i];
	 if(cumsum_w > lb){
	    cumsum_x += x_cpy[i] * w_cpy[i];
	    wsumwin += w_cpy[i];
	 }else{
	    lastlow++;
	    sumlow = cumsum_w;
	 }
	 i++;
      }
      *ptr_mean = (x[lastlow] * sumlow + cumsum_x + (total - sumlow - wsumwin) *
	 x[i]) / total;
   /*housekeeping*/
   Free(w_cpy); Free(x_cpy); Free(iarray);
   }
}



