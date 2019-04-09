/* include some header files*/
#include <R.h>
#include <Rmath.h>
/* the following header files come with R.h, but we list them separately for 
 * reasons of transparency*/
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

/*functions that are accessible to others*/
void rwlslm(double*, double*, double*, double*, double*, double*, int*, int*, 
   double*, double*, double*, int*, double*, int*); 
void wmad(double*, double*, double*, int*);
void wquantile(double*, double*, double*, double*, int*); 
void wmeantrimmed(double*, double*, double*, double*, double*, int*);
void wmeanwinsorized(double*, double*, double*, double*, double*, int*);
void med(double*, int*, double*);
void wtselect(double*, double*, double*, int*, double*); 
void wmed(double*, double*, int*, double*);

void sel0(double*, double*, int*, int*);

/*some definitions: weighted median*/
#define BIN_SIZE_SELECT 5 
#define SMALL_N 128
















