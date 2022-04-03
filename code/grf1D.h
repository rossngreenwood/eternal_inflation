#ifndef GRF_1D
#define GRF_1D

#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>

double* gaussian_random_field_eval_c(const mxArray *prhs1, const double prhs0, const int nd)
{
    //declare variables
    mxArray *a_in_m, *nd_in_m, *V_out_m, *qq;
    const mwSize *dims;
    double *a, *b, *c, *d, *V;
    int dimx, kmax, numdims;
    int k, j;
    double q, x, sqk, V_local;

    // Associate inputs
    a_in_m = mxDuplicateArray(prhs1);

    //figure out dimensions
    dims = mxGetDimensions(prhs1);
    numdims = mxGetNumberOfDimensions(prhs1);
    kmax = (int)dims[0]-1; dimx = (int)dims[1];

    // Potential output
    V_out_m = mxCreateDoubleScalar(mxREAL);

    // Associate pointers
    a = mxGetPr(a_in_m);
    x = prhs0;

    V = mxGetPr(V_out_m);

    V_local = 0;

    sqk = sqrt(kmax);
    switch(nd) {
        case 0 :
            for(k=0;k<=kmax;k++) {
                q = k/sqk;
                V_local = V_local + a[k]*cos(q*x) + a[k+(kmax+1)]*sin(q*x); }
                break;
        case 1 :
            for(k=0;k<=kmax;k++) {
                q = k/sqk;
                V_local = V_local - a[k]*q*sin(q*x) + a[k+(kmax+1)]*q*cos(q*x); }
                break;
        case 2 :
            for(k=0;k<=kmax;k++) {
                q = k/sqk;
                V_local = V_local - a[k]*q*q*cos(q*x) - a[k+(kmax+1)]*q*q*sin(q*x); }
                break;
        case 3 :
            for(k=0;k<=kmax;k++) {
                q = k/sqk;
                V_local = V_local + a[k]*pow(q,nd)*sin(q*x) - a[k+(kmax+1)]*pow(q,nd)*cos(q*x); }
                break;
    }

    V[0] = V_local;
    return V;

}

#endif
