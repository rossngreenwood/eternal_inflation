#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>
#include "grf1D.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Inputs
    mxArray *phiend_in_m, *a_in_m, *Vscale_in_m, *rho_Lambda_in_m, *phiscale_in_m;
    mxArray *lambdascreenmode_in_m, *precisionmode_in_m;
    double phiend, *a, Vscale, rho_Lambda, phiscale;
    bool lambdascreenmode, precisionmode;

    // Outputs
    mxArray *phipeak_out_m, *Vpeak_out_m;
    double phipeak, phipeak_last, *Vppeak, *phi, *phipeak_out;

    // Local variables
    double *Vend, *Vpend;
    double dphi,phimin,phitmp,phimax;
    double step;
    const mwSize *dims;
    int kmax, sgn, r;

    // Associate inputs
    phiend_in_m           = mxDuplicateArray(prhs[0]);
    a_in_m                = mxDuplicateArray(prhs[1]);
    Vscale_in_m           = mxDuplicateArray(prhs[2]);
    rho_Lambda_in_m       = mxDuplicateArray(prhs[3]);
    phiscale_in_m         = mxDuplicateArray(prhs[4]);

    // Figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    kmax = (int)dims[0]-1;

    // Potential output
    phipeak_out_m = plhs[0] = mxCreateDoubleScalar(mxREAL);
    Vpeak_out_m   = plhs[1] = mxCreateDoubleScalar(mxREAL);

    // Associate pointers
    phiend = mxGetScalar(phiend_in_m);
    a = mxGetPr(a_in_m);
    Vscale = mxGetScalar(Vscale_in_m);
    rho_Lambda = mxGetScalar(rho_Lambda_in_m);
    phiscale = mxGetScalar(phiscale_in_m);

    phipeak = mxGetScalar(phipeak_out_m);
    phipeak_out = mxGetPr(phipeak_out_m);

    Vend  = gaussian_random_field_eval_c(a_in_m,phiend,0);
    Vpend = gaussian_random_field_eval_c(a_in_m,phiend,1);

    if (*Vpend > 0) {
      sgn = 1;
    } else {
      sgn = -1;
    }
    phimin = pow(phiscale,2) * *Vpend / (*Vend-rho_Lambda) / phiscale;
    if (phimin < 0) { phimin *= -1; }
    dphi = ( (phimin > phiscale) ? phiscale : phimin);
    dphi = 0.001*sgn*((0.01 * phiscale) > dphi ? (0.01*phiscale) : dphi);

    step = 1e2;
    phipeak = phiend;
    while (1) {

        phipeak_last = phipeak;
        phipeak += (dphi*step);

        Vppeak = gaussian_random_field_eval_c(a_in_m,phipeak,1);
        if (sgn*(*Vppeak) < 0) {
            break;
        }
    }

    phimin = phipeak;
    phimax = phipeak_last;

    *Vppeak = 0;

    dphi = 1e-7;

    phipeak = 0.5*(phimin + phimax);
    Vppeak = gaussian_random_field_eval_c(a_in_m,phipeak,1);
    while (abs((phimax-phimin)/dphi) >= 1 || sgn*(*Vppeak) < 0) {
        if (sgn*(*Vppeak) > 0) {
            phimax = phipeak;
        } else {
            phimin = phipeak;
        }
        phipeak = 0.5*(phimin + phimax);
        Vppeak = gaussian_random_field_eval_c(a_in_m,phipeak,1);
    }

    *phipeak_out = phipeak;

}
