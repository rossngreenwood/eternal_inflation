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
    mxArray *phistop_out_m, *Vstop_out_m;
    double phistop, phistop_last, *Vstop, *Vpstop, *phi, *phistop_out;

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
    lambdascreenmode_in_m = mxDuplicateArray(prhs[5]);

    // Figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    kmax = (int)dims[0]-1;

    // Potential output
    phistop_out_m = plhs[0] = mxCreateDoubleScalar(mxREAL);
    Vstop_out_m   = plhs[1] = mxCreateDoubleScalar(mxREAL);

    // Associate pointers
    phiend = mxGetScalar(phiend_in_m);
    a = mxGetPr(a_in_m);
    Vscale = mxGetScalar(Vscale_in_m);
    rho_Lambda = mxGetScalar(rho_Lambda_in_m);
    phiscale = mxGetScalar(phiscale_in_m);
    lambdascreenmode = mxGetScalar(lambdascreenmode_in_m);

    phistop = mxGetScalar(phistop_out_m);
    phistop_out = mxGetPr(phistop_out_m);

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
    dphi = -0.001*sgn*((0.01 * phiscale) > dphi ? (0.01*phiscale) : dphi);

    step = 1e2;
    phistop = phiend;
    while (1) {

        phistop_last = phistop;
        phistop += (dphi*step);

        if(lambdascreenmode == 1) {
            Vstop = gaussian_random_field_eval_c(a_in_m,phistop,0);
            if (*Vstop < rho_Lambda) {
              *phistop_out = mxGetNaN();
              return;
            }
        }

        Vpstop = gaussian_random_field_eval_c(a_in_m,phistop,1);
        if (sgn*(*Vpstop) < 0) {
            break;
        }
    }

    phimin = phistop;
    phimax = phistop_last;

    phistop = 0.5*(phimin + phimax);
    while (abs((phimax-phimin)/dphi) >= 1) {
        Vpstop = gaussian_random_field_eval_c(a_in_m,phistop,1);
        if (sgn*(*Vpstop) > 0) {
            phimax = phistop;
        } else {
            phimin = phistop;
        }
        phistop = 0.5*(phimin + phimax);
    }

    *phistop_out = phistop;

}
