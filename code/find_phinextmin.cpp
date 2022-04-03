#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>
#include "grf1D.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Inputs
    mxArray *phiend_in_m, *a_in_m, *Vscale_in_m, *rho_Lambda_in_m, *phiscale_in_m;
    mxArray *lambdascreenmode_in_m, *search_direction_in_m, *also_peak_in_m;
    double phiend, *a, Vscale, rho_Lambda, phiscale;
    bool lambdascreenmode, also_peak;
    int search_direction;

    mxArray *phipeak_out_m;
    double *phipeak_out, phipeak, phipeakmin, phipeakmax;

    // Outputs
    mxArray *phinextmin_out_m, *Vstop_out_m;
    double phinextmin, phinextmin_last, *Vstop, *Vpstop, *phi, *phinextmin_out;

    // Local variables
    double *Vend, *Vpend, *Vppstop, *Vnextmin;
    double dphi,phimin,phitmp,phimax;
    double step;
    const mwSize *dims;
    int kmax, sgn, r, ind;
    bool flag_crossed_maximum;

    // Associate inputs
    phiend_in_m           = mxDuplicateArray(prhs[0]);
    a_in_m                = mxDuplicateArray(prhs[1]);
    Vscale_in_m           = mxDuplicateArray(prhs[2]);
    rho_Lambda_in_m       = mxDuplicateArray(prhs[3]);
    phiscale_in_m         = mxDuplicateArray(prhs[4]);
    lambdascreenmode_in_m = mxDuplicateArray(prhs[5]);
    search_direction_in_m = mxDuplicateArray(prhs[6]);
    also_peak_in_m        = mxDuplicateArray(prhs[7]);

    // Figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    kmax = (int)dims[0]-1;

    // Potential output
    phinextmin_out_m = plhs[0] = mxCreateDoubleScalar(mxREAL);
    Vstop_out_m    = mxCreateDoubleScalar(mxREAL);
    phipeak_out_m = plhs[1] = mxCreateDoubleScalar(mxREAL);

    // Associate pointers
    phiend = mxGetScalar(phiend_in_m);
    a = mxGetPr(a_in_m);
    Vscale = mxGetScalar(Vscale_in_m);
    rho_Lambda = mxGetScalar(rho_Lambda_in_m)/Vscale;
    phiscale = mxGetScalar(phiscale_in_m);
    lambdascreenmode = mxGetScalar(lambdascreenmode_in_m);
    search_direction = mxGetScalar(search_direction_in_m);
    also_peak = mxGetScalar(also_peak_in_m);

    phinextmin = mxGetScalar(phinextmin_out_m);
    phinextmin_out = mxGetPr(phinextmin_out_m);

    phipeak_out = mxGetPr(phipeak_out_m);
    *phipeak_out = mxGetNaN();

    Vend  = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,0);
    Vpend = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,1);

    sgn = search_direction;
    phimin = pow(phiscale,2) * *Vpend / (*Vend-rho_Lambda) / phiscale;
    if (phimin < 0) { phimin *= -1; }
    dphi = ( (phimin > phiscale) ? phiscale : phimin);
    dphi = 0.001*sgn*((0.01 * phiscale) > dphi ? (0.01*phiscale) : dphi);

    flag_crossed_maximum = 0;
    step = pow(10.0,0.0625);
    phinextmin = phiend;
    ind = 1;
    while (1) {

        phinextmin_last = phinextmin;
        phinextmin += dphi*pow(step,ind);

        Vpstop  = gaussian_random_field_eval_c(a_in_m,phinextmin/phiscale,1);

        if(flag_crossed_maximum == 1) {
            Vstop = gaussian_random_field_eval_c(a_in_m,phinextmin/phiscale,0);
            if (*Vstop < rho_Lambda) {
                *phinextmin_out = mxGetNaN();
                return;
            } else {
                if (search_direction*(*Vpstop) > 0) {
                  break;
                }
            }
        } else {
            Vppstop = gaussian_random_field_eval_c(a_in_m,phinextmin/phiscale,2);
            if (search_direction*(*Vpstop) < 0 && *Vppstop < 0) {
                flag_crossed_maximum = 1;
                if (also_peak == 1) {
                    phipeakmax = phinextmin;
                    phipeakmin = phinextmin_last;
                }
            }
        }

        ind = ind + 1;
    }

    // Nominclature assumes search_direction = 1
    phimin = phinextmin;
    phimax = phinextmin_last;

    phinextmin = 0.5*(phimin + phimax);
    while (abs((phimax-phimin)/dphi) >= 1) {
        Vpstop = gaussian_random_field_eval_c(a_in_m,phinextmin/phiscale,1);
        if (search_direction*(*Vpstop) < 0) {
            phimax = phinextmin;
        } else {
            phimin = phinextmin;
        }
        phinextmin = 0.5*(phimin + phimax);
    }

    *phinextmin_out = phinextmin;

    if (also_peak == 1) {
        Vnextmin = gaussian_random_field_eval_c(a_in_m,phinextmin/phiscale,0);
        if (*Vnextmin > *Vend) {
            return;
        }
        phipeak = 0.5*(phipeakmin + phipeakmax);
        Vpstop = gaussian_random_field_eval_c(a_in_m,phipeak/phiscale,1);
        while (abs((phipeakmax-phipeakmin)/dphi) >= 1 || search_direction*(*Vpstop) > 0) {
            if (search_direction*(*Vpstop) > 0) {
                phipeakmin = phipeak;
            } else {
                phipeakmax = phipeak;
            }
            phipeak = 0.5*(phipeakmin + phipeakmax);
            Vpstop = gaussian_random_field_eval_c(a_in_m,phipeak/phiscale,1);
        }

        *phipeak_out = phipeak;
    }

}
