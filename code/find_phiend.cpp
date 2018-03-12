#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>

double* gaussian_random_field_eval_c(const mxArray *prhs1, const double prhs0, const int nd)//const mxArray *prhs[], const int nd)
{
    //declare variables
    mxArray *a_in_m, *nd_in_m, *V_out_m, *qq;
    const mwSize *dims;
    double *a, *b, *c, *d, *V;
    // const double x_in_m;
    int dimx, kmax, numdims;
    int k, j;
    double q, x, sqk, V_local;

    // Associate inputs
    a_in_m = mxDuplicateArray(prhs1);
    // x_in_m = mxDuplicateArray(prhs0);
    // x_in_m = prhs0;

    //figure out dimensions
    dims = mxGetDimensions(prhs1);
    numdims = mxGetNumberOfDimensions(prhs1);
    kmax = (int)dims[0]-1; dimx = (int)dims[1];

    // Potential output
    V_out_m = mxCreateDoubleScalar(mxREAL);

    // Associate pointers
    a = mxGetPr(a_in_m);
    // x = mxGetScalar(x_in_m);
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Inputs
    mxArray *phistart_in_m, *a_in_m, *Vscale_in_m, *rho_Lambda_in_m, *phiscale_in_m;
    mxArray *lambdascreenmode_in_m, *search_direction_in_m, *also_peak_in_m;
    double phistart, *a, Vscale, rho_Lambda, phiscale;
    bool lambdascreenmode, also_peak;
    int search_direction;

    mxArray *phipeak_out_m;
    double *phipeak_out, phipeak, phipeakmin, phipeakmax;

    // Outputs
    mxArray *phiend_out_m, *Vstop_out_m;
    double phiend, phiend_last, *Vstop, *Vpstop, *phi, *phiend_out;

    // Local variables
    double *Vstart, *Vpstart, *Vppstop, *Vend;
    double dphi,phimin,phitmp,phimax;
    double step;
    const mwSize *dims;
    int kmax, sgn, r, ind;
    bool flag_crossed_maximum;

    // Associate inputs
    phistart_in_m         = mxDuplicateArray(prhs[0]);
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
    phiend_out_m = plhs[0] = mxCreateDoubleScalar(mxREAL);
    Vstop_out_m    = mxCreateDoubleScalar(mxREAL);
    phipeak_out_m = plhs[1] = mxCreateDoubleScalar(mxREAL);

    // Associate pointers
    phistart = mxGetScalar(phistart_in_m);
    a = mxGetPr(a_in_m);
    Vscale = mxGetScalar(Vscale_in_m);
    rho_Lambda = mxGetScalar(rho_Lambda_in_m)/Vscale;
    phiscale = mxGetScalar(phiscale_in_m);
    lambdascreenmode = mxGetScalar(lambdascreenmode_in_m);
    search_direction = mxGetScalar(search_direction_in_m);
    also_peak = mxGetScalar(also_peak_in_m);

    phiend = mxGetScalar(phiend_out_m);
    phiend_out = mxGetPr(phiend_out_m);

    phipeak_out = mxGetPr(phipeak_out_m);
    *phipeak_out = mxGetNaN();

    Vstart  = gaussian_random_field_eval_c(a_in_m,phistart/phiscale,0);
    Vpstart = gaussian_random_field_eval_c(a_in_m,phistart/phiscale,1);

    // mexPrintf("V_start = %2.4f\n",*Vstart);
    // mexPrintf("Vp_start = %2.4f\n",*Vpstart);

    // mexPrintf("rho_Lambda = %2.4f\n",rho_Lambda);

    sgn = search_direction;
    phimin = pow(phiscale,2) * *Vpstart / (*Vstart-rho_Lambda) / phiscale;
    if (phimin < 0) { phimin *= -1; }
    dphi = ( (phimin > phiscale) ? phiscale : phimin);
    dphi = 0.001*sgn*((0.01 * phiscale) > dphi ? (0.01*phiscale) : dphi);

    // mexPrintf("dphi = %2.8f\n",dphi);

    flag_crossed_maximum = 0;
    step = pow(10.0,0.0625);
    phiend = phistart;
    ind = 1;
    while (0 < 1) {

        phiend_last = phiend;
        phiend += dphi*pow(step,ind);

        Vpstop  = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,1);

        if(flag_crossed_maximum == 1) {
            Vstop = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,0);
            if (*Vstop < rho_Lambda) {
                *phiend_out = mxGetNaN();
                return;
            } else {
                if (search_direction*(*Vpstop) > 0) {
                  break;
                }
            }
        } else {
            Vppstop = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,2);
            if (search_direction*(*Vpstop) < 0 && *Vppstop < 0) {
                flag_crossed_maximum = 1;
                if (also_peak == 1) {
                    phipeakmin = phiend;
                    phipeakmax = phiend_last;
                }
            }
        }

        // if (ind == 10000) {
        //     *phiend_out = mxGetNaN();
        //     return;
        // }
        ind = ind + 1;
    }

    // Nominclature assumes search_direction = 1
    phimin = phiend;
    phimax = phiend_last;

    phiend = 0.5*(phimin + phimax);
    while (abs((phimax-phimin)/dphi) >= 1) {
        Vpstop = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,1);
        if (search_direction*(*Vpstop) < 0) {
            phimax = phiend;
        } else {
            phimin = phiend;
        }
        phiend = 0.5*(phimin + phimax);
    }

    *phiend_out = phiend;

    if (also_peak == 1) {
        Vend = gaussian_random_field_eval_c(a_in_m,phiend/phiscale,0);
        if (*Vend > *Vstart) {
            return;
        }
        phipeak = 0.5*(phipeakmin + phipeakmax);
        while (abs((phipeakmax-phipeakmin)/dphi) >= 1) {
            Vpstop = gaussian_random_field_eval_c(a_in_m,phipeak/phiscale,1);
            if (search_direction*(*Vpstop) > 0) {
                phipeakmax = phipeak;
            } else {
                phipeakmin = phipeak;
            }
            phipeak = 0.5*(phipeakmin + phipeakmax);
        }

        *phipeak_out = phipeak;
    }

}
