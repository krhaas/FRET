/*
 *Testing such that 
 *[Anorm, alpha] = runAlpha(atemp,DTs,iTad,PA,PD,LQ)
 *It also modifies atemp for final time step output
 */

#if defined(_WIN32)
#define dsymv_ dsymv
#define dgemv_ dgemv
#endif

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "blas.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double *atemp = mxGetPr(prhs[0]);
    double *DTs = mxGetPr(prhs[1]);
    bool *iTad = mxGetLogicals(prhs[2]);
    double *PA = mxGetPr(prhs[3]);
    double *PD = mxGetPr(prhs[4]);
    double *LQ = mxGetPr(prhs[5]);
    ptrdiff_t Nv = mxGetM(prhs[0]);
    int Times = mxGetM(prhs[1]);
    
    double *Anorm;
    plhs[0] = mxCreateDoubleMatrix(1, Times+2, mxREAL);
    Anorm = mxGetPr(plhs[0]);
    Anorm[0]=1;
    
    double *alpha;
    plhs[1] = mxCreateDoubleMatrix(Nv, Times+1, mxREAL);
    alpha = mxGetPr(plhs[1]);
    
    char *chn = "N";
    char *chu = "U";
    
    char buffer[50];
            
    double zero = 0.0;
    double one = 1.0;
    long ione = 1;
    int t,i,j,tp1;
    double dt;
    double tempA[Nv];
    
    t=0;
    Anorm[t] = sqrt(ddot_( &Nv, atemp, &ione, atemp, &ione ));
    
    for( i=0; i<Nv; i++) {
        atemp[i] = atemp[i]/Anorm[t];
        alpha[t*Nv+i] = atemp[i];
    }
    
    for( t=0; t < Times; t++ )
    {   
        /* Get Time Step*/
        /*dt = Tad[t]-Tad[t-1];*/
        
        /* Time propagate with eigenvalue*/
        for( i=0; i<Nv; i++)  
            atemp[i] *= exp(-LQ[i]*DTs[t]);
        
        /* Apply photon matrix*/
        if( iTad[t] ) 
            dsymv_(chu, &Nv, &one, PD, &Nv, atemp, &ione, &zero, tempA, &ione);
                
        else 
            dsymv_(chu, &Nv, &one, PA, &Nv, atemp, &ione, &zero, tempA, &ione);
        
        tp1=t+1;
        /* Get Normalization*/
        Anorm[tp1] = sqrt(ddot_( &Nv, tempA, &ione, tempA, &ione ));
        
        /* Normalize and save alpha state*/
        for( i=0; i<Nv; i++) {
            atemp[i] = tempA[i]/Anorm[tp1];
            alpha[tp1*Nv+i] = atemp[i];
        }
        
    }
}
