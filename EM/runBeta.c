/*
 *Testing such that 
 *[Gamma] = runAlpha(btemp,Tad,iTad,PA,PD,LQ,alpha,Anorm,iLQ)
 *It also modifies btemp for final time step output
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
    
    double *btemp = mxGetPr(prhs[0]);
    double *DTs = mxGetPr(prhs[1]);
    bool *iTad = mxGetLogicals(prhs[2]);
    double *PA = mxGetPr(prhs[3]);
    double *PD = mxGetPr(prhs[4]);
    double *LQ = mxGetPr(prhs[5]);
    double *alpha = mxGetPr(prhs[6]);
    double *Anorm = mxGetPr(prhs[7]);
    double *iLQ = mxGetPr(prhs[8]);
    /*size_t Nv = mxGetM(prhs[0]);*/
    ptrdiff_t Nv = mxGetM(prhs[0]);
    int Times = mxGetM(prhs[1]);
    
    double *Gamma;
    plhs[0] = mxCreateDoubleMatrix(Nv, Nv, mxREAL);
    Gamma = mxGetPr(plhs[0]);
    
    
    char *chn = "N";
    char *chu = "U";
    
    char buffer[50];
            
    double zero = 0.0;
    double one = 1.0;
    long ione = 1;
    int t,i,j,iNv;
    double dt;
    double tempEXP[Nv];
    double tempB[Nv];
    double talpha,tEXP;
    
    for( t=Times-1; t >= 0; t-- )
    {   
        /* Get Time Step*/
        dt = DTs[t];
        
        /* Time propagate with eigenvalue*/
        for( i=0; i<Nv; i++)  
            tempEXP[i] = exp(-LQ[i]*dt);
        
        /* Apply photon matrix*/
        if( iTad[t] ) 
            dsymv_(chu, &Nv, &one, PD, &Nv, btemp, &ione, &zero, tempB, &ione);
                
        else 
            dsymv_(chu, &Nv, &one, PA, &Nv, btemp, &ione, &zero, tempB, &ione);
        
        
        /* Get Gamma transfer matrix*/
        for( i=0; i<Nv; i++) {
            talpha = alpha[t*Nv+i];
            iNv=i*Nv;
            tEXP=tempEXP[i];
            for( j=0; j<Nv; j++) {
                if( i == j)
                    Gamma[iNv+j] += talpha*tempB[j]*dt*tempEXP[j];
                else
                    Gamma[iNv+j] += talpha*tempB[j]*(tEXP-tempEXP[j])*iLQ[iNv+j];
            }
        }
        
        /* Time propagate and Normalize*/
        for( i=0; i<Nv; i++) 
            btemp[i] = tempB[i]*tempEXP[i]/Anorm[t];
        
    }
}