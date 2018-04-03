

#include "mex.h"


/* forward-backward, posterior probability */
void cBackward(double *prio, double *tran, double *emis, double *obsv, double *stat, double *outMatrix, mwSize nst, mwSize nob, mwSize nem)
{
    mwSize i;
    mwSize j;
    mwSize k;
    double sum;

    
    // first dimension first

    sum = 0;
    for (j=0; j<nst; j++)
    {
        outMatrix[nob-1+j*nob] = 1;
        sum = sum + outMatrix[nob-1+j*nob];
    }
     for (j=0; j<nst; j++)
    {
        outMatrix[nob-1+j*nob] = outMatrix[nob-1+j*nob]/sum;
    }


    
    for (i=nob-2; i>=0; i--) 
    {
        sum = 0;
        for (j=0; j<nst; j++)
        {
            outMatrix[i+j*nob] = 0;
            for (k=0; k<nst; k++)
            {
                outMatrix[i+j*nob] = outMatrix[i+j*nob]+outMatrix[(i+1)+k*nob]*tran[i+1+j*nob+k*nob*nst]* emis[i+1+k*nob+((int)obsv[i+1]-1)*nob*nst];
            }
            sum = sum + outMatrix[i+j*nob];
        }
        for (j=0; j<nst; j++)
        {
            outMatrix[i+j*nob] = outMatrix[i+j*nob] / sum;
        }        
    }


}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *prio;
    double *tran;
    double *emis;
    double *obsv;
    double *stat;
    mwSize nob;
    mwSize nst;                   /* size of matrix */
    mwSize nem;
    double *outMatrix;              /* output matrix */
    double *error;

    /* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("", "Five inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("", "Two outputs required.");
    }


    /* create a pointer to the real data in the input matrix  */
    prio = mxGetPr(prhs[0]);
    tran = mxGetPr(prhs[1]);
    emis = mxGetPr(prhs[2]);
    obsv = mxGetPr(prhs[3]);
    stat = mxGetPr(prhs[4]);

    /* get dimensions of the input matrix */
    nob = mxGetN(prhs[3]);
    nst = mxGetN(prhs[4]);
    nem = mxGetN(prhs[2])/nst;
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(nob,nst,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
    error = mxGetPr(plhs[1]);
    error[0] = 0;
    
  
    /* call the computational routine */
    cBackward(prio, tran, emis, obsv, stat, outMatrix, nst, nob, nem);
    
}






