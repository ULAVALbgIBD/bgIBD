

#include "mex.h"



/* ibd inference: viterbi decoding */
void cViterbi(double *prio, double *tran, double *emis, double *obsv, double *outVit, double *outLog, double *outPre, mwSize nst, mwSize nob, mwSize nem)
{
    mwSize i;
    mwSize j;
    mwSize k;
    
    double temp;
    double max;


    for (j=0; j<nst; j++)
    {
        outLog[0+j*nob] = prio[j] + emis[0+j*nob+((int)obsv[0]-1)*nob*nst];
    }

    
    for (i=1; i<nob; i++) 
    {
        for (j=0; j<nst; j++)
        {
            max = outLog[i-1+j*nob] + tran[i+0*nob+j*nob*nst] - 1;
            outPre[i+nob*j] = 0;
            for (k=0; k<nst; k++)
            {
                temp = outLog[i-1+k*nob] + tran[i+k*nob+j*nob*nst];
                
                if( temp > max )
                {
                    max = temp;
                    outPre[i+j*nob] = k;
                    outLog[i+j*nob] = max;
                }
            }
        }
        for (j=0; j<nst; j++)
        {
            outLog[i+j*nob] = outLog[i+j*nob] + emis[i+j*nob+((int)obsv[i]-1)*nob*nst];
        }
    }
    
    max = outLog[nob-1+0*nob] - 1;
    outVit[nob-1] = 0;
    for (j=0; j<nst; j++)
    {
        if( outLog[nob-1+j*nob] > max )
        {
            max = outLog[nob-1+j*nob];
            outVit[nob-1] = j;
        }
    }
    for (i=nob-2; i>=0; i--)
    {
        outVit[i] = outPre[i+1+((int)outVit[i+1])*nob];
    }
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *prio;
    double *tran;
    double *emis;
    double *obsv;
    mwSize nob;
    mwSize nst;                   /* size of matrix */
    mwSize nem;
    mwSize i;
    mwSize j;
    mwSize k;
    
    double *outVit;              /* output matrix */
    double *outLog;
    double *outPre;
    double *error;

    /* check for proper number of arguments */
    if(nrhs!=4) 
    {
        mexErrMsgIdAndTxt("", "Four inputs required.");
    }
    if(nlhs!=4) 
    {
        mexErrMsgIdAndTxt("", "Four outputs required.");
    }


    /* create a pointer to the real data in the input matrix  */
    prio = mxGetPr(prhs[0]);
    tran = mxGetPr(prhs[1]);
    emis = mxGetPr(prhs[2]);
    obsv = mxGetPr(prhs[3]);

    /* get dimensions of the input matrix */
    nst = mxGetN(prhs[0]);
    nob = mxGetN(prhs[3]);
    nem = mxGetN(prhs[2])/nst;
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(nob, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nob, nst, mxREAL);    
    plhs[2] = mxCreateDoubleMatrix(nob, nst, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    outVit = mxGetPr(plhs[0]);
    outLog = mxGetPr(plhs[1]);
    outPre = mxGetPr(plhs[2]);
    error = mxGetPr(plhs[3]);
    error[0] = 0;   
    

    /* call the computational routine */
    cViterbi(prio, tran, emis, obsv, outVit, outLog, outPre, nst, nob, nem);
    for( i=0; i<nob; i++ )
    {
        outVit[i] = outVit[i] + 1;
    }
}






