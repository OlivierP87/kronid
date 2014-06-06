#include "mex.h"
#include <stdlib.h>

/* prototype */
void vecPsvec(unsigned int n, int *I, int *J, double *K);

/* Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *K,*Ip,*Jp,*Kp;
	int	*I, *J;
    int n,ir,ic,sa1,sa2,mi,nzmax,nrc,ncc;
    /* Check Input-Output */
    if (nrhs!=1)
    {
        mexErrMsgTxt("one input required: n, size of the matrix");
    }
    else if (nlhs>3)
    {
        mexErrMsgTxt("Too many outputs arguments.");
    }
	/* read input data */
    n=(int)mxGetScalar(prhs[0]);
   
   	/* Create output sparse matrix */
	
	/* compute size of the output matrix  */
	nrc=n*n;
    ncc=n*(n+1)/2;

	/* allocation size  */
	nzmax = n*n;
    plhs[0]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
	
    Ip=mxGetPr(plhs[0]);
    Jp=mxGetPr(plhs[1]);
    Kp=mxGetPr(plhs[2]);
	
	/* pointer allocation */
	I = malloc(nzmax*sizeof(int));
	J = malloc(nzmax*sizeof(int));
	K = malloc(nzmax*sizeof(double));

	vecPsvec(n,I,J,K);

	for (ir=0;ir<nzmax;ir++)
	{
		Ip[ir] = (double) I[ir]+1;
		Jp[ir] = (double) J[ir]+1;
		Kp[ir] = K[ir];
	}

	/* free memory */
	free(I);
	free(J);
	free(K);
  }
