/* 

a = conjSolve(a0,emap,msk,numiter,tol);

Very little error checking, use only as directed!! 

--Ayan Chakrabarti <ayanc@eecs.harvard.edu>
*/

#include <string.h>
#include <omp.h>
#include "mex.h"

#define MIN(a,b) (a < b ? a : b)

void doAx(double * Ap,double * p,
	  double * wx,double * wy,
	  double * wd1,double * wd2,
	  double * msk,int w, int h) {
	
	int i,j,idx;
	double tmp;

	for(i=0;i<w*h;i++)
		Ap[i] = 0;

	idx = 0;
	for(j = 0; j < w;j++)
		for(i = 0; i < h;i++) {
			if(j < w-1) {
				tmp= p[idx]-p[idx+h]; tmp *= wx[idx];
				Ap[idx] += tmp; Ap[idx+h] -= tmp;
			}
			if(i < h-1) {
				tmp= p[idx]-p[idx+1]; tmp *= wy[idx];
				Ap[idx] += tmp; Ap[idx+1] -= tmp;
			}
			if(i < h-1 && j < w-1) {
				tmp= p[idx]-p[idx+h+1]; tmp *= wd1[idx];
				Ap[idx] += tmp; Ap[idx+h+1] -= tmp;
			}
			if(i < h-1 && j > 0) {
				tmp= p[idx]-p[idx+1-h]; tmp *= wd2[idx];
				Ap[idx] += tmp; Ap[idx+1-h] -= tmp;
			}
			idx++;
		}

	for(i = i;i<w*h;i++)
		Ap[i] *= msk[i];
}

void conjSolve(double *a0, double *msk,
	       int numiter,double tol,int w,int h,int n,
	       double *emap1, double *emap2, double *emap3, double *emap4) {
	double *wx, *wy, *wd1, *wd2, *r0, *p0, *Ap0;
	int i,j,idx;

	tol *= ((double)(w*h));
	
	wx = mxCalloc(w*h,sizeof(double));
	wy = mxCalloc(w*h,sizeof(double));
	wd1 = mxCalloc(w*h,sizeof(double));
	wd2 = mxCalloc(w*h,sizeof(double));

	r0 = mxCalloc(w*h*n,sizeof(double));
	p0 = mxCalloc(w*h*n,sizeof(double));
	Ap0 = mxCalloc(w*h*n,sizeof(double));

	idx = 0;
	for(j = 0; j < w;j++)
		for(i = 0; i < h;i++) {
			if(j < w-1)
				wx[idx] = MIN(emap1[idx],emap1[idx+h]);
			if(i < h-1)
				wy[idx] = MIN(emap2[idx],emap2[idx+1]);
			if(i < h-1 && j < w-1)
				wd1[idx] = MIN(emap3[idx],emap3[idx+h+1])/2.0;
			if(i < h-1 && j > 0)
				wd2[idx] = MIN(emap4[idx],emap4[idx+1-h])/2.0;
			idx++;
		}

	omp_set_num_threads(n);
#pragma omp parallel for private(idx) schedule(static)
	for(idx = 0; idx < n; idx++) {
		double alpha, beta, tmp;
		double * a, * r, * p, * Ap;
		int i,iters;

		a = &a0[w*h*idx]; r = &r0[w*h*idx];
		p = &p0[w*h*idx]; Ap = &Ap0[w*h*idx];

		doAx(r,a,wx,wy,wd1,wd2,msk,w,h);
		for(i = 0;i < w*h;i++)
			r[i] = -r[i];
		memcpy(p,r,sizeof(double)*w*h);

		for(iters = 0;iters < numiter;iters++) {
			doAx(Ap,p,wx,wy,wd1,wd2,msk,w,h);
			alpha = 0; tmp = 0; beta = 0;
			for(i = 0; i < w*h; i++) {
				tmp = tmp + r[i]*r[i];
				alpha = alpha + p[i]*Ap[i];
			}
		
			alpha = tmp / alpha;
			for(i = 0; i < w*h; i++) {
				a[i] = a[i] + alpha*p[i];
				r[i] = r[i] - alpha*Ap[i];
				beta = beta + r[i]*r[i];
			}
			if(tmp < tol)
				break;
			beta = beta/tmp;

			for(i = 0; i < w*h;i++) {
				p[i] = r[i] + beta*p[i];
			}
		}

		/* Clip just in case */
		for(i = 0; i < w*h; i++) {
			if(a[i] < 0.0) a[i] = 0.0;
			if(a[i] > 1.0) a[i] = 1.0;
		}
	}

	mxFree(wx); mxFree(wy); mxFree(r0); mxFree(p0); mxFree(Ap0);
}


/* a = conjSolve(a0,emap,msk,numiter,tol); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	double *a, *emap1, *emap2, *emap3, *emap4, *msk;
	double tol;
	int w,h,n, numiter;

	const mwSize * dims;
	mwSize odims[3];

	dims = mxGetDimensions(prhs[0]);
	h = dims[0]; w = dims[1]; n = dims[2];
	odims[0] = h; odims[1] = w; odims[2] = n;

	plhs[0] = mxCreateNumericArray(3, odims, mxDOUBLE_CLASS, mxREAL);
	a = mxGetPr(plhs[0]);
	memcpy(a,mxGetPr(prhs[0]),w*h*n*sizeof(double));
	
	emap1 = mxGetPr(mxGetCell(prhs[1],0));
	emap2 = mxGetPr(mxGetCell(prhs[1],1));
	emap3 = mxGetPr(mxGetCell(prhs[1],2));
	emap4 = mxGetPr(mxGetCell(prhs[1],3));

	msk = mxGetPr(prhs[2]);

	numiter = (int) mxGetScalar(prhs[3]);
	tol = mxGetScalar(prhs[4]);
	conjSolve(a,msk,numiter,tol,w,h,n,emap1,emap2,emap3,emap4);
}
