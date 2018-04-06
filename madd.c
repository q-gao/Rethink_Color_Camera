/* 

madd(wts,imref,di,dj);

Very little error checking, use only as directed!! 

--Ayan Chakrabarti <ayanc@eecs.harvard.edu>
*/

#include "mex.h"
#include <omp.h>

void dmadd(double * img, double * wts,
	   double * imref,double * twts,
	   int di, int dj, int w, int h) {

	int lx,rx;
	int c,h2;
	
	h2 = h-di;

	if(dj > 0) {
		lx = 0; rx = w-dj;
	} else {
		lx = -dj; rx = w;
	}

	omp_set_num_threads(3);

#pragma omp parallel for private(c) schedule(static)
	for(c = 0; c < 3; c++) {
		double wt;
		int m,n,cidx,id1,id2;
		cidx = c*w*h;
		for(n = lx; n < rx; n++) {
			id1 = cidx+n*h; id2 = id1+dj*h+di;
			for(m = 0; m < h-di; m++) {
				wt = wts[(n-lx)*h2+m];
				if(c == 0) {
					twts[id1] += wt; twts[id2] += wt;
				}
				img[id1] += wt*imref[id2]; img[id2] += wt*imref[id1];
				id1++; id2++;
			}
		}
	}
}


/* madd(wts,imref,di,dj); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	double *img, *imref, *wts, *twts;
	int di, dj;
	int w,h;
	const mwSize * dims;
	mwSize odims[3];

	dims = mxGetDimensions(prhs[1]);
	h = dims[0]; w = dims[1];

	wts = mxGetPr(prhs[0]);
	imref = mxGetPr(prhs[1]);

	di = (int) mxGetScalar(prhs[2]);
	dj = (int) mxGetScalar(prhs[3]);

	img = mxGetPr(mexGetVariablePtr("caller","img"));
	twts = mxGetPr(mexGetVariablePtr("caller","wts"));

	dmadd(img,wts,imref,twts,di,dj,w,h);
}
