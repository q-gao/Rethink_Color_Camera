/* 

y = spConv(x,k,f);

% Matlab equivalent
function y = spConv(x,k,f)

y = zeros(size(x)-size(k)+1);
i = [1:f:size(y,1) 2:f:size(y,1)];
j = [1:f:size(y,2) 2:f:size(y,2)];
for m = 1:size(k,1)
  for n = 1:size(k,2)
    y(i,j) = y(i,j) + k(end-m+1,end-n+1)*x(i+m-1,j+n-1);
  end;
end;

Very little error checking, use only as directed!! 

--Ayan Chakrabarti <ayanc@eecs.harvard.edu>
*/

#include "mex.h"
#include <omp.h>

void conv2(double * y,double * x,double *k,
	   int yi, int yj, int ksz, int f) {

	int xi,jj;
	xi = yi+ksz-1;
	
	omp_set_num_threads(3);
#pragma omp parallel for private(jj) schedule(static)
	for(jj = 0; jj < yj; jj+=f) {
		int i,j,m,n,id;
		for(j = jj; j < jj+2 && j < yj;j++) {
			id = 0;
			for(i = 0; i < yi; i+=id) {
				for(m = 0; m < ksz; m++)
					for(n = 0; n < ksz; n++) {
						y[yi*j+i] += k[(ksz-1-n)*ksz+ksz-1-m]
							* x[(j+n)*xi+i+m];
					}
				id = (id == 1) ? (f-1) : 1;
			}
		}
	}
}



/* y = spConv(x,k,f); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	double * y, * x, * k;
	int f, yi, yj, ksz;

	const mwSize * dims;
	mwSize odims[2];

	dims = mxGetDimensions(prhs[1]);
	ksz = dims[0];
	dims = mxGetDimensions(prhs[0]);
	yi = dims[0]-ksz+1; yj = dims[1]-ksz+1;

	x = mxGetPr(prhs[0]);
	k = mxGetPr(prhs[1]);
	f = (int) mxGetScalar(prhs[2]);

	odims[0] = yi; odims[1] = yj;

	plhs[0] = mxCreateNumericArray(2, odims, mxDOUBLE_CLASS, mxREAL);
	y = mxGetPr(plhs[0]);
	
	conv2(y,x,k,yi,yj,ksz,f);
}
