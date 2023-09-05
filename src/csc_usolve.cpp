#include "../include/csc_smatrix.h"

double* CSC_SMatrix::csc_usolve(double* b)const
{
    if (empty() || b==NULL) return NULL;
    double* x = (double*)sm_calloc(ncol,sizeof(double));
    for (smi j=0;j<ncol;++j) x[j] = b[j];
    for (smi j=ncol-1;j>=0;--j) {
        x[j] /= value[pcol[j+1]-1];
        for (smi i=pcol[j];i<pcol[j+1]-1;++i) {
            x[irow[i]] -= x[j] * value[irow[i]];
        }
    }
    return x;
}