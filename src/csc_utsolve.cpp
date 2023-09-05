#include "../include/csc_smatrix.h"

double* CSC_SMatrix::csc_utsolve(double* b)const
{
    if (empty() || b==NULL) return NULL;
    double* x = (double*)sm_calloc(ncol,sizeof(double));
    for (smi j=0;j<ncol;++j) x[j] = b[j];
    for (smi j=0;j<ncol;++j) {
        for (smi i=pcol[j];i<pcol[j+1]-1;++j) {
            x[j] -= x[irow[i]] * value[irow[i]];
        }
        x[j] /= value[pcol[j+1]-1];
    }
    return x;
}