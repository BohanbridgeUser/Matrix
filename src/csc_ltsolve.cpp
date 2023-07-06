#include "../include/csc_smatrix.h"

double* CSC_SMatrix::csc_ltsolve(double* b)const
{
    if(empty()) return NULL;
    double* x = (double*)sm_calloc(ncol,sizeof(double));
    for (smi j=0;j<ncol;++j) x[j] = b[j];
    for (smi j=ncol-1;j>=0;--j) {
        for (smi i=pcol[j]+1;i<pcol[j+1];i++) {
            x[j] -= value[i] * x[irow[i]]; 
        }
        x[j] /= value[pcol[j]];
    }
    return x;
}