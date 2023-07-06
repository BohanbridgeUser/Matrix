#include "../include/csc_smatrix.h"

double* CSC_SMatrix::csc_lsolve(double* b)const
{
    if(empty()) return NULL;
    double* x = (double*)sm_calloc(ncol,sizeof(double));
    for (smi i=0;i<ncol;++i) x[i] = b[i];
    for (smi j=0;j<ncol;++j) {
        x[j] = x[j] / value[pcol[j]];
        for (smi i=pcol[j]+1;i<pcol[j+1];++i) {
            x[irow[i]] -= value[i]*x[j];
        }
    }
    return x;
}