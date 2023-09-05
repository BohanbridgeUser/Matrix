#include "../include/csc_smatrix.h"

smi CSC_SMatrix::csc_spsolve(CSC_SMatrix& B, smi k, smi* xi, double* x, const smi* pinv, smi lo)
{
    smi top = ncol;
    top = csc_reach(B,k,xi,pinv);
    for (smi i=top;i<ncol;++i) x[xi[i]] = 0;
    for (smi i=pcol[k];i<pcol[k+1];++i) x[irow[i]] = B.value[i];
    for (smi j=top;j<ncol;++j) {
        smi i = xi[j];
        smi I = pinv? pinv[j]:j;
        x[i] /= value[lo? pcol[I]:pcol[I+1]-1];
        smi p = lo? pcol[I]+1 : pcol[I];
        smi q = lo? pcol[I+1] : pcol[I+1]-1;
        for(;p<q;++p) {
            x[irow[p]] -= value[irow[p]]*x[i];
        }
    }
    return top;
}