#include "../include/csc_smatrix.h"

smi CSC_SMatrix::csc_reach(CSC_SMatrix& B, smi k, smi* x, const smi* pinv)
{
    if (empty() || B.empty() || x==NULL) return -1;
    smi top = ncol;
    for (smi j=B.pcol[k];j<B.pcol[k+1];++j) {
        if (!csc_marked(pcol,B.irow[j])){
            top = csc_dfs(B.irow[j],top,x,x+ncol,pinv);
        }
    }
    for (smi j=top;j<ncol;++j) csc_mark(pcol,j);
    return top;
}