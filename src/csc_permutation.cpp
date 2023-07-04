#include "../include/csc_smatrix.h"

CSC_SMatrix CSC_SMatrix::csc_permutation(const smi* pinv, const smi* q)
{
    if (empty() || pinv==NULL || q==NULL) {
        return CSC_SMatrix();
    }
    CSC_SMatrix C(ncol,nrow,nentries);
    smi nz = 0;
    for (smi j=0;j<ncol;++j) {
        smi i = q[j];
        C.pcol[j] = nz;
        for (smi p=pcol[i];p<pcol[i+1];++p) {
            C.value[nz] = value[p];
            C.irow[nz++] = pinv[irow[p]];
        }
    }
    C.pcol[ncol] = nz;
    return C;
}