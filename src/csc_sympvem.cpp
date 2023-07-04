#include "../include/csc_smatrix.h"

CSC_SMatrix CSC_SMatrix::csc_sympvem(const smi* pinv)const
{
    if (empty() || pinv==NULL) return CSC_SMatrix();
    CSC_SMatrix C(ncol,nrow,nentries);
    smi* w = (smi*)sm_calloc(ncol,sizeof(smi));
    for (smi j=0;j<ncol;++j) {
        smi j2 = pinv[j];
        for (smi i=pcol[j];i<pcol[j+1];++i) {
            smi i2 = pinv[i];
            if (irow[i]>j) continue;
            w[CSC_MAX(i2,j2)]++;
        }
    }
    C.csc_cumsum(C.pcol,w,ncol);
    for (smi j=0;j<ncol;++j) {
        smi j2 = pinv[j];
        for (smi i=pcol[j];i<pcol[j+1];++i) {
            if (irow[i] > j) continue;
            smi i2 = pinv[i];
            C.irow[w[CSC_MAX(i2,j2)]] = CSC_MIN(i2,j2);
            C.value[w[CSC_MAX(i2,j2)]++] = value[i];
        }
    }
    sm_free(w);
    return C;
}