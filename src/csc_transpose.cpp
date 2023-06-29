#include "../include/csc_smatrix.h"

CSC_SMatrix CSC_SMatrix::csc_transpose()
{
    smi *Cp,*Ci,*Cx;
    if (empty()) return CSC_SMatrix();
    CSC_SMatrix ret(ncol,nrow,nentries);
    smi* p = (smi*)sm_calloc(nrow,sizeof(smi));
    for (int i=0;i<nentries;++i) p[irow[i]]++;
    csc_cumsum(ret.pcol,p,nrow);
    for (int i=0;i<ncol;++i) {
        for (int j=pcol[i];j<pcol[i+1];++j) {
            smi rrow = p[irow[j]]++;
            ret.irow[rrow] = i;
            ret.value[rrow] = value[j];
        }
    }
    sm_free(p);
    return ret;
}
CSC_SMatrix CSC_SMatrix::csc_sort()
{
    if(empty()) return CSC_SMatrix();
    else {
        return (this->csc_transpose().csc_transpose());
    }
}