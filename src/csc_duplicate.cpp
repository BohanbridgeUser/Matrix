#include "../include/csc_smatrix.h"


CSC_SMatrix CSC_SMatrix::csc_duplicate()
{
    if (empty()) return CSC_SMatrix();
    smi* w = (smi*)sm_calloc(nrow,sizeof(smi));
    for (int i=0;i<nrow;++i) w[i] = -1;
    smi nz = 0, q;
    for (int j=0;j<ncol;++j) {
        q = nz;
        for (int i=pcol[j];i<pcol[j+1];++i) {
            if (w[irow[i]] >= q) {
                value[w[irow[i]]] += value[i];
            }else {
                w[irow[i]] = nz;
                irow[nz] = irow[irow[i]];
                value[nz++] = value[i];
            }
        }
        pcol[j] = nz;
    }
    pcol[ncol] = nz;
    sm_free(w);
    this->sm_sprealloc(0);
    return *this;
}
