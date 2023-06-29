#include "../include/csc_smatrix.h"

bool CSC_SMatrix::csc_sm_fkeep(bool(*fkeep)(int, int, double, void*), void* other)
{
    smi nz = 0, i;
    if (empty()) return false;
    for (i=0;i<ncol;++i) {
        smi ri = nz;
        pcol[i] = nz;
        for (;ri<pcol[i];++ri) {
            if (fkeep(ri,i,value[irow[ri]],other)){
                value[nz] = value[ri];
                irow[nz++] = irow[ri];
            }
        }
    }
    nentries = nz;
    pcol[i] = nz;
    sm_sprealloc();
    return true;
}
bool csc_nonzero(int i, int j, double aij, void* other)
{
    return (aij != 0.0);
}
bool CSC_SMatrix::csc_dropzero()
{
    return (csc_sm_fkeep(csc_nonzero,NULL));
}