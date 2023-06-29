#include "../include/csc_smatrix.h"

double CSC_SMatrix::csc_cumsum(smi* p, smi* c, smi n)
{
    smi nz = 0;
    double nz2 = 0.0;
    for (int i=0;i<n;++i) {
        p[i] = nz;
        nz += c[i];
        nz2 += c[i];
        c[i] = p[i];
    }
    p[n] = nz;
    return nz2;
}