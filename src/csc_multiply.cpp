#include "../include/csc_smatrix.h"

/* 
C = A * B
This algorithm is based on a fact that
Cij is sum of Aik * Bkj which k from 0 to m.
So, the inner cycle had calculated the jth column of C.
*/
CSC_SMatrix CSC_SMatrix::csc_multiply(const CSC_SMatrix& another)
{
    if (empty() || another.empty()) {
        std::cerr << "Empty Matrix!\n";
        return CSC_SMatrix();
    }

    CSC_SMatrix C(nrow,another.ncol,nentries+another.nentries);
    smi* w = (smi*)sm_calloc(nrow,sizeof(smi));
    double* x = (double*)sm_calloc(nrow,sizeof(double));
    smi nz = 0;
    for (smi j=0;j<another.ncol;++j) {
        if (nz + nrow > C.entries() || !C.sm_sprealloc(2*(C.entries()+nrow))){
            std::cerr << "Failed allocate memory!\n";
            exit(0);
        };
        C.pcol[j] = nz;
        for (smi i=another.pcol[j];i<another.pcol[j+1];++i) {
            nz = this->csc_scatter(another.irow[i],another.value[i],w,x,j+1,C,nz);
        }
        for (smi i=C.pcol[j];i<nz;++i) C.value[i] = x[C.irow[i]];
    }
    C.pcol[another.cols()] = nz;
    C.sm_sprealloc(0);
    sm_free(w);
    sm_free(x);
    return C;
}