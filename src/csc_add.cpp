#include "../include/csc_smatrix.h"

CSC_SMatrix csc_add(const CSC_SMatrix& A, const CSC_SMatrix& B, const double& alph, const double& beta)
{
    if (A.empty() || B.empty()) {
        std::cerr << "Emtyp Matrix!\n";
        return CSC_SMatrix();
    }
    CSC_SMatrix C(A.rows(),A.cols(),A.entries()+B.entries());
    smi *w = (smi*)sm_calloc(A.rows(),sizeof(smi));
    double *x = (double*)sm_calloc(A.rows(),sizeof(smi));
    smi nz = 0;
    std::cout << A.cols() <<std::endl;
    for (smi p=0;p<A.cols();++p) {
        C.pcol[p] = nz;
        nz = A.csc_scatter(p,alph,w,x,p+1,C,nz);
        nz = B.csc_scatter(p,beta,w,x,p+1,C,nz);
        for (smi i=C.pcol[p];i<nz;++i) C.value[i] = x[C.irow[i]];
    }
    C.pcol[A.cols()] = nz;
    std::cout << C << std::endl;
    C.sm_sprealloc(0);
    return C;
}