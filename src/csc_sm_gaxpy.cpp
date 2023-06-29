#include "../include/csc_smatrix.h"

double* CSC_SMatrix::sm_gaxpy(const double* x, const double* y)const
{
    // Return a vector that Ax + y;
    if (empty()) {
        std::cerr << "Empty Matrix!\n";
        return nullptr;
    }else if (x == nullptr){
        std::cerr << "Empty Vector!\n";
        return nullptr;
    }else if (y == nullptr){
        double* ret = (double*)(sm_calloc(nrow,sizeof(double)));
        for (int j=0;j<ncol;++j) {
            for (int p=pcol[j];p<pcol[j+1];++p) {
                ret[irow[p]] += (value[irow[p]] * x[j]);
            }
        }
        return ret;
    }else {
        double* ret = (double*)(sm_calloc(nrow,sizeof(double)));
        for (int k=0;k<nrow;++k) ret[k] = y[k];
        for (int j=0;j<ncol;++j) {
            for (int p=pcol[j];p<pcol[j+1];++p) {
                ret[irow[p]] += (value[irow[p]] * x[j] + y[irow[p]]);
            }
        }
        return ret;
    }
}