#include "../include/csc_smatrix.h"

// @ sm_gaxy: return Ax + y where x and y are dense vector

double* CSC_SMatrix::sm_gaxpy(const double* x, const double* y)const
{
    // Return a vector that Ax + y;
    if (empty()) {
        std::cerr << "Empty Matrix!\n";
        return nullptr;
    }
    if (x == nullptr){
        std::cerr << "Empty Vector!\n";
        return nullptr;
    } 
    __restrict_arr double* ret = (double*)(sm_calloc(nrow,sizeof(double)));
    if (y!=nullptr) {
        for (smi i=0;i<nrow;++i) ret[i] = y[i];
    } 

    smi num_cores =omp_get_num_procs(), num_threads;
    // omp_set_num_threads(4);
    #pragma omp parallel for 
    for (smi k=0;k<nrow;++k) ret[k] = y[k];
    
    #pragma omp parallel
    {
        smi count = 0;
        #pragma omp for reduction(+ : ret[:nrow])
        for (smi j=0;j<ncol;++j) { 
            // std::cout << "Thread_Num:" << omp_get_num_threads() << std::endl;
            #pragma omp simd
            for (int p=pcol[j];p<pcol[j+1];++p) {
                ret[irow[p]] += (value[p] * x[j]);
            }
        }
    }
    return ret;
}