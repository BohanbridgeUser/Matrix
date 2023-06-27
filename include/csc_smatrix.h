#ifndef _CSC_SPARSEMATRIX_H_
#define _CSC_SPARSEMATRIX_H_
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <optional>
#include "triple_smatrix.h"
#include "basicfun.h"
class CSC_SMatrix{
    private:
        smi nentries;
        smi nrow;
        smi ncol;
        smi* pcol;
        smi* irow;
        double* value;
        bool empty()const;
    public:
        /* Basic */
        CSC_SMatrix();
        CSC_SMatrix(const Triple_SMatrix& A);
        CSC_SMatrix(const smi& ne, smi* rows, smi* cols, double* values);
        ~CSC_SMatrix();

        /* Memory */
        void* sm_malloc(smi i, size_t size)const;
        void* sm_calloc(smi i, size_t size)const;
        void sm_free(void* p)const;
        void* sm_realloc(void* p, size_t size);
        //CSC_SMatrix* sm_smalloc(const smi& nr, const smi& nc, const smi& ne, const smi& values);
        //CSC_SMatrix* sm_resmalloc(CSC_SMatrix& A, const smi& ne);

        /* Mrithmetic */
        double* sm_gaxpy(const double* x, const double* y = nullptr)const;
        // friend void scatter(const CSC_SMatrix& oriM, const smi& colj, smi* record, smi* result_col, 
        //                     const smi& mark, CSC_SMatrix& objectM, smi& rowindex);
        
        /* Utilization */
        int cols()const;
        int rows()const;
        int entries()const;
        double operator()(const smi& i, const smi& j)const;
        friend std::ostream& operator<<(std::ostream& os, const CSC_SMatrix& csc);
};

#endif