#ifndef _CSC_SPARSEMATRIX_H_
#define _CSC_SPARSEMATRIX_H_
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <optional>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include <ctime>
#include "triple_smatrix.h"
#include "basicfun.h"
class CSC_SMatrix{
    public:
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
        CSC_SMatrix(const int row, const int col, const int ne);
        CSC_SMatrix(const Triple_SMatrix& A);
        CSC_SMatrix(const CSC_SMatrix& A);
        CSC_SMatrix(CSC_SMatrix&& A);
        CSC_SMatrix(const smi& ne, smi* rows, smi* cols, double* values);
        ~CSC_SMatrix();
        
        /* Memory */
        friend void* sm_malloc(smi i, size_t size);
        friend void* sm_calloc(smi i, size_t size);
        friend void sm_free(void* p);
        friend void* sm_realloc(void* p, size_t size);
        bool sm_sprealloc(const smi& ne);
        //CSC_SMatrix* sm_smalloc(const smi& nr, const smi& nc, const smi& ne, const smi& values);
        //CSC_SMatrix* sm_resmalloc(CSC_SMatrix& A, const smi& ne);

        /* Mrithmetic */
        double* sm_gaxpy(const double* x, const double* y = nullptr)const;
        CSC_SMatrix csc_transpose();
        CSC_SMatrix csc_sort();
        smi csc_scatter(smi j, const double& beta, smi* w, double* x, smi mark, CSC_SMatrix& C, smi ne)const;
        CSC_SMatrix csc_multiply(const CSC_SMatrix& another);
        friend CSC_SMatrix csc_add(const CSC_SMatrix& A, const CSC_SMatrix& B, const double& alph, const double& beta);
        friend smi* csc_pvec(const smi* p, smi* x, smi* b, smi& n);
        friend smi* csc_ipvec(const smi* p, smi* x, smi* b, smi& n);
        friend smi* csc_pinv(const smi* p, smi n);
        CSC_SMatrix csc_permutation(const smi* pinv, const smi* q);
        CSC_SMatrix csc_sympvem(const smi* pinv)const;
        double csc_norm();
        double* csc_lsolve(double* b)const;
        double* csc_ltsolve(double* b)const;
        double* csc_usolve(double* b)const;
        double* csc_utsolve(double* b)const;
        #define csc_flip(i) ((-i)-2)
        #define csc_unflip(i) ((i<0)? csc_flip(i):(i))
        #define csc_mark(w,j) {w[j] = csc_flip(w[j]);}
        #define csc_marked(w,j) (w[j] < 0)
        smi csc_reach(CSC_SMatrix& B, smi k, smi* xi, const smi* pinv);
        smi csc_dfs(smi j, smi top, smi* x, smi* pstack, const smi* pinv);
        smi csc_spsolve(CSC_SMatrix& B, smi k, smi* xi, double* x, const smi* pinv, smi lo);

        /* Utilization */
        smi cols()const;
        smi rows()const;
        smi entries()const;
        smi& cols();
        smi& rows();
        smi& entries();
        CSC_SMatrix csc_duplicate();
        double csc_cumsum(int* p, int* c, int n);
        bool csc_sm_fkeep(bool(*fkeep)(int ,int, double, void*), void* other);
        friend bool csc_nonzero(int i, int j, double aij, void *other);
        bool csc_dropzero();
        double operator()(const smi& i, const smi& j)const;
        friend std::ostream& operator<<(std::ostream& os, const CSC_SMatrix& csc);
        friend std::fstream& output(std::fstream& os, const CSC_SMatrix& csc);
};
void* sm_malloc(smi i, size_t size);
void* sm_calloc(smi i, size_t size);
void sm_free(void* p);
void* sm_realloc(void* p, size_t size);
smi* csc_pvec(const smi* p, smi* x, smi* b, smi& n);
smi* csc_ipvec(const smi* p, smi* x, smi* b, smi& n);
#endif