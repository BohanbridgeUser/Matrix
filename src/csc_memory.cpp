#include "../include/csc_smatrix.h"

void* sm_malloc(smi i,size_t size)
{
    void* p = (malloc(CSC_MAX<uint32_t>(i,1)*size));
    if ( p == NULL ) {
        std::cerr << "Malloc failed!\n";
    }
    return p;
}
void* sm_calloc(smi i,size_t size)
{
    void* p = (calloc(CSC_MAX<uint32_t>(i,1), size));
    if ( p == NULL ) {
        std::cerr << "Calloc failed!\n";
    }
    return p;
}
void* sm_realloc(void* p, size_t size)
{
    p = (realloc(p,CSC_MAX<size_t>(size,1)));
    if ( p == NULL ) {
        std::cerr << "realloc failed!\n";
    }
    return p;
}
void sm_free(void* p)
{
    free(p);
}
bool CSC_SMatrix::sm_sprealloc(const smi& ne)
{
    if (ne <= 0) nentries = (empty())? 0:(pcol[ncol]);
    else nentries = ne;
    pcol = (smi*)sm_realloc(pcol,(ncol+1)*sizeof(smi));
    irow = (smi*)sm_realloc(irow,nentries*sizeof(smi));
    value = (double*)sm_realloc(value,nentries*sizeof(double));
    return true;
}