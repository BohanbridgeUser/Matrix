#include "../include/csc_smatrix.h"

void* CSC_SMatrix::sm_malloc(smi i,size_t size)const
{
    void* p = (malloc(CSC_MAX<uint32_t>(i,1)*size));
    if ( p == NULL ) {
        std::cerr << "Malloc failed!\n";
    }
    return p;
}
void* CSC_SMatrix::sm_calloc(smi i,size_t size)const
{
    void* p = (calloc(CSC_MAX<uint32_t>(i,1), size));
    if ( p == NULL ) {
        std::cerr << "Calloc failed!\n";
    }
    return p;
}
void CSC_SMatrix::sm_free(void* p)const
{
    free(p);
}
void* CSC_SMatrix::sm_realloc(void* p, size_t size)
{
    p = (realloc(p,CSC_MAX<size_t>(size,1)));
    if ( p == NULL ) {
        std::cerr << "realloc failed!\n";
    }
    return p;
}
bool CSC_SMatrix::sm_sprealloc()
{
    pcol = (smi*)sm_realloc(pcol,ncol+1);
    irow = (smi*)sm_realloc(irow,nentries);
    value = (double*)sm_realloc(value,nentries);
    return true;
}