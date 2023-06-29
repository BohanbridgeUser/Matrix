#include "../include/csc_smatrix.h"

CSC_SMatrix::CSC_SMatrix():nentries(0),nrow(0),ncol(0),pcol(nullptr),irow(nullptr),value(nullptr)
{

}
CSC_SMatrix::CSC_SMatrix(const int row, const int col, const int ne):nentries(ne),nrow(row),ncol(col)
{
    pcol = (smi*)sm_calloc(col+1,sizeof(smi));
    irow = (smi*)sm_calloc(ne,sizeof(smi));
    value = (double*)sm_calloc(ne,sizeof(double));
}
CSC_SMatrix::CSC_SMatrix(const Triple_SMatrix& A)
{
    /* A is arranged by col values */
    nentries = A.nentries;
    nrow = A.nrow;
    ncol = A.ncol;
    pcol = (smi*)sm_calloc(ncol+1,sizeof(smi));
    irow = (smi*)sm_calloc(nentries,sizeof(smi));
    value = (double*)sm_calloc(nentries,sizeof(double));
    int i=0,j=0,k=0;
    smi temp = 0, count = 0;
    while (i<nentries) {
        while (j<=A.col[i]) {pcol[j++] = count;}
        temp = A.col[i];
        while (temp == A.col[i] && i < nentries){
            irow[k] = A.row[i];
            value[k] = A.value[i];
            i++,k++;
            count++;
        }
        pcol[j] = count;
    }
}
CSC_SMatrix::CSC_SMatrix(const CSC_SMatrix& A):nentries(A.nentries),nrow(A.nrow),ncol(A.ncol)
{
    pcol = (smi*)sm_calloc(ncol+1,sizeof(smi));
    irow = (smi*)sm_calloc(nentries,sizeof(smi));
    value = (double*)sm_calloc(nentries,sizeof(double));
    for (int i=0;i<ncol+1;++i) pcol[i] = A.pcol[i];
    for (int i=0;i<nentries;++i) {
        irow[i] = A.irow[i];
        value[i] = A.value[i];
    }
}
CSC_SMatrix::CSC_SMatrix(CSC_SMatrix&& A):nentries(A.nentries),nrow(A.nrow),ncol(A.ncol)
{
    pcol = A.pcol;
    A.pcol = NULL;
    irow = A.irow;
    A.irow = NULL;
    value = A.value;
    A.value = NULL;
}
bool CSC_SMatrix::empty()const
{
    return (nentries == 0);
}
CSC_SMatrix::~CSC_SMatrix()
{
    sm_free(pcol);
    sm_free(irow);
    sm_free(value);
}
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
int CSC_SMatrix::cols()const
{
    return ncol;
}
int CSC_SMatrix::rows()const
{
    return nrow;
}
int CSC_SMatrix::entries()const
{
    return nentries;
}
std::ostream& operator<<(std::ostream& os, const CSC_SMatrix& csc)
{
    using namespace std;
    os << "CSC_Matrix:\n";
    os << "Pcol:  ";
    for (int i=0;i<csc.ncol+1;++i) {
        os << csc.pcol[i] << ' ';
    }
    os << "\nRow:   ";
    for (int i=0;i<csc.nentries;++i) {
        os << csc.irow[i] << ' ';
    }
    os << "\nValue: ";
    for (int i=0;i<csc.nentries;++i) {
        os << csc.value[i] << ' ';
    }
    os << endl;
    return os;
}
double CSC_SMatrix::operator()(const smi& i, const smi& j)const
{
    if (j < 0 || j > ncol || i <0 || i > nrow){
        std::cerr << "Index out of range!\n";
        exit(0);
    }
    smi count = pcol[j+1] - pcol[j];
    if (count == 0) return 0.0;
    else {
        for (int x = pcol[j];x < pcol[j+1];++x) {
            if (irow[x] == i) return value[x];
        }
        return 0.0;
    }
}
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
CSC_SMatrix CSC_SMatrix::csc_transpose()
{
    smi *Cp,*Ci,*Cx;
    if (empty()) return CSC_SMatrix();
    CSC_SMatrix ret(ncol,nrow,nentries);
    smi* p = (smi*)sm_calloc(nrow,sizeof(smi));
    for (int i=0;i<nentries;++i) p[irow[i]]++;
    csc_cumsum(ret.pcol,p,nrow);
    for (int i=0;i<ncol;++i) {
        for (int j=pcol[i];j<pcol[i+1];++j) {
            smi rrow = p[irow[j]]++;
            ret.irow[rrow] = i;
            ret.value[rrow] = value[j];
        }
    }
    sm_free(p);
    return ret;
}