#include "../include/csc_smatrix.h"

CSC_SMatrix::CSC_SMatrix():nentries(0),nrow(0),ncol(0),pcol(nullptr),irow(nullptr),value(nullptr)
{

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
CSC_SMatrix::~CSC_SMatrix()
{
    sm_free(pcol);
    sm_free(irow);
    sm_free(value);
}
void* CSC_SMatrix::sm_malloc(smi i,size_t size)
{
    return (malloc(CSC_MAX<u_int32_t>(i,1)*size));
}
void* CSC_SMatrix::sm_calloc(smi i,size_t size)
{
    return (calloc(CSC_MAX<u_int32_t>(i,1), size));
}
void CSC_SMatrix::sm_free(void* p)
{
    free(p);
}
void* CSC_SMatrix::sm_realloc(void* p, size_t size)
{
    return (realloc(p,CSC_MAX<size_t>(size,1)));
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
