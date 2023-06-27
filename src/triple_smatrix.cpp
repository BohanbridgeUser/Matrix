#include "../include/triple_smatrix.h"

Triple_SMatrix::Triple_SMatrix(const smi& ne, const smi* rows, const smi* cols,
                                const double* values)
    :nentries(ne),col(nullptr),row(nullptr),value(nullptr),ncol(0),nrow(0)
{   
    col = new smi[ne];
    row = new smi[ne];
    value = new double[ne];
    for (int i=0;i<ne;++i) {
        row[i] = rows[i];
        col[i] = cols[i];
        value[i] = values[i];
        (ncol < col[i])? ncol = col[i]:ncol=ncol;
        (nrow < row[i])? nrow = row[i]:nrow=nrow;
    }
    ncol++, nrow++;
}
Triple_SMatrix::~Triple_SMatrix()
{
    delete [] col;
    delete [] row;
    delete [] value;
}
std::ostream& operator<<(std::ostream& os, const Triple_SMatrix& tm)
{
    using namespace std;
    os << "Triple_Matrix:\n";
    os << "Cols = " << tm.ncol << " Row = " << tm.nrow << " Entries = " << tm.nentries << std::endl;
    os << "Col:   ";
    for (int i=0;i<tm.nentries;++i) {
        os << tm.col[i] << ' ';
    }
    os << "\nRow:   ";
    for (int i=0;i<tm.nentries;++i) {
        os << tm.row[i] << ' ';
    }
    os << "\nValue: ";
    for (int i=0;i<tm.nentries;++i) {
        os << tm.value[i] << ' ';
    }
    os << endl;
    return os;
}