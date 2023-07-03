#include "../include/csc_smatrix.h"

/* 
w use mark to record whether each row in the vector has ever been assigned value.
For ease of distinction, mark uses the column number of ther next column.
*/
smi CSC_SMatrix::csc_scatter(smi j, const double& beta, smi* w, 
                                double* x, smi mark, CSC_SMatrix& C, smi nz)const 
{
    if (empty() || w == NULL || x == NULL) {
        std::cerr << "NULL INPUT!\n";
        return nz;
    } 
    for (smi p=pcol[j];p<pcol[j+1];++p) {
        if (w[irow[p]] < mark) {
            w[irow[p]] = mark;
            C.irow[nz++] = irow[p];
            x[irow[p]] = beta * value[p];
        } else {
            x[irow[p]] += beta * value[p];
        }
    }
    return nz;
}