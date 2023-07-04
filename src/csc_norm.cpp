#include "../include/csc_smatrix.h"

double CSC_SMatrix::csc_norm()
{
    if (empty()) return 0;
    double norm = 0, sum = 0;
    for (smi j=0;j<ncol;++j) {
        for (smi p=pcol[j];p<pcol[j+1];++p) {
            sum += fabs(value[p]);
        }
        norm = CSC_MAX(norm,sum);
    }
    return norm;
}