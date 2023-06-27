#ifndef _TRIPLE_SMATRIX_H_
#define _TRIPLE_SMATRIX_H_
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#define smi int32_t

class Triple_SMatrix{
        friend class CSC_SMatrix;
    private:
        smi nentries;
        smi nrow;
        smi ncol;
        smi* col;
        smi* row;
        double* value;
    public:
        Triple_SMatrix(const smi& ne, const smi* nr, const smi* nc, const double* values);
        virtual ~Triple_SMatrix();
        void input(const smi& r, const smi& c, const double& value);
        friend std::ostream& operator<<(std::ostream& os, const Triple_SMatrix& tm);
};

#endif