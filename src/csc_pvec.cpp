#include "../include/csc_smatrix.h"

smi* csc_pvec(const smi* p, smi* x, smi* b, const smi& n)
{
    if (!p || !x || !b) {
        return NULL;
    }
    for (smi i=0;i<n;++i) {
        x[i] = b[p[i]];
    }
    return x;
}

smi* csc_ipvec(const smi* p, smi* x, smi* b, smi& n)
{
    if (!p || !x || !b) {
        return NULL;
    }
    for (smi i=0;i<n;++i) {
        x[p[i]] = b[i];
    }
    return x; 
}

smi* csc_pinv(const smi* p, smi& n)
{
    if (!p) return NULL;
    smi* invp = (smi*)sm_calloc(n,sizeof(smi));
    for (smi k=0;k<n;++k) {
        invp[p[k]] = k;
    }
    return invp;
}