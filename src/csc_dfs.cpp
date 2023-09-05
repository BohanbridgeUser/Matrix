#include "../include/csc_smatrix.h"

smi CSC_SMatrix::csc_dfs(smi j, smi top, smi* x, smi* pstack, const smi* pinv)
{
    if (empty() || x ==NULL || pstack == NULL) return -1;
    smi head = 0;
    x[0] = j;
    while (head >=0 ) {
        j = x[head];
        smi jnew = ((pinv)? pinv[j]:j) ;
        if (!csc_marked(pcol,j)){
            csc_mark(pcol,j);
            pstack[head] = (jnew<0)? 0:csc_unflip(pcol[jnew]);
        }
        smi done = 1;
        smi p2 = (jnew<0)? 0: csc_unflip(pcol[jnew+1]);
        for (smi i=pstack[head];i<p2;++i) {
            if (!csc_marked(pcol,irow[i])) continue;
            csc_mark(pcol,irow[i]);
            x[head] = irow[i];
            pstack[++head] = i;
            done = 0;
            break;
        }
        if (done) {
            head--;
            x[--top] = irow[j];
        }
    }
    return top;
}