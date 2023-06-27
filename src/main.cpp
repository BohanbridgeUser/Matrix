#include "../include/triple_smatrix.h"
#include "../include/csc_smatrix.h"
#include <iostream>

int main(int argv, char* argc[])
{
    /* 
    Matrix 1:
            0 1 3 4 0 9
            0 2 1 0 0 0
            0 0 1 2 0 0
            0 1 0 0 0 3
            0 0 0 1 0 4
    */
   std::cout << "*************Constructor test*******************\n"; 
    std::cout << "Matrix1:\n";
    smi ent = 12;
    smi c[] = {1,1,1,2,2,2,3,3,3,5,5,5};
    smi r[] = {0,1,3,0,1,2,0,2,4,0,3,4};
    double v[] = {1,2,1,3,1,1,4,2,1,9,3,4};
    Triple_SMatrix tm(ent,r,c,v);
    std::cout << tm;
    CSC_SMatrix csc(tm);
    std::cout << csc;

    /*
    Matrix 2:
    0
    0
    1
    2
    3
    4
    */
    std::cout << "Matrix2:\n";
    ent = 4;
    smi c2[] = {0,0,0,0};
    smi r2[] = {2,3,4,5};
    double v2[] = {1,2,3,4};
    Triple_SMatrix tm2(ent,r2,c2,v2);
    std::cout << tm2;
    CSC_SMatrix csc2(tm2);
    std::cout << csc2;

    /*
    Matrix 3:
    0 1 2 3 4 5 0
    */
    std::cout << "Matrix3:\n";
    ent = 5;
    smi c3[] = {1,2,3,4,5};
    smi r3[] = {0,0,0,0,0};
    double v3[] = {1,2,3,4,5};
    Triple_SMatrix tm3(ent,r3,c3,v3);
    std::cout << tm3;
    CSC_SMatrix csc3(tm3);
    std::cout << csc3;

    /*
    Matrix 4:
    0 0 0 0 0 0
    0 0 0 1 2 3
    0 0 0 0 5 0
    0 0 0 0 8 0
    */
    std::cout << "Matrix4:\n";
    ent = 5;
    smi c4[] = {3,4,4,4,5};
    smi r4[] = {1,1,2,3,2};
    double v4[] = {1,2,3,5,8};
    Triple_SMatrix tm4(ent,r4,c4,v4);
    std::cout << tm4;
    CSC_SMatrix csc4(tm4);
    std::cout << csc4;

    /*
    Matrix 5:
    0 0 0 0 0 0
    0 1 2 3 0 3
    0 0 5 0 0 0
    0 6 0 7 0 0
    */
    std::cout << "Matrix5:\n";
    ent = 7;
    smi c5[] = {1,1,2,2,3,3,5};
    smi r5[] = {1,3,1,2,1,3,1};
    double v5[] = {1,6,2,5,3,7,3};
    Triple_SMatrix tm5(ent,r5,c5,v5);
    std::cout << tm5;
    CSC_SMatrix csc5(tm5);
    std::cout << csc5;

    /* Random accese test */
    std::cout << "*************Random access test*******************\n"; 
    std::cout << "Random accese test1:\n";
    for (int i=0;i<5;++i) {
        for (int j=0;j<6;++j) {
            std::cout << csc(i,j) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "Random accese test2:\n";
    for (int i=0;i<6;++i) {
        for (int j=0;j<1;++j) {
            std::cout << csc2(i,j) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "Random accese test3:\n";
    for (int i=0;i<1;++i) {
        for (int j=0;j<7;++j) {
            std::cout << csc3(i,j) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "Random accese test4:\n";
    for (int i=0;i<4;++i) {
        for (int j=0;j<6;++j) {
            std::cout << csc4(i,j) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << "Random accese test5:\n";
    for (int i=0;i<4;++i) {
        for (int j=0;j<6;++j) {
            std::cout << csc5(i,j) << ' ';
        }
        std::cout << std::endl;
    }

    /* Gaxpy accese test */
    std::cout << "*************Gaxpy accese test*******************\n"; 
    std::cout << "Gaxpy accese test1:\n";
    double x[] = {2,9,8,6,4,5};
    double y[] = {0,9,7,3,4};
    double* t=csc.sm_gaxpy(x,y);
    for (int i=0;i<csc.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    csc.sm_free(t);
    std::cout << "Gaxpy accese test2:\n";
    double x2[] = {10};
    double y2[] = {30,20,15,23,13,20};
    t=csc2.sm_gaxpy(x2,y2);
    for (int i=0;i<csc2.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    csc.sm_free(t);
    std::cout << "Gaxpy accese test3:\n";
    double x3[] = {10,11,18,19,20,21,22};
    double y3[] = {30};
    t=csc3.sm_gaxpy(x3,y3);
    for (int i=0;i<csc3.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    csc.sm_free(t);
    std::cout << "Gaxpy accese test4:\n";
    double x4[] = {10,11,18,19,20,21};
    double y4[] = {30,23,42,35};
    t=csc4.sm_gaxpy(x4,y4);
    for (int i=0;i<csc4.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    csc.sm_free(t);
    std::cout << "Gaxpy accese test5:\n";
    double x5[] = {10,11,18,19,20,21};
    double y5[] = {30,23,42,35};
    t=csc5.sm_gaxpy(x5,y5);
    for (int i=0;i<csc5.rows();++i) {
        std::cout << t[i] << ' ';
    }
    csc.sm_free(t);
    std::cout << std::endl;
    return 0;
}