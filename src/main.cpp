#include "../include/triple_smatrix.h"
#include "../include/csc_smatrix.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>

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
    smi r4[] = {1,1,2,3,1};
    double v4[] = {1,2,5,8,3};
    Triple_SMatrix tm4(ent,r4,c4,v4);
    std::cout << tm4;
    CSC_SMatrix csc4(tm4);
    CSC_SMatrix csc4_1(tm4);
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
    CSC_SMatrix csc5_1(tm5);
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

    /* Gaxpy test */
    std::cout << "*************Gaxpy test*******************\n"; 
    std::cout << "Gaxpy test1:\n";
    double x[] = {2,9,8,6,4,5};
    double y[] = {0,9,7,3,4};
    double* t=csc.sm_gaxpy(x,y);
    for (int i=0;i<csc.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    sm_free(t);
    std::cout << "Gaxpy test2:\n";
    double x2[] = {10};
    double y2[] = {30,20,15,23,13,20};
    t=csc2.sm_gaxpy(x2,y2);
    for (int i=0;i<csc2.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    sm_free(t);
    std::cout << "Gaxpy test3:\n";
    double x3[] = {10,11,18,19,20,21,22};
    double y3[] = {30};
    t=csc3.sm_gaxpy(x3,y3);
    for (int i=0;i<csc3.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    sm_free(t);
    std::cout << "Gaxpy test4:\n";
    double x4[] = {10,11,18,19,20,21};
    double y4[] = {30,23,42,35};
    t=csc4.sm_gaxpy(x4,y4);
    for (int i=0;i<csc4.rows();++i) {
        std::cout << t[i] << ' ';
    }
    std::cout << std::endl;
    sm_free(t);
    std::cout << "Gaxpy test5:\n";
    double x5[] = {10,11,18,19,20,21};
    double y5[] = {30,23,42,35};
    t=csc5.sm_gaxpy(x5,y5);
    for (int i=0;i<csc5.rows();++i) {
        std::cout << t[i] << ' ';
    }
    sm_free(t);
    std::cout << std::endl;

    /* Transpose test */
    std::cout << "*************Transpose test*******************\n"; 
    std::cout << "Transpose test1:\n";
    CSC_SMatrix csct = csc.csc_transpose();
    std::cout << csct << std::endl;

    std::cout << "Transpose test2:\n";
    CSC_SMatrix csct2 = csc2.csc_transpose();
    std::cout << csct2 << std::endl;

    std::cout << "Transpose test3:\n";
    CSC_SMatrix csct3 = csc3.csc_transpose();
    std::cout << csct3 << std::endl;

    std::cout << "Transpose test4:\n";
    CSC_SMatrix csct4 = csc4.csc_transpose();
    std::cout << csct4 << std::endl;

    std::cout << "Transpose test5:\n";
    CSC_SMatrix csct5 = csc5.csc_transpose();
    std::cout << csct5 << std::endl;

    std::cout << "Sort test6:\n";
    CSC_SMatrix csct6 = csc5.csc_sort();
    std::cout << csct6 << std::endl;

    /* Multiply test */
    std::cout << "*************Multiply test*******************\n"; 
    std::cout << "Multiply test1:\n";
    CSC_SMatrix csct7 = csc.csc_multiply(csc2);
    std::cout << csct7 << std::endl;

    std::cout << "Multiply test2:\n";
    std::cout << "Matrix6:\n";
    ent = 12;
    smi c7[] = {1,1,1,2,2,2,2,3,3,4,4,4};
    smi r7[] = {0,2,4,0,1,3,4,2,3,0,2,3};
    double v7[] = {1,3,5,2,2,4,4,1,2,1,1,2};
    Triple_SMatrix tm7(ent,r7,c7,v7);
    std::cout << tm7;
    CSC_SMatrix csc7(tm7);
    std::cout << csc7;

    std::cout << "Matrix7:\n";
    ent = 7;
    smi c8[] = {0,0,0,0,2,2,2};
    smi r8[] = {0,1,3,4,0,1,3};
    double v8[] = {1,3,5,7,7,9,5};
    Triple_SMatrix tm8(ent,r8,c8,v8);
    std::cout << tm8;
    CSC_SMatrix csc8(tm8);
    std::cout << csc8;

    std::cout << "M7*M8:\n";
    CSC_SMatrix cscm = csc7.csc_multiply(csc8);
    std::cout << cscm << std::endl;
    std::cout << cscm.csc_sort() << std::endl;

    /* Add test */
    std::cout << "*************Add test*******************\n"; 
    std::cout << "Add test1:\n";
    std::cout <<"csc4_1:\n" << csc4_1 << std::endl;
    std::cout << csc5_1 << std::endl; 
    CSC_SMatrix addm = csc_add(csc4_1,csc5_1,1,1);
    std::cout << addm << std::endl;
    std::cout << addm.csc_sort() << std::endl;
    
    // /* 1138_bus matrix test 4054 nonzeros */
    // using namespace std;
    // std::cout << "*************1138_bus matrix test*******************\n"; 
    // std::cout << "1138_bus matrix:\n";
    // std::fstream file;
    // file.open("G:/COMPUTER/VS_Matrix/Matrix/example/1138_bus/1138_bus.mtx", std::ios::in);
    // if (file.is_open()) std::cout << "file opened!\n";
    // smi rn = 0;
    // char str[100];
    // while (rn != 13 && file.getline(str, 100)) {
    //     rn++;
    //     cout << str << endl;
    // }
    // std::cout << str << std::endl;

    // smi rows, cols, values;
    // file >> rows >> cols >> values;
    // smi* c11 = (smi*)calloc(values, sizeof(smi));
    // smi* r11 = (smi*)calloc(values, sizeof(smi));
    // double* v11 = (double*)calloc(values, sizeof(double));
    // for (smi i = 0; i < values; ++i) {
    //     file >> r11[i];
    //     file >> c11[i];
    //     file >> v11[i];
    // }
    // /*for (smi i = 0; i < values; ++i) {
    //     cout << c11[i]-- << ' ';
    //     cout << r11[i]-- << ' ';
    //     cout << v11[i] << endl;
    // }*/
    // Triple_SMatrix tm11(values, r11, c11, v11);
    // //std::cout << tm11;
    // CSC_SMatrix csc11(tm11);
    // //std::cout << csc11;

    // // fstream ofile;
    // // ofile.open("cpp1138.txt",ios::out);
    // // if (ofile.is_open()) cout << 1 << endl;
    // // output(ofile,csc11);
    // double* x12 = (double*)calloc(rows,sizeof(double));
    // double* y12 = (double*)calloc(rows,sizeof(double));

    // #pragma omp parallel for
    // for (smi i=0;i<rows;++i) x12 [i] = 2, y12[i] = 3;

    // clock_t tbegin = 0,tend = 0;
    // tbegin = clock();
    // printf("%lfd\n",tbegin);
    // double* r12;
    // for (smi i=0;i<1e3;++i) 
    //     r12 = csc11.sm_gaxpy(x12,y12);
    // tend = clock();

    // //for (smi i=0;i<rows;++i) std::cout << r12[i] <<' ';
    // std::cout << "\ntbegin :" << tbegin  <<std::endl;
    // std::cout << "tend :" << tend  << std::endl;
    // std::cout <<  tend - tbegin << std::endl;


    /* mawi_201512020330 matrix test 4054 nonzeros */
    // @ Input Matrix{
    using namespace std;
    std::cout << "*************Goodwin_030 matrix test*******************\n"; 
    std::cout << "Goodwin_030 matrix:\n";
    std::fstream file;
    file.open("G:/COMPUTER/Matrix/example_matrix/Goodwin_030/Goodwin_030.mtx", std::ios::in);
    if (file.is_open()) std::cout << "file opened!\n";
    else std::cout << "Failed Open file\n";
    smi rn = 0;
    char str[100];
    while (rn != 17 && file.getline(str, 100)) {
        rn++;
        cout << str << endl;
    }
    std::cout << str << std::endl;
    smi rows, cols, values;
    file >> rows >> cols >> values;
    smi* c11 = (smi*)calloc(values, sizeof(smi));
    smi* r11 = (smi*)calloc(values, sizeof(smi));
    double* v11 = (double*)calloc(values, sizeof(double));
    for (smi i = 0; i < values; ++i) {
        file >> r11[i];
        file >> c11[i];
        file >> v11[i];
        r11[i]--;
        c11[i]--;
    }
    Triple_SMatrix tm11(values, r11, c11, v11);
    std::cout << tm11;
    CSC_SMatrix csc11(tm11);
    std::cout << csc11;
    // }

    // @ Output Matrix Check Sparse Matrix { 
    fstream ofile1;
    ofile1.open("../matlab/Goodwin_030cpp.txt",ios::out);
    if (ofile1.is_open()) cout << 1 << endl;
    output(ofile1,csc11);
    ofile1.close();
    //}

    // @ Algebra Operation { 
    // @ Initialization
    double* x12 = (double*)calloc(rows,sizeof(double));
    double* y12 = (double*)calloc(rows,sizeof(double));
    #pragma omp parallel for
    for (smi i=0;i<rows;++i) x12 [i] = 2, y12[i] = 3;

    // @ Operation
    clock_t tbegin = 0,tend = 0;
    tbegin = clock();
    double* r12;
    for (smi i=0;i<1e4;++i) 
        r12 = csc11.sm_gaxpy(x12,y12); 
     
    tend = clock();
    std::cout << "\ntbegin :" << tbegin  <<std::endl;
    std::cout << "tend :" << tend  << std::endl;
    std::cout <<  tend - tbegin << std::endl;

    // Output Result
    fstream ofile;
    ofile.open("../matlab/cppGoodwin_030.txt",ios::out);
    if (ofile.is_open()) cout << "File opened!\n";
    else cout << "Failed open file!\n";
    for (smi i=0;i<rows;++i) ofile << fixed << std::setprecision(8) << r12[i] << endl;
    // }
    return 0;
}