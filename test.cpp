#include <iostream>
#include "linalg.hpp"
using namespace std;

main () {
/*
    //symmetric matrix 4x4 with det=0
    int r = 4, c = 4;
    double** mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
    }
    mat[0][0] = 1;
    mat[0][1] = -2;
    mat[0][2] = 1;
    mat[0][3] = 3;
    mat[1][0] = -2;
    mat[1][1] = 0;
    mat[1][2] = 2;
    mat[1][3] = -6;
    mat[2][0] = 1;
    mat[2][1] = 2;
    mat[2][2] = -2;
    mat[2][3] = 3;
    mat[3][0] = 3;
    mat[3][1] = -6;
    mat[3][2] = 3;
    mat[3][3] = 9;
    Matrix S (r, c, mat);
*/

    srand((unsigned int)time(NULL));
    int r = 5, c = 2;
    double** mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 4 - 2;
        }
    }
    Matrix A (r, c, mat);
    cout << "\nRectangular matrix A =\n" << A;
    Matrix* USV = A.SVdecomposition();

    cout << "\nU =\n" << USV[0];
    cout << "\nsigma =\n" << USV[1];
    cout << "\nV =\n" << USV[2];
    Matrix B = USV[0] * USV[1] * USV[2].T();
    cout << "\nB =\n" << B;

/*
    Matrix triu = A.Triu();
    Matrix triuT = triu.T();
    Matrix diag = A.Diag();
    Matrix S = triu + triuT - diag;

    cout << "\nSymmetric matrix S =\n" << S;

    //double norm = S.Norm('f');
    //cout << "\nFrobenius Norm = " << norm;

    double tr = S.Trace();
    cout << "\nTrace = " << tr;

    double det = S.Determinant();
    cout << "\nDeterminant = " << det;

    Matrix* QL = S.Eigendecomposition();
    Matrix Q = QL[0], L = QL[1];
    cout << "\n\nQ =\n" << Q;
    cout << "\nL =\n" << L;
    Matrix S_again = Q * L * Q.T();
    cout << "\nQ * L * Qt =\n" << S_again;
*/
}