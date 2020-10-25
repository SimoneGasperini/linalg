#include "linalg.hpp"

main () {

    srand((unsigned int)time(NULL));
    int r = 7, c = 5;
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
    cout << "\nSigma =\n" << USV[1];
    cout << "\nV =\n" << USV[2];
    Matrix A_again = USV[0] * USV[1] * USV[2].T();
    cout << "\nU * Sigma * V.T =\n" << A_again;

    r = 8, c = 8;
    mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 4 - 2;
        }
    }
    Matrix B (r, c, mat);
    Matrix S = B.Triu() + (B.Triu()).T() - Diag(Diag(B));
    cout << "\n\n\nSymmetric matrix S =\n" << S;

    double norm = S.Norm('f');
    cout << "Frobenius Norm = " << norm;

    double tr = S.Trace();
    cout << "\nTrace = " << tr;

    double det = S.Determinant();
    cout << "\nDeterminant = " << det << "\n";

    Matrix* QL = S.Eigendecomposition();
    cout << "\nQ =\n" << QL[0];
    cout << "\nLambda =\n" << QL[1];
    Matrix S_again = QL[0] * QL[1] * QL[0].T();
    cout << "\nQ * Lambda * Q.T =\n" << S_again;

    r = 5, c = 5;
    mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 4 - 2;
        }
    }
    mat[0][0] = 0;
    mat[3][3] = 0;
    Matrix C (r, c, mat);
    cout << "\n\n\nSquare matrix C =\n" << C;

    Matrix* dec = C.LUdecomposition();
    cout << "\nL =\n" << dec[0];
    cout << "\nU =\n" << dec[1];
    cout << "\nP =\n" << dec[2];
    Matrix C_again = dec[2].T() * dec[0] * dec[1];
    cout << "\nP.T * L * U =\n" << C_again;

}