#include <iostream>
#include "linalg.hpp"
using namespace std;

main () {

    srand((unsigned int)time(NULL));
    int r = 12, c = 6;
    double** mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 4 - 2;
        }
    }
    Matrix rect (r, c, mat);
    cout << "\nRectangular matrix generated R =\n" << rect;
    cout << "\nRank(R) = " << rect.Rank();
    double* svalues = rect.SingularValues();
    int n = rect.GetRows() < rect.GetCols() ? rect.GetRows() : rect.GetCols();
    Vector singvalues(n, svalues);
    cout << "\nSingular values = " << singvalues << "\n\n";

    r = 5; c = 5;
    mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 4 - 2;
        }
    }
    Matrix A (r, c, mat);
    cout << "\nSquare matrix generated A =\n" << A;
    cout << "\nRank(A) = " << A.Rank();

    Matrix G = A.Gauss();
    cout << "\n\nUpper triangular matrix (Gauss algorithm) G =\n" << G;

    Matrix I = A.I();
    cout << "\n\nInverse matrix I =\n" << I;

    Matrix triu = A.Triu();
    Matrix triuT = triu.T();
    Matrix diag = A.Diag();
    Matrix S = triu + triuT - diag;
    cout << "\n\nSymmetric matrix (from the initial one) S =\n" << S;

    double tr = S.Trace();
    cout << "\nTrace = " << tr;

    double det = S.Determinant();
    cout << "\nDeterminant = " << det;

    double* ev = S.Eigenvalues();
    Vector evals(S.GetCols(), ev);
    cout << "\nEigenvalues = " << evals;

    double* sv = S.SingularValues();
    Vector singvals(S.GetCols(), sv);
    cout << "\nSingular values = " << singvals << "\n";

    Polynomial charpol = S.CharacteristicPol ();
    cout << "Pol(S) = " << charpol << "\n\n";

}