#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double APPROX = 10e-9;

class Polynomial {
    int degree;
    double* coefficients;

    public:
    Polynomial () {};
    Polynomial (int);
    Polynomial (int, double*);
    Polynomial (const Polynomial&);
    ~Polynomial ();
    int GetDegree () const;
    double* GetCoefficients () const;
    void SetCoefficient (int, double);
    double Evaluate (double);
    double EvaluateDerivative (double);
    double ComputeZero ();
    double* ComputeZeros ();

    friend ostream& operator << (ostream&, Polynomial&);
    friend istream& operator >> (istream&, Polynomial&);
    Polynomial& operator = (const Polynomial&);
    Polynomial operator + (const Polynomial&);
    Polynomial operator - (const Polynomial&);
    Polynomial operator * (const Polynomial&);
    Polynomial operator / (const Polynomial&);
};

class Vector {
    int size;
    double* array;

    public:
    Vector () {};
    Vector (int);
    Vector (int, double*);
    Vector (const Vector&);
    ~Vector ();
    int GetSize () const;
    double* GetArray () const;
    double GetElement (int) const;
    void SetElement (int, double);
    double Norm (int p = 2);
    Vector Normalized (int p = 2);
    Vector ProjectedOnto (Vector);

    friend ostream& operator << (ostream&, Vector&);
    friend istream& operator >> (istream&, Vector&);
    Vector& operator = (const Vector&);
    Vector operator + (const Vector&);
    Vector operator - (const Vector&);
    double operator * (const Vector&);
    Vector operator * (double);
    Vector operator / (double);
    bool operator == (const Vector&);
    bool operator != (const Vector&);
};

class Matrix {
    int rows, cols;
    double** matrix;

    bool IsContained (int, int*, int);
    Matrix GetMatrix (int);
    Matrix GetMatrix (int*, int);
    Matrix GetCombs (int, int, int);
    Matrix Merge (const Matrix&);
	Matrix SwapRows ();
    Matrix GetPermutation (Matrix&);

    public:
    Matrix () {};
    Matrix (int, int);
    Matrix (int, int, double**);
    Matrix (const Matrix&);
    ~Matrix ();
    int GetRows () const;
    int GetCols () const;
    Vector GetRowVector (int) const;
    Vector GetColVector (int) const;
    double GetElement (int, int) const;
    void SetElement (int, int, double);
    double Norm (char p = 'f');
    Matrix T ();
    Matrix I ();
    Matrix Triu ();
    Matrix Gauss ();
    int Rank ();
    double Trace ();
    double Determinant (char method = 'g');
    Matrix* Eigendecomposition ();
    Matrix* SVdecomposition ();
    Matrix* QRdecomposition ();
    Matrix* Choleskydecomposition (bool diag = false);
    Matrix* LUdecomposition (bool diag = false);
    Matrix HouseholderReflection (int);
    Polynomial CharacteristicPol ();
    bool IsSquare ();
    bool IsSymmetric ();
    bool IsPositiveDefinite ();

    friend ostream& operator << (ostream&, Matrix&);
    friend istream& operator >> (istream&, Matrix&);
    Matrix& operator = (const Matrix&);
    Matrix operator + (const Matrix&);
    Matrix operator - (const Matrix&);
    Matrix operator * (const Matrix&);
    Matrix operator * (double);
    Matrix operator / (double);
    bool operator == (const Matrix&);
    bool operator != (const Matrix&);
};

Matrix Eye (int);
Vector Diag (Matrix);
Matrix Diag (Vector);
double Dot (Vector, Vector);
Vector Dot (Matrix, Vector);
Vector Dot (Vector, Matrix);
Matrix Dot (Matrix, Matrix);
Matrix Outer (Vector, Vector);
Vector Hadamard (Vector, Vector);
Vector* GramSchmidt (Vector*, int, int p = 2);
double WilkinsonShift (Matrix);
Vector Convolution (Vector, Vector);
Matrix Convolution (Matrix, Matrix);