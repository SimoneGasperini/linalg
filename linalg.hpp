#include <iostream>
using namespace std;

const double APPROX = 10e-9;

class Polynomial {
    int degree;
    double* coefficients;

    void Approx ();

    public:
    Polynomial () {};
    Polynomial (int);
    Polynomial (int, double*);
    Polynomial (const Polynomial&);
    int GetDegree ();
    double* GetCoefficients ();
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

class Matrix {
    int rows, cols;
    double** matrix;

    bool IsContained (int, int*, int);
    Matrix GetMatrix (int);
    Matrix GetMatrix (int*, int);
    Matrix GetCombs (int, int, int);

    Matrix Merge (const Matrix&);
	void SwapRows ();
    void Approx ();

    public:
    Matrix () {};
    Matrix (int, int);
    Matrix (int, int, double**);
    Matrix (const Matrix&);
    int GetRows ();
    int GetCols ();
    double GetElement (int, int);
    void SetElement (int, int, double);
    double Norm(char p = '2');
    Matrix T ();
    Matrix I ();
    Matrix Triu ();
    Matrix Diag ();
    Matrix Gauss ();
    int Rank();
    double Trace ();
    double Determinant (char method = 'g');
    Matrix* Eigen();
    double* SingularValues();
    Polynomial CharacteristicPol ();
    Matrix* QRdecomposition ();

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

class Vector {
    int size;
    double* array;

    void Approx ();

    public:
    Vector () {};
    Vector (int);
    Vector (int, double*);
    Vector (const Vector&);
    int GetSize();
    double* GetArray();
    double GetElement (int);
    void SetElement (int, double);
    double Norm (int p = 2);
    Vector Normalized (int p = 2);
    Vector ProjectedOnto (Vector);

    friend ostream& operator << (ostream&, Vector&);
    friend istream& operator >> (istream&, Vector&);
    Vector& operator = (const Vector&);
    Vector operator + (const Vector&);
    Vector operator - (const Vector&);
    Vector operator * (double);
    Vector operator / (double);
    bool operator == (const Vector&);
    bool operator != (const Vector&);
};

Matrix Eye(int);
double Dot (Vector, Vector);
Vector Dot (Matrix, Vector);
Vector Dot (Vector, Matrix);
Matrix Outer (Vector, Vector);
Vector Hadamard (Vector, Vector);
Vector* GramSchmidt (Vector*, int, int p = 2);
Vector Convolution (Vector, Vector);
Matrix Convolution (Matrix, Matrix);