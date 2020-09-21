#include <cmath>
#include "linalg.hpp"

Matrix Eye (int n) {
    Matrix id(n,n);
    for (int i = 0; i < n; i++) {
        id.SetElement(i,i,1);
    }
    return id;
}

double Dot (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    if (s1 != s2) {throw 5;}
    double result = 0;
    for (int i = 0; i < s1; i++) {
        result += v1.GetElement(i) * v2.GetElement(i);
    }
    return result;
}

Vector Dot (Matrix mat, Vector vec) {
    int rows = mat.GetRows();
    int cols = mat.GetCols();
    int size = vec.GetSize();
    if (cols != size) {throw 5;}
    double* arr = new double[rows];
    for (int i = 0; i < rows; i++) {
        arr[i] = 0;
        for (int j = 0; j < cols; j++) {
            arr[i] += mat.GetElement(i,j) * vec.GetElement(j);
        }
    }
    Vector result(rows, arr);
    return result;
}

Vector Dot (Vector vec, Matrix mat) {
    int size = vec.GetSize();
    int rows = mat.GetRows();
    int cols = mat.GetCols();
    if (size != rows) {throw 5;}
    double* arr = new double[cols];
    for (int i = 0; i < cols; i++) {
        arr[i] = 0;
        for (int j = 0; j < rows; j++) {
            arr[i] += vec.GetElement(j) * mat.GetElement(i,j);
        }
    }
    Vector result(cols, arr);
    return result;
}

Matrix Outer (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    double** mat = new double*[s1];
    for (int i = 0; i < s1; i++) {
        mat[i] = new double[s2];
        for (int j = 0; j < s2; j++) {
            mat[i][j] = v1.GetElement(i) * v2.GetElement(j);
        }
    }
    Matrix result(s1, s2, mat);
    return result;
}

Vector Hadamard (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    if (s1 != s2) {throw 5;}
    double* arr = new double[s1];
    for (int i = 0; i < s1; i++) {
        arr[i] = v1.GetElement(i) * v2.GetElement(i);
    }
    Vector result(s1, arr);
    return result;
}

Vector* GramSchmidt (Vector* vecs, int num, int p) {
    Vector* orthogonal = new Vector[num];
    for (int i = 0; i < num; i++) {
        orthogonal[i] = vecs[i];
    }
    for (int i = 1; i < num; i++) {
        for (int k = 0; k < i; k++) {
            Vector proj = vecs[i].ProjectedOnto(orthogonal[k]);
            orthogonal[i] = orthogonal[i] - proj;
        }
    }
    Vector* orthonormal = new Vector[num];
    for (int i = 0; i < num; i++) {
        orthonormal[i] = orthogonal[i].Normalized(p);
    }
    return orthonormal;
}

Vector Convolution (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    int size = s1+s2-1;
    double* arr = new double[size];
    for (int i = 0; i < s1; i++) {
        for(int j = 0; j < s2; j++) {
            arr[i+j] += v1.GetElement(i) * v2.GetElement(j);
        }
    }
    Vector result(size, arr);
    return result;
}

Matrix Convolution (Matrix m1, Matrix m2) {
    int r1 = m1.GetRows(), c1 = m1.GetCols();
    int r2 = m2.GetRows(), c2 = m2.GetCols();
    int rows = r1+r2-1, cols = c1+c2-1;
    double** arr = new double*[rows];
    for (int i = 0; i < rows; i++) {
        arr[i] = new double[cols];
    }
    for (int i1 = 0; i1 < r1; i1++) {
        for (int j1 = 0; j1 < c1; j1++) {
            for (int i2 = 0; i2 < r2; i2++) {
                for (int j2 = 0; j2 < c2; j2++) { 
                    arr[i1+i2][j1+j2] += m1.GetElement(i1,j1) * m2.GetElement(i2,j2);
                }
            }
        }
    }
    Matrix result(rows, cols, arr);
    return result;
}