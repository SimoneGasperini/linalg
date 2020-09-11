#include <cmath>
#include <iomanip>
#include "linalg.hpp"

Vector::Vector (int s) {
    size = s;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = 0;
    }
}

Vector::Vector (int s, double* arr) {
    size = s;
    array = arr;
}

Vector::Vector (const Vector& vec) {
    size = vec.size;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = vec.array[i];
    }
}

int Vector::GetSize () {return size;}

double* Vector::GetArray () {return array;}

double Vector::GetElement (int i) {return array[i];}

void Vector::SetElement (int i, double val) {array[i] = val;}

double Vector::Norm (int p) {
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += pow(abs(array[i]), p);
    }
    result = pow(result, 1/(double)p);
    return result;
}

Vector Vector::Normalized (int p) {
    Vector vec = *this / Norm(p);
    return vec;
}

Vector Vector::ProjectedOnto (Vector vec) {
    double k = Dot(*this,vec) / Dot(vec,vec);
    Vector proj = vec*k;
    return proj;
}

ostream& operator << (ostream& os, Vector& vec) {
    os << "[";
    for (int i = 0; i < vec.size; i++) {
        os << setprecision(3) << vec.array[i];
        if (i != vec.size-1) {os << "  ";}
    }
    os << "]";
}

Vector& Vector::operator = (const Vector& vec) {
    size = vec.size;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = vec.array[i];
    }
    return *this;
}

Vector Vector::operator + (const Vector& vec) {
    if (size != vec.size) {throw 5;}
    Vector sum(size);
    for (int i = 0; i < size; i++) {
        sum.array[i] = array[i] + vec.array[i];
    }
    return sum;
}

Vector Vector::operator - (const Vector& vec) {
    if (size != vec.size) {throw 5;}
    Vector diff(size);
    for (int i = 0; i < size; i++) {
        diff.array[i] = array[i] - vec.array[i];
    }
    return diff;
}

Vector Vector::operator * (double k) {
    Vector res(size);
    for (int i = 0; i < size; i++) {
        res.array[i] = array[i] * k;
    }
    return res;
}

Vector Vector::operator / (double k) {
    Vector res(size);
    for (int i = 0; i < size; i++) {
        res.array[i] = array[i] / k;
    }
    return res;
}