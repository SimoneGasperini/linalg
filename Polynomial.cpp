#include <cmath>
#include "linalg.hpp"

Polynomial::Polynomial (int k) {
    degree = k;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = 0;
    }
}

Polynomial::Polynomial (int k, double *coeff) {
    degree = k;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = coeff[i];
    }
}

Polynomial::Polynomial (const Polynomial& pol) {
    degree = pol.degree;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = pol.coefficients[i];
    }
}

int Polynomial::GetDegree () {return degree;}

double* Polynomial::GetCoefficients () {return coefficients;}

void Polynomial::SetCoefficient (int indice, double valore) {
    coefficients[indice] = valore;
}

void Polynomial::Approx () {
    for (int i = 0; i < degree+1; i++) {
        if (abs(coefficients[i]) < APPROX) {
            coefficients[i] = 0;
        }
    }
}

double Polynomial::Evaluate (double x) {
    double valore = 0;  
    for (int i = degree; i >= 0; i--) {
        valore += pow(x, i)*coefficients[i];
    }
    return valore;
}

double Polynomial::EvaluateDerivative (double x) {
    double valore = 0;
    for (int i = degree; i > 0; i--) {
        valore += i*pow(x, i-1)*coefficients[i];
    }
    return valore;
}

double Polynomial::ComputeZero () {
    double rad = static_cast <double> (rand());
    int i = 0;
    int MAXITER = 1000;
    while (abs(Evaluate(rad)) > APPROX) {
        if (i > MAXITER) {
            throw 1;
        }
        rad -= Evaluate(rad)/EvaluateDerivative(rad);
        i++;
    }
    return rad;
}

double* Polynomial::ComputeZeros () {
    double* radici = new double[degree];
    Polynomial div, temp;
    double rad;
    if (coefficients[0] == 0) {
        rad = 0;
    } else {      
        rad = ComputeZero();
    }
    radici[0] = rad;
    int curIndex = 1;
    if (degree > 1) {
        Polynomial div(1);
        div.SetCoefficient(1, 1);
        div.SetCoefficient(0, -rad);
        temp = *this / div;
        double* tempRadici = temp.ComputeZeros();
        for (int i = 0; i < temp.GetDegree(); i++) {
            radici[curIndex++] = tempRadici[i];
        }
    }
    return radici;
}

ostream& operator << (ostream& os, Polynomial& pol) {
    for (int k = pol.degree; k > 0; k--) {   
        double corrente = abs(pol.coefficients[k]);
        if (pol.coefficients[k] != 0) {
            os << (pol.coefficients[k] > 0 ? '+' : '-');
            if (corrente != 1) {
                os << corrente;
            }
            if (k != 1) {
                os << "x^" << k;
            } else {
                os << "x";
            }
            os << ' ';
        }
    }
    if (pol.coefficients[0] != 0) {
        os << (pol.coefficients[0] > 0 ? '+' : '-') << abs(pol.coefficients[0]);
    }
    return os;
}

istream& operator >> (istream& is, Polynomial& pol) {
    int n;
    if (cin) {
        clog << "Digita il degree del polinomio: ";
        is >> n;
    }
    if (n > pol.degree) {
        pol.coefficients = new double[n+1];
        pol.degree = n;
    }
    if (cin) {
        clog << "Digita i " << n+1 << " coefficients del polinomio: \n";
        for (int k = n; k >= 0; k--) { 
            is >> pol.coefficients[k];
        }
    }
    return is;
}

Polynomial& Polynomial::operator = (const Polynomial& pol) {
    degree = pol.degree;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = pol.coefficients[i];
    }
    return *this;
}

Polynomial Polynomial::operator + (const Polynomial& pol) {
    Polynomial max, min;
    if (degree > pol.degree) {
        max = *this;
        min = pol;
    } else {
        max = pol;
        min = *this;  
    }
    Polynomial somma(max.degree);
    for (int i = 0; i <= somma.degree; i++) { 
        somma.coefficients[i] = max.coefficients[i] + (i <= min.degree ? min.coefficients[i] : 0);
    }       
    return somma;
}

Polynomial Polynomial::operator - (const Polynomial& pol) {
    Polynomial diff(pol.degree <= degree ? degree : pol.degree);
    for (int i = 0; i <= diff.degree; i++) { 
        diff.coefficients[i] = coefficients[i] - pol.coefficients[i];
    }       
    return diff;
}

Polynomial Polynomial::operator * (const Polynomial& pol) {
    Vector v1(degree+1, coefficients);
    Vector v2(pol.degree+1, pol.coefficients);
    Vector conv = Convolution(v1, v2);
    Polynomial prod(degree+pol.degree, conv.GetArray());
    return prod;
}

Polynomial Polynomial::operator / (const Polynomial& pol) {
    Polynomial ret(degree-pol.degree);
    ret.coefficients[ret.degree] = coefficients[degree] / pol.coefficients[pol.degree];
    for (int i = 1; i <= ret.degree; i++) {
        double* u = new double[degree-i+2];
        for (int k = 0; k <= degree-i+1; k++) {
            u[k] = coefficients[k];
        }
        Polynomial temp1(degree-i);
        Polynomial temp2(degree-i+1, u);
        temp1.coefficients[degree-i] = ret.coefficients[ret.degree-i+1];
        temp1 = temp1 * pol,
        temp2 = temp2 - temp1;
        ret.coefficients[ret.degree-i] = temp2.coefficients[degree-i] / pol.coefficients[pol.degree];
    }
    return ret;
}