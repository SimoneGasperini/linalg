#include <cmath>
#include <iomanip>
#include "linalg.hpp"

Matrix::Matrix (int m, int n) {
    rows = m;
    cols = n;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0;
        }
    }
}

Matrix::Matrix (int m, int n, double** arr) {
    rows = m;
    cols = n;
    matrix = arr;
}

Matrix::Matrix (const Matrix& mat) {
    rows = mat.rows;
    cols = mat.cols;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = mat.matrix[i][j];
        }
    }
}

int Matrix::GetRows () {return rows;}

int Matrix::GetCols () {return cols;}

double Matrix::GetElement (int i, int j) {return matrix[i][j];}

void Matrix::SetElement (int i, int j, double val) {matrix[i][j] = val;}

bool Matrix::IsContained (int i, int* indici, int dim) {
    bool result = false;
    for (int k = 0; k < dim; k++) {
        if (i == indici[k]) {
            result = true;
        }
    }
    return result;
}

Matrix Matrix::GetMatrix (int k) {
    int riga = 0, colonna;
    Matrix mat(rows-1, cols-1);
    for (int i = 0; i < rows; i++) {
        colonna = 0;
        if (i != k) {
            for (int j = 1; j < cols; j++) {
                    mat.matrix[riga][colonna] = matrix[i][j];
                    colonna++;
            }
            riga++;
        }
    }
    return mat;
}

Matrix Matrix::GetMatrix (int* indici, int dim) {
    Matrix mat(rows-dim, cols-dim);
    int riga = 0, colonna;
    for (int i = 0; i < rows; i++) {
        colonna = 0;
        if (!IsContained(i, indici, dim)) {
            for (int j = 0; j < cols; j++) {
                if (!IsContained(j, indici, dim)) {
                    mat.matrix[riga][colonna] = matrix[i][j];
                    colonna++;
                }
            }
            riga++;
        }
    }
    return mat;
}

Matrix Matrix::GetCombs (int dim, int m, int n) {
    Matrix result(m, n);
    int* vett = new int[n];
    for (int i = 0; i < n; i++) {
        vett[i] = i + 1;
    }
    int i = 0;
    while (i < m) {
        int j = 0;
        while (j < n) {
            result.matrix[i][j] = vett[j]; 
            j++;      
        }
        int z = n - 1;
        while (z >= 0) {
            int max = dim - n + z + 1;  
            if (vett[z] < max) {     
                vett[z]++;
                int k = z;
                while (k < n - 1) {
                    vett[k + 1] = vett[k] + 1;
                    k++;
                }
                break;
            }
            z--;
        }
        i++;
    }  
    return result;
}

Matrix Matrix::Merge (const Matrix& right) {
    Matrix composta(rows, cols*2);
    for (int i = 0; i < composta.rows; i++) {
        for (int j = 0; j < composta.cols; j++) {
            if (j < cols) {
                composta.matrix[i][j] = matrix[i][j];
            } else {
                composta.matrix[i][j] = right.matrix[i][j-composta.rows];
            }
        }
    }
    return composta;
}

void Matrix::SwapRows () {
    double *temp = new double[rows];
    for (int i = 0; i < rows-1; i++) {
		if (matrix[i][i] == 0) {
			for (int p = i+1; p < rows; p++) {
				if (matrix[p][i] != 0) {
					for (int k = 0; k < cols; k++) {
						temp[k] = matrix[i][k];
						matrix[i][k] = matrix[p][k];
						matrix[p][k] = temp[k];
					}
				}
				break;
			}
		}
    }
}

void Matrix::Approx () {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (abs(matrix[i][j]) < APPROX) {
                matrix[i][j] = 0;
            }
        }
    }
}

double Matrix::Norm (char p) {
    double result;
    //maximum of the singular values
    if (p == '2') {
        int n = rows < cols ? rows : cols;
        double* sigmas = SingularValues();
        result = sigmas[0];
        for (int i = 1; i < n; i++) {
            if (sigmas[i] > result) {
                result = sigmas[i];
            }
        }
    }
    //frobenius norm
    if (p == 'f') {
        result = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result += matrix[i][j]*matrix[i][j];
            }
        }
    }
    //nuclear norm (works fine only for positive semi-definite matrices)
    if (p == 'n') {
        result = 0;
        int n = rows < cols ? rows : cols;
        double* sigmas = SingularValues();
        for (int i = 0; i < n; i++) {
            result += sigmas[i];
        }
    }
    return result;
}

Matrix Matrix::T () {
    Matrix trasp(cols, rows);
    for (int i = 0; i < trasp.rows; i++) {
        for (int j = 0; j < trasp.cols; j++) {
            trasp.matrix[i][j] = matrix[j][i];
        }
    }
    return trasp;
}

Matrix Matrix::I () {
    if (Determinant() == 0) {throw 4;}
    Matrix id = Eye(rows);
    Matrix composta = this->Merge(id);
    Matrix gauss = composta.Gauss();
    double a, g;
    for (int i = gauss.rows-1; i > 0; i--) {
		if (gauss.matrix[i][i] != 0) {
			for (int q = i-1; q >= 0; q--) {
				a = gauss.matrix[q][i];
				g = gauss.matrix[i][i];
				for (int j = 0; j < gauss.cols; j++) {
                    gauss.matrix[q][j] -= (gauss.matrix[i][j] / g) * a ;
				}
			}
		}
	}
    double q;
    for (int i = 0; i < gauss.rows; i++) {
        q = gauss.matrix[i][i];
		if (q != 1 && q != 0) {
			for (int j = 0; j < gauss.cols; j++) {
			    gauss.matrix[i][j] /= q;
			}
		}
	}
    Matrix inv(rows, cols);
    for (int i = 0; i < inv.rows; i++) {
        for (int j = 0; j < inv.cols; j++) {
            inv.matrix[i][j] = gauss.matrix[i][j+cols];
        }
    }
    return inv;
}

Matrix Matrix::Triu () {
    Matrix triu = *this;
    for (int i = 1; i < triu.rows; i++) {
        for (int j = 0; j < i; j++) {
            triu.matrix[i][j] = 0;
        }
    }
    return triu;
}

Matrix Matrix::Diag () {
    Matrix diag(rows, cols);
    for (int i = 0; i < diag.rows; i++) {
        diag.matrix[i][i] = matrix[i][i];
    }
    return diag;
}

Matrix Matrix::Gauss () {
    Matrix gauss = *this;
    gauss.SwapRows();
    double a, g;
    for (int i = 0; i < gauss.rows-1; i++) {
		for (int q = i+1; q < gauss.rows; q++) {
			a = gauss.matrix[q][i];
			g = gauss.matrix[i][i];
			for (int j = 0; j < gauss.cols; j++) {
				gauss.matrix[q][j] -= (gauss.matrix[i][j] / g) * a ;
			}
		}
	}
    return gauss;
}

int Matrix::Rank () {
    Matrix gauss = Gauss();
    int rank = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (abs(gauss.matrix[i][j]) > APPROX) {
                rank++;
                break;
            }
        }
    }
    return rank;
}

double Matrix::Trace () {
    if (rows != cols) {throw 2;}
    double result = 0;
    for (int i = 0; i < rows; i++) {
        result += matrix[i][i];
    }
    return result;
}

double Matrix::Determinant (char method) {
    if (rows != cols) {throw 2;}
    double result;
    if (method == 'g') {
        result = 1;
        Matrix triang = Gauss();
        for (int i = 0; i < triang.rows; i++) {
            result *= triang.matrix[i][i];
        }
    }
    if (method == 'l') {
        if (rows == 1) {
            result = matrix[0][0];
        } else {
            result = 0;
            for (int k = 0; k < rows; k++) {
                Matrix temp = GetMatrix(k);
                result += pow(-1,k)*matrix[k][0]*temp.Determinant();
            }
        }
    }
    return result;
}

double* Matrix::Eigenvalues () {
    //works fine only for symmetric matrices
    double* lambdas = new double[cols];
    Matrix mat = *this;
    for (int iter = 0; iter < 1000; iter++) {
        Matrix* QR = mat.QRdecomposition();
        mat = QR[1] * QR[0];
    }
    for (int i = 0; i < cols; i++) {
        lambdas[i] = mat.matrix[i][i];
    }
    return lambdas;
}

double* Matrix::SingularValues () {
    //works fine apart from the sign ambiguity
    int n = rows < cols ? rows : cols;
    Matrix S = n == rows ? (*this)*this->T() : this->T()*(*this);
    double* sigmas = S.Eigenvalues();
    for (int i = 0; i < n; i++) {
        sigmas[i] = pow(sigmas[i], 0.5);
    }
    return sigmas;
}

Polynomial Matrix::CharacteristicPol () {
    auto Fattoriale = [] (int n) {
        int res = 1;
        for (int i = 2; i <= n; i++) {
            res = res * i;
        }
        return res;
    };
    Polynomial car (rows);
    int iCoeff = car.GetDegree();
    car.SetCoefficient(iCoeff, pow(-1, iCoeff));
    car.SetCoefficient(iCoeff-1, pow(-1, iCoeff-1)*Trace());
    for (int k = iCoeff-2; k > 0; k--) {
        int segno = pow(-1 ,k);   
        int comb = Fattoriale(iCoeff)/(Fattoriale(iCoeff-k)*Fattoriale(k));    
        Matrix temp1 = GetCombs(rows, comb, k);
        double minore = 0;
        for (int t = 0; t < comb; t++) {
            int* indici = new int[k];
            for (int z = 0; z < k; z++) {
                indici[z] = (int)temp1.matrix[t][z]-1;
            }
            Matrix temp2 = GetMatrix(indici, k);
            minore += temp2.Determinant();
        }
        car.SetCoefficient(k, segno*minore);
    }
    car.SetCoefficient(0, Determinant());
    return car;
}

Matrix* Matrix::QRdecomposition () {
    if (rows != cols) {throw 2;}
    Matrix* QR = new Matrix[2];
    Vector* columns = new Vector[cols];
    for (int i = 0; i < cols; i++) {
        double* array = new double[rows];
        for (int j = 0; j < rows; j++) {
            array[j] = matrix[j][i];
        }
        columns[i] = Vector(rows, array);
    }
    Vector* ortho = GramSchmidt(columns, cols);
    double** q = new double*[rows];
    for (int i = 0; i < rows; i++) {
        q[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            q[i][j] = ortho[j].GetArray()[i];
        }
    }
    QR[0] = Matrix(rows, cols, q);
    double** r = new double*[rows];
    for (int i = 0; i < rows; i++) {
        r[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            if (i > j) {
                r[i][j] = 0;
            } else {
                r[i][j] = Dot(ortho[i], columns[j]);
            }
        }
    }
    QR[1] = Matrix(rows, cols, r);
    return QR;
}

ostream& operator << (ostream& os, Matrix& mat) {
    mat.Approx();
    for (int i = 0; i < mat.rows; i++) {
        os << "[";
        for (int j = 0; j < mat.cols; j++) {
            os << setprecision(3) << setw(8) << mat.matrix[i][j] << ' ';
        }    
        os << "]\n";
    }
}

istream& operator >> (istream& is, Matrix& mat) {
    if (cin) {
        clog << "\nDigita il numero di righe: ";
        is >> mat.rows;
    }
    if (cin) {
        clog << "\nDigita il numero di colonne: ";
        is >> mat.cols;
    }
    clog << "\n";
    mat.matrix = new double*[mat.rows];
    for (int i = 0; i < mat.rows; i++) {
        mat.matrix[i] = new double[mat.cols];
    }
    if (cin) {
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                clog << "Mat[" << i+1 << "][" << j+1 << "] = ";
                is >> mat.matrix[i][j];
            }
            clog << "\n";
        }
    }
}

Matrix& Matrix::operator = (const Matrix& mat) {
    rows = mat.rows;
    cols = mat.cols;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = mat.matrix[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator + (const Matrix& mat) {
    if (rows != mat.rows || cols != mat.cols) {throw 3;}
    Matrix somma(rows, cols);
    for (int i = 0; i < somma.rows; i++) {
        for (int j = 0; j < somma.cols; j++) {
            somma.matrix[i][j] = matrix[i][j] + mat.matrix[i][j];
        }
    }
    return somma;
}

Matrix Matrix::operator - (const Matrix& mat) {
    if (rows != mat.rows || cols != mat.cols) {throw 3;}
    Matrix diff(rows, cols);
    for (int i = 0; i < diff.rows; i++) {
        for (int j = 0; j < diff.cols; j++) {
            diff.matrix[i][j] = matrix[i][j] - mat.matrix[i][j];
        }
    }
    return diff;
}

Matrix Matrix::operator * (const Matrix& mat) {
    if (cols != mat.rows) {throw 3;}
    Matrix prod(rows, mat.cols);
    for (int i = 0; i < rows; i++) {
        for (int k = 0; k < cols; k++) {
            for (int j = 0; j < mat.cols; j++) {
            prod.matrix[i][j] += matrix[i][k] * mat.matrix[k][j];
            }
        }
    }
    return prod;
}

Matrix Matrix::operator * (double k) {
    Matrix res(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res.matrix[i][j] = matrix[i][j] * k;
        }
    }
    return res;
}

Matrix Matrix::operator / (double k) {
    Matrix res(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res.matrix[i][j] = matrix[i][j] / k;
        }
    }
    return res;
}