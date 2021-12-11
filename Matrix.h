#ifndef _MAT_H
#define _MAT_H

#include <vector>
#include <stdexcept>

// Holds Matrix class and some std::vector operations/overloads
// NOTE: move std::vector operations to another header possibly

class Matrix {
public:

    Matrix(int nr) : nr{nr} {}  // set local nr with provided value of nr
    virtual ~Matrix() {}    // virtual empty destructor

    int getNr() {return nr;}    // gets value of protected variable

    // pure virtual function, to be implemented by specific matrix types
    // the "const" after function name tells the compiler that the function will
    // not be modifying any member (class Matrix) data, thus making it possible
    // to be used with a "const Matrix &A" parameter in operator*
    virtual double dotRow(int r, const std::vector<double>&x) const = 0;

    // operator called for A*x expression
    friend std::vector<double> operator*(const Matrix &A, const std::vector<double>&x) {

        // make a new empty vector for the result
        std::vector<double> y(A.nr);
        for (int r = 0; r < A.nr; r++) {
            y[r] = A.dotRow(r,x);
        }
        return y;   // this is where a "move" can come in handy
    }

    // overload () operator
    virtual double& operator()(int r, int c) = 0;   // for setting and write access
    virtual double operator()(int r, int c) const = 0;  // for read-only access

protected:
    int nr = 0; // number of rows, in protected block to limit access
};

// dense (square) matrix
class DenseMatrix : public Matrix {
public:
    // constructor chaining by calling Matrix::Matrix(int nr)
    DenseMatrix(int nr) : Matrix(nr) {
        // allocate memory
        a = new double*[nr];        // rows
        for (int r = 0; r < nr; r++)
            a[r] = new double[nr];  // columns
        
        // set all to zero
        for (int r = 0; r < nr; r++)
            for (int c = 0; c < nr; c++)
                a[r][c] = 0.0;
    }

    // destructor, free memory
    ~DenseMatrix() {
        for (int r = 0; r < nr; r++) delete[] a[r];
        delete[] a;
        a = nullptr;    // clear stored address
    }

    // overload () operator
    double& operator()(int r, int c) {return a[r][c];}      // reference
    double operator()(int r, int c) const {return a[r][c];} // value


    double dotRow(int r, const std::vector<double>&x) const {
        double sum = 0;
        for (int c = 0; c < nr; c++) sum += a[r][c]*x[c];
        return sum;
    }

    double **a = nullptr;   //coefficients
};

// identity matrix
class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int nr) : Matrix(nr) {}

    // overload () operator
    double& operator()(int r, int c) {
        throw std::runtime_error("Unsupported operation!");
    }
    double operator()(int r, int c) const {
        if (r == c) return 1; else return 0;
    }

    double dotRow(int r, const std::vector<double>&x) const {
        return x[r];
    }
};

// sparse 5-banded matrix
class FiveBandMatrix : public Matrix {
public:
    FiveBandMatrix(int nr, int d1) : Matrix(nr), d1{d1} {
        a = newAndClear(nr);
        b = newAndClear(nr);
        c = newAndClear(nr);
        d = newAndClear(nr);
        e = newAndClear(nr);
    }

    ~FiveBandMatrix() {
        delete[] a; a = nullptr;
        delete[] b; b = nullptr;
        delete[] c; c = nullptr;
        delete[] d; d = nullptr;
        delete[] e; e = nullptr;
    }

    /*
    // operator called for c*A expression
    FiveBandMatrix operator*(const FiveBandMatrix &A, const double c) {

        // make a new empty vector for the result
        FiveBandMatrix B(A.getNr(),A.getd1());
        for (int r = 0; r < A.getNr(); r++) {
            a[r] *= c;
        }
        return B;   // this is where a "move" can come in handy
    }
    */

    int getd1() {return d1;}    // gets value of protected variable

    // overloading () operator
    double& operator()(int r, int c) {
        if (c-r == -d1) return a[r];
        else if (c-r == -1) return b[r];
        else if (c-r == 0) return this->c[r]; // to distinguish argument from member
        else if (c-r == 1) return d[r];
        else if (c-r == d1) return e[r];
        // modifications for periodic boundary conditions
        else throw std::runtime_error("Unsupported operation!");
    }

    double operator()(int r, int c) const {
        if (c-r == -d1) return a[r];
        else if (c-r == -1) return b[r];
        else if (c-r == 0) return this->c[r]; // to distinguish argument from member
        else if (c-r == 1) return d[r];
        else if (c-r == d1) return e[r];
        else throw std::runtime_error("Unsupported operation!");
    }

    double dotRow(int r, const std::vector<double>&x) const {
        double sum = 0;
        
        if (r-d1 >= 0) sum += a[r]*x[r-d1];
        if (r-1 >= 0) sum += b[r]*x[r-1];
        sum += c[r]*x[r];
        if (r+1 < nr) sum += d[r]*x[r+1];
        if (r+d1 < nr) sum += e[r]*x[r+d1];

        // modifications for periodic boundary conditions
        // TODO CHECK THIS IS CORRECT
        if (r-d1 < 0) sum += a[r]*x[nr+(r-d1)];
        if (r-1 < 0) sum += b[r]*x[nr+(r-1)];
        if (r+1 >= nr) sum += d[r]*x[r+1-nr];
        if (r+d1 >= nr) sum += e[r]*x[r+d1-nr];
        return sum;
    }

protected:
    // helper function to allocate and clear an array
    static double *newAndClear(int nr) {
        double *p = new double[nr];
        for (int i = 0; i < nr; i++) p[i] = 0;
        return p;
    }


    double *a = nullptr;    //coefficients
    double *b = nullptr;    //coefficients
    double *c = nullptr;    //coefficients
    double *d = nullptr;    //coefficients
    double *e = nullptr;    //coefficients

    int d1;     // offset for a and e
                // assuming b and d offset by 1
};

// vector addition, returns c = a + b
std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b) {
    int n = a.size();
    std::vector<double> c(n);    // same size as a
    for (int i = 0; i < n; i++) c[i] = a[i] + b[i];
    return c;
}

// vector subtraction, returns c = a - b
std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b) {
    int n = a.size();
    std::vector<double> c(n);    // same size as a
    for (int i = 0; i < n; i++) c[i] = a[i] - b[i];
    return c;
}

// returns dot product with itself, magnitude squared
double mag2(std::vector<double> &a) {
    double sum = 0;
    for (double d:a) sum += d*d;
    return sum;
}




#endif