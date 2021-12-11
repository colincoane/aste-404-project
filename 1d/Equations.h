#ifndef _EQN_H
#define _EQN_H

#include <vector>
#include <math.h>
#include "Matrix.h"

// Include various math equations to help with computations
// TODO fix up functions and make cleaner

// functions for solving equations of state and key equations over fields
std::vector<double> iglPres(std::vector<double> &rho, double R, std::vector<double> &T) {
    // Solve Ideal Gas Equation for Pressure
    int sr = rho.size();
    int sT = T.size();

    std::vector<double> P(sr); // solution vector, pressure

    if (sr == sT) { // solve equation
        for (int i = 0; i < sr; i++) {
            P[i] = rho[i]*R*T[i]; // P = rho*R*T
        }
    } else {
        for (int i = 0; i < sr; i++) {
            P[i] = 0; // return zeroes if vectors are incorrect sizes
        }
    }

    // return solution vector, hint to compiler to try moving instead of copying
    return std::move(P);
}

std::vector<double> iglTemp(std::vector<double> &P, std::vector<double> &rho, double R) {
    // Solve Ideal Gas Equation for Temperature
    int sr = rho.size();
    int sP = P.size();

    std::vector<double> T(sr); // solution vector, temperature

    if (sr == sP) { // solve equation
        for (int i = 0; i < sr; i++) {
            T[i] = P[i]/(rho[i]*R); // T = P/rho*r
        }
    } else {
        for (int i = 0; i < sr; i++) {
            T[i] = 0; // return zeroes if vectors are incorrect sizes
        }
    }

    // return solution vector, hint to compiler to try moving instead of copying
    return std::move(T);
}

std::vector<double> EtoP(std::vector<double> &E, std::vector<double> &rho, std::vector<double> &u,
                            std::vector<double> &v, std::vector<double> &Y, double gamma,
                            double q)
{
    // Calculate pressure P from energy E
    int sr = rho.size();

    std::vector<double> P(sr);

    for (int i = 0; i < sr; i++) {
        P[i] = (E[i] - rho[i]*q*Y[i] - 0.5*(u[i]*u[i] + v[i]*v[i]))*(gamma - 1);
    }

    return std::move(P);
}

std::vector<double> PtoE(std::vector<double> &P, std::vector<double> &rho, std::vector<double> &u,
                            std::vector<double> &v, std::vector<double> &Y, double gamma,
                            double q)
{
    // Calculate energy E from pressure P
    int sr = rho.size();

    std::vector<double> E(sr);

    for (int i = 0; i < sr; i++) {
        E[i] = P[i]/(gamma - 1) + 0.5*(u[i]*u[i] + v[i]*v[i]) + rho[i]*q*Y[i];
    }

    return std::move(E);
}

std::vector<double> wCalc(double K, double Ti, std::vector<double> &rho, std::vector<double> &Y,
                            std::vector<double> &T)
{
    // Calculate source term w
    int sr = rho.size();

    std::vector<double> w(sr);

    for (int i = 0; i < sr; i++) {
        w[i] = -K*rho[i]*Y[i]*exp(-Ti / T[i]);
    }

    return std::move(w);
}


#endif