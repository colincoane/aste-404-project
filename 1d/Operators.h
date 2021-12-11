#ifndef _OPER_H
#define _OPER_H

#include <vector>
#include <math.h>
#include "Matrix.h"
#include "Equations.h"

// Set up operators for solving each equation using the Crank Nicolson scheme

// Set vector initial conditions
void setInitialConds(double dt, int ni, int nj, std::vector<double> &rho, std::vector<double> &P,
                        std::vector<double> &T, std::vector<double> &u, std::vector<double> &v,
                        std::vector<double> &E, std::vector<double> &Y)
{
    // Set initial conditions for all vectors

    double coeff = dt/2; // coefficient for crank nicolson

    /* set matrix values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0 && i < 15) { // Dirichlet on inlet
                
                // TODO set values

                rho[n] = 1000.0;
                P[n] = 100000.0;
                T[n] = 600.0;
                u[n] = 0.8;
                v[n] = 0.2;
                E[n] = 10.0;
                Y[n] = 1.0;
            } else if (j == nj-1) { // Initial conditions on outlet
                
                // TODO set values

                rho[n] = 10.0;
                P[n] = 10000.0;
                T[n] = 300.0;
                u[n] = 0.8;
                v[n] = 0.2;
                E[n] = 1.0;
                Y[n] = 1.0;
            } else {    // regular internal node
                rho[n] = 10.0;
                P[n] = 10000.0;
                T[n] = 300.0;
                u[n] = 0.8;
                v[n] = 0.2;
                E[n] = 1.0;
                Y[n] = 1.0;
            }
            // Set high density shockwave
            if (i >= 5 && i <= 7) {
                rho[n] = 1000.0;
            }
        }
    }
}

// Set matrix values here

// Set Rho matrix values
void setRho(double dt, double dx, double dy, int ni, int nj, FiveBandMatrix &Arho,
                        FiveBandMatrix &Brho, std::vector<double> &u, std::vector<double> &v)
{
    // Set initial conditions for all vectors

    double coeff = dt/2; // coefficient for crank nicolson
    double D = 0.001;         // false diffsion coefficient
    double difc2 = (D*dt/2.0);  // set fake diffusion for convergence purposes

    /* set matrix values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on inlet
                
                // TODO set values
                Arho(n,n) = 1.0;
                Brho(n,n) = 1.0;
            
            } else if (j == nj-1) { // Backward difference on outlet

                // set A coefficients, A = I - dt/2 * f
                Arho(n,n-ni) = -difc2*(2.0/(dy*dy)) - coeff*v[n-ni]/(dy);
                Arho(n,n-1) = -difc2*(2.0/(dx*dx)) - coeff*u[n-1]/(dx);
                Arho(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) + coeff*v[n]/(dy) + coeff*u[n]/(dx);
                
                // set B coefficients, B = I + dt/2 * f
                Brho(n,n-ni) = difc2*(2.0/(dy*dy)) + coeff*v[n-ni]/(dy);
                Brho(n,n-1) = difc2*(2.0/(dx*dx)) + coeff*u[n-1]/(dx);
                Brho(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) - coeff*v[n]/(dy) - coeff*u[n]/(dx);
            } else if (i == 0) { // Periodic on left side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Arho(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*v[n-ni]/(2*dy);
                Arho(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*u[n-1+ni]/(2*dx);
                Arho(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Arho(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*u[n+1]/(2*dx);
                Arho(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Brho(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*v[n-ni]/(2*dy);
                Brho(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*u[n-1+ni]/(2*dx);
                Brho(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Brho(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*u[n+1]/(2*dx);
                Brho(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*v[n+ni]/(2*dy);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Arho(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*v[n-ni]/(2*dy);
                Arho(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*u[n-1]/(2*dx);
                Arho(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Arho(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*u[n+1-ni]/(2*dx);
                Arho(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Brho(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*v[n-ni]/(2*dy);
                Brho(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*u[n-1]/(2*dx);
                Brho(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Brho(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*u[n+1-ni]/(2*dx);
                Brho(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*v[n+ni]/(2*dy);

            } else {    // regular internal node
                
                // TODO check values

                // set A coefficients, A = I - dt/2 * f
                Arho(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*v[n-ni]/(2*dy);
                Arho(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*u[n-1]/(2*dx);
                Arho(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Arho(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*u[n+1]/(2*dx);
                Arho(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Brho(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*v[n-ni]/(2*dy);
                Brho(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*u[n-1]/(2*dx);
                Brho(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Brho(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*u[n+1]/(2*dx);
                Brho(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*v[n+ni]/(2*dy);
            }
        }
    }
}


// Set u matrix values
void setUvel(double dt, double dx, double dy, int ni, int nj, FiveBandMatrix &Au,
                        FiveBandMatrix &Bu, std::vector<double> &rho, std::vector<double> &P,
                        std::vector<double> &u_old, std::vector<double> &v)
{
    // Set initial conditions for all vectors

    double coeff = dt/2; // coefficient for crank nicolson
    double D = 0.001;         // false diffsion coefficient
    double difc2 = (D*dt/2.0);  // set fake diffusion for convergence purposes

    /* set matrix values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on inlet
                
                // TODO set values
                Au(n,n) = 1.0;
                Bu(n,n) = 1.0;
            
            } else if (j == nj-1) { // Backward difference on outlet
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(2.0/(dy*dy)) -coeff*rho[n-ni]*v[n-ni]/(dy);
                Au(n,n-1) = -difc2*(2.0/(dx*dx)) -coeff*rho[n-1]*u_old[n-1]/(dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) + coeff*rho[n]*v[n]/(dy) + coeff*rho[n]*u_old[n]/(dx);
                
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(2.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(dy);
                Bu(n,n-1) = difc2*(2.0/(dx*dx)) + coeff*rho[n-1]*u_old[n-1]/(dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) - coeff*rho[n]*v[n]/(dy) - coeff*rho[n]*u_old[n]/(dx);;

            } else if (i == 0) { // Periodic on left side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(1.0/(dy*dy)) -coeff*rho[n-ni]*v[n-ni]/(2*dy);
                Au(n,n-1) = -difc2*(1.0/(dx*dx)) -coeff*rho[n-1+ni]*u_old[n-1+ni]/(2*dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Au(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u_old[n+1]/(2*dx);
                Au(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(2*dy);
                Bu(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1+ni]*u_old[n-1+ni]/(2*dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Bu(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u_old[n+1]/(2*dx);
                Bu(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v[n+ni]/(2*dy);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(1.0/(dy*dy)) -coeff*rho[n-ni]*v[n-ni]/(2*dy);
                Au(n,n-1) = -difc2*(1.0/(dx*dx)) -coeff*rho[n-1]*u_old[n-1]/(2*dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Au(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1-ni]*u_old[n+1-ni]/(2*dx);
                Au(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(2*dy);
                Bu(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u_old[n-1]/(2*dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Bu(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1-ni]*u_old[n+1-ni]/(2*dx);
                Bu(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v[n+ni]/(2*dy);

            } else {    // regular internal node
                
                // TODO check values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(1.0/(dy*dy)) -coeff*rho[n-ni]*v[n-ni]/(2*dy);
                Au(n,n-1) = -difc2*(1.0/(dx*dx)) -coeff*rho[n-1]*u_old[n-1]/(2*dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Au(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u_old[n+1]/(2*dx);
                Au(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(2*dy);
                Bu(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u_old[n-1]/(2*dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Bu(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u_old[n+1]/(2*dx);
                Bu(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v[n+ni]/(2*dy);
            }
        }
    }
}


// Set v matrix values
void setVvel(double dt, double dx, double dy, int ni, int nj, FiveBandMatrix &Au,
                        FiveBandMatrix &Bu, std::vector<double> &rho, std::vector<double> &P,
                        std::vector<double> &u, std::vector<double> &v_old)
{
    // Set initial conditions for all vectors

    double coeff = dt/2; // coefficient for crank nicolson
    double D = 0.001;         // false diffsion coefficient
    double difc2 = (D*dt/2.0);  // set fake diffusion for convergence purposes

    /* set matrix values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on inlet
                
                // TODO set values
                Au(n,n) = 1.0;
                Bu(n,n) = 1.0;
            
            } else if (j == nj-1) { // Backward difference on outlet
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(2.0/(dy*dy)) -coeff*rho[n-ni]*v_old[n-ni]/(dy);
                Au(n,n-1) = -difc2*(2.0/(dx*dx)) -coeff*rho[n-1]*u[n-1]/(dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) + coeff*rho[n]*v_old[n]/(dy) + coeff*rho[n]*u[n]/(dx);
                
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(2.0/(dy*dy)) + coeff*rho[n-ni]*v_old[n-ni]/(dy);
                Bu(n,n-1) = difc2*(2.0/(dx*dx)) + coeff*rho[n-1]*u[n-1]/(dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) - coeff*rho[n]*v_old[n]/(dy) - coeff*rho[n]*u[n]/(dx);;

            } else if (i == 0) { // Periodic on left side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(1.0/(dy*dy)) -coeff*rho[n-ni]*v_old[n-ni]/(2*dy);
                Au(n,n-1) = -difc2*(1.0/(dx*dx)) -coeff*rho[n-1+ni]*u[n-1+ni]/(2*dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Au(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u[n+1]/(2*dx);
                Au(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v_old[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v_old[n-ni]/(2*dy);
                Bu(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1+ni]*u[n-1+ni]/(2*dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Bu(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u[n+1]/(2*dx);
                Bu(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v_old[n+ni]/(2*dy);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(1.0/(dy*dy)) -coeff*rho[n-ni]*v_old[n-ni]/(2*dy);
                Au(n,n-1) = -difc2*(1.0/(dx*dx)) -coeff*rho[n-1]*u[n-1]/(2*dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Au(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1-ni]*u[n+1-ni]/(2*dx);
                Au(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v_old[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v_old[n-ni]/(2*dy);
                Bu(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u[n-1]/(2*dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Bu(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1-ni]*u[n+1-ni]/(2*dx);
                Bu(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v_old[n+ni]/(2*dy);

            } else {    // regular internal node
                
                // TODO check values

                // set A coefficients, A = I - dt/2 * f
                Au(n,n-ni) = -difc2*(1.0/(dy*dy)) -coeff*rho[n-ni]*v_old[n-ni]/(2*dy);
                Au(n,n-1) = -difc2*(1.0/(dx*dx)) -coeff*rho[n-1]*u[n-1]/(2*dx);
                Au(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Au(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u[n+1]/(2*dx);
                Au(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v_old[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                Bu(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v_old[n-ni]/(2*dy);
                Bu(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u[n-1]/(2*dx);
                Bu(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                Bu(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u[n+1]/(2*dx);
                Bu(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v_old[n+ni]/(2*dy);
            }
        }
    }
}


// Set E matrix values
void setE(double dt, double dx, double dy, int ni, int nj, FiveBandMatrix &AE,
                        FiveBandMatrix &BE, std::vector<double> &P,
                        std::vector<double> &u, std::vector<double> &v)
{
    // Set initial conditions for all vectors

    double coeff = dt/2; // coefficient for crank nicolson
    double D = 0.001;         // false diffsion coefficient
    double difc2 = (D*dt/2.0);  // set fake diffusion for convergence purposes

    /* set matrix values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on inlet
                
                // TODO set values
                AE(n,n) = 1.0;
                BE(n,n) = 1.0;
            
            } else if (j == nj-1) { // Backward difference on outlet

                // set A coefficients, A = I - dt/2 * f
                AE(n,n-ni) = -difc2*(2.0/(dy*dy)) - coeff*v[n-ni]/(dy);
                AE(n,n-1) = -difc2*(2.0/(dx*dx)) - coeff*u[n-1]/(dx);
                AE(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) + coeff*v[n]/(dy) + coeff*u[n]/(dx);
                
                // set B coefficients, B = I + dt/2 * f
                BE(n,n-ni) = difc2*(2.0/(dy*dy)) + coeff*v[n-ni]/(dy);
                BE(n,n-1) = difc2*(2.0/(dx*dx)) + coeff*u[n-1]/(dx);
                BE(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) - coeff*v[n]/(dy) - coeff*u[n]/(dx);
            } else if (i == 0) { // Periodic on left side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                AE(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*v[n-ni]/(2*dy);
                AE(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*u[n-1+ni]/(2*dx);
                AE(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                AE(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*u[n+1]/(2*dx);
                AE(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                BE(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*v[n-ni]/(2*dy);
                BE(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*u[n-1+ni]/(2*dx);
                BE(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                BE(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*u[n+1]/(2*dx);
                BE(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*v[n+ni]/(2*dy);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                AE(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*v[n-ni]/(2*dy);
                AE(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*u[n-1]/(2*dx);
                AE(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                AE(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*u[n+1-ni]/(2*dx);
                AE(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                BE(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*v[n-ni]/(2*dy);
                BE(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*u[n-1]/(2*dx);
                BE(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                BE(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*u[n+1-ni]/(2*dx);
                BE(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*v[n+ni]/(2*dy);

            } else {    // regular internal node
                
                // TODO check values

                // set A coefficients, A = I - dt/2 * f
                AE(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*v[n-ni]/(2*dy);
                AE(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*u[n-1]/(2*dx);
                AE(n,n) = 1.0 + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                AE(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*u[n+1]/(2*dx);
                AE(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                BE(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*v[n-ni]/(2*dy);
                BE(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*u[n-1]/(2*dx);
                BE(n,n) = 1.0 - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                BE(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*u[n+1]/(2*dx);
                BE(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*v[n+ni]/(2*dy);
            }
        }
    }
}



// Set Y matrix values
void setY(double dt, double dx, double dy, int ni, int nj, FiveBandMatrix &AY,
                        FiveBandMatrix &BY, std::vector<double> &rho,
                        std::vector<double> &u, std::vector<double> &v)
{
    // Set initial conditions for all vectors

    double coeff = dt/2; // coefficient for crank nicolson
    double D = 0.001;         // false diffsion coefficient
    double difc2 = (D*dt/2.0);  // set fake diffusion for convergence purposes

    /* set matrix values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on inlet
                
                // TODO set values
                AY(n,n) = 1.0;
                BY(n,n) = 1.0;
            
            } else if (j == nj-1) { // Backward difference on outlet

                // set A coefficients, A = I - dt/2 * f
                AY(n,n-ni) = -difc2*(2.0/(dy*dy)) - coeff*rho[n-ni]*v[n-ni]/(dy);
                AY(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*rho[n-1]*u[n-1]/(dx);
                AY(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) + coeff*rho[n]*v[n]/(dy) + coeff*rho[n]*u[n]/(dx);
                
                // set B coefficients, B = I + dt/2 * f
                BY(n,n-ni) = difc2*(2.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(dy);
                BY(n,n-1) = difc2*(2.0/(dx*dx)) + coeff*rho[n-1]*u[n-1]/(dx);
                BY(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy)) - coeff*rho[n]*v[n]/(dy) - coeff*rho[n]*u[n]/(dx);
                
            } else if (i == 0) { // Periodic on left side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                AY(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*rho[n-ni]*v[n-ni]/(2*dy);
                AY(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*rho[n-1]*u[n-1+ni]/(2*dx);
                AY(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                AY(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u[n+1]/(2*dx);
                AY(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                BY(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(2*dy);
                BY(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u[n-1+ni]/(2*dx);
                BY(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                BY(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u[n+1]/(2*dx);
                BY(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v[n+ni]/(2*dy);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values

                // set A coefficients, A = I - dt/2 * f
                AY(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*rho[n-ni]*v[n-ni]/(2*dy);
                AY(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*rho[n-1]*u[n-1]/(2*dx);
                AY(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                AY(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u[n+1-ni]/(2*dx);
                AY(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                BY(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(2*dy);
                BY(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u[n-1]/(2*dx);
                BY(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                BY(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u[n+1-ni]/(2*dx);
                BY(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v[n+ni]/(2*dy);

            } else {    // regular internal node
                
                // TODO check values

                // set A coefficients, A = I - dt/2 * f
                AY(n,n-ni) = -difc2*(1.0/(dy*dy)) - coeff*rho[n-ni]*v[n-ni]/(2*dy);
                AY(n,n-1) = -difc2*(1.0/(dx*dx)) - coeff*rho[n-1]*u[n-1]/(2*dx);
                AY(n,n) = rho[n] + difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                AY(n,n+1) = -difc2*(1.0/(dx*dx)) + coeff*rho[n+1]*u[n+1]/(2*dx);
                AY(n,n+ni) = -difc2*(1.0/(dy*dy)) + coeff*rho[n+ni]*v[n+ni]/(2*dy);
                // set B coefficients, B = I + dt/2 * f
                BY(n,n-ni) = difc2*(1.0/(dy*dy)) + coeff*rho[n-ni]*v[n-ni]/(2*dy);
                BY(n,n-1) = difc2*(1.0/(dx*dx)) + coeff*rho[n-1]*u[n-1]/(2*dx);
                BY(n,n) = rho[n] - difc2*(2.0/(dx*dx) + 2.0/(dy*dy));
                BY(n,n+1) = difc2*(1.0/(dx*dx)) - coeff*rho[n+1]*u[n+1]/(2*dx);
                BY(n,n+ni) = difc2*(1.0/(dy*dy)) - coeff*rho[n+ni]*v[n+ni]/(2*dy);
            }
        }
    }
}

// derivative of P with respect to y
std::vector<double> diffPy(double dx, double dy, int ni, int nj,std::vector<double> &P)
{
    // create vector
    int nu = ni*nj;
    std::vector<double> dPy(nu);

    /* set vector values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on bottom
                
                // TODO set values
                dPy[n] = 0;

            } else if (j == nj-1) { // Dirichlet on top
                
                // TODO set values
                dPy[n] = 0;

            } else {    // regular internal node
                
                // TODO check values
                dPy[n] = (P[n+ni]-P[n-ni])/(2*dy);
            }
        }
    }
    return std::move(dPy);
}


// derivative of P with respect to x
std::vector<double> diffPx(double dx, double dy, int ni, int nj,std::vector<double> &P)
{
    // create vector
    int nu = ni*nj;
    std::vector<double> dPx(nu);

    /* set vector values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (i == 0) { // Periodic on left side
                
                // TODO set values
                dPx[n] = (P[n+1]-P[n+ni-1])/(2*dx);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values
                dPx[n] = (P[n-ni+1]-P[n-1])/(2*dx);

            } else {    // regular internal node
                
                // TODO check values
                dPx[n] = (P[n+1]-P[n-1])/(2*dx);
            }
        }
    }

    return std::move(dPx);
}


// derivative of P*u with respect to x
std::vector<double> diffPUx(double dx, double dy, int ni, int nj,std::vector<double> &P, std::vector<double> &u)
{
    // create vector
    int nu = ni*nj;
    std::vector<double> dPu(nu);

    /* set vector values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (i == 0) { // Periodic on left side
                
                // TODO set values
                dPu[n] = (P[n+1]*u[n+1]-P[n+ni-1]*u[n+ni-1])/(2*dx);

            } else if (i == ni-1) { // Periodic on right side
                
                // TODO set values
                dPu[n] = (P[n-ni+1]*u[n-ni+1]-P[n-1]*u[n-1])/(2*dx);

            } else {    // regular internal node
                
                // TODO check values
                dPu[n] = (P[n+1]*u[n+1]-P[n-1]*u[n-1])/(2*dx);
            }
        }
    }

    return std::move(dPu);
}


// derivative of P*v with respect to y
std::vector<double> diffPVy(double dx, double dy, int ni, int nj,std::vector<double> &P,std::vector<double> &v)
{
    // create vector
    int nu = ni*nj;
    std::vector<double> dPv(nu);

    /* set vector values */
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            if (j == 0) { // Dirichlet on bottom
                
                // TODO set values
                dPv[n] = 0;

            } else if (j == nj-1) { // Dirichlet on top
                
                // TODO set values
                dPv[n] = 0;

            } else {    // regular internal node
                
                // TODO check values
                dPv[n] = (P[n+ni]*v[n+ni]-P[n-ni]*v[n-ni])/(2*dy);
            }
        }
    }
    return std::move(dPv);
}



#endif