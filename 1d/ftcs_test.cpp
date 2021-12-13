/*
ASTE-404 Project Test Bed
*/

#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<math.h>
#include <iomanip>
#include "Matrix.h"
#include "Vec.h"


using namespace std;

struct World { // world structure
  double3 x0;	// origin
  double3 dh;   // cell spacing
  int3 nn;  	// number of nodes
  int U(int i, int j, int k=0) {return k*nn[1]*nn[0]+j*nn[0]+i;}
};


// function prototypes
vector<double> solveGS(Matrix &A, vector<double> &g, bool verbose=false);
void saveVTI(int ts, const World &world, const vector<double> &den);


int main() {

    //************** INITIALIZATION **************//

    // set world variables
    World world;
    world.x0 = {0.0,0.0,0.0};		// origin for plotting
    world.dh = {0.02,0.02,0.02};	// cell spacing
    world.nn = {161,26,1};           // number of cells
	
	int ni = world.nn[0];
	int nj = world.nn[1];
	int nu = ni*nj;

    double L = 2.16*M_PI;   // Domain length
	
	// set temporal / spatial constants
    world.dh[0] = L*1/(ni-1);
	double dx = world.dh[0];    // x spacing
	double dy = world.dh[1];    // y spacing
    double dt = 1e-4;           // time step
	int nt = 350001;			// number of timesteps

    // set physical constants
    //double L = 2*M_PI;
    double q = 1.0;
    double alpha = 0.3;
    double u_c = 1.1;
    double s = 4.5;
    double k0 = 1;
    double eps = 0.8;
    double r = 5;
    double nu1 = 0.0075;
    double nu2 = 0.0075;
    double u_p = 0.8;

    //************** END INITIALIZATION **************//



    //************** CREATE VECTORS AND MATRICES **************//

    // set solution vectors
    vector<double> u(nu);           // internal energy
    vector<double> lambda(nu);      // reactant mass fraction

    // initialize pulse characteristics
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;
            u[n] = 1.5*pow((1/cosh(1*(i*dx-1))),20);
            lambda[n] = 0.5;
        }
    }

    // Coefficients to save space
    // TODO

    //************** SOLVE FTCS **************//

    // initial pressure or energy (possibly needed, arbitrary rn)
    //P = EtoP(E,rho,u,v,Y,gamma,q);
    //E = PtoE(P,rho,u,v,Y,gamma,q);

    for (int t = 0; t < nt; t++) {  // FTCS time marching

        // create copy of previous step properties
        vector<double> u_old = u;
        vector<double> lam_old = lambda;

        // Apply FTCS to internal nodes
        for (int j = 0; j < nj; j++) {
            for (int i = 0; i < ni; i++) {
                int n = j*ni + i;   // indexing
                // Apply periodic boundary conditions
                if (i == 0) {
                    // solve for energy
                    u[n] = u_old[n] + (nu1*dt/(dx*dx))*(u_old[n+1] - 2*u_old[n] + u_old[n+ni-1])
                        - (dt/(2.0*dx))*u_old[n]*(u_old[n+1] - u_old[n+ni-1])
                        + dt*k0*q*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*eps*u_old[n]*u_old[n];
                    // solve for mass fraction
                    lambda[n] = lam_old[n] + (nu2*dt/(dx*dx))*(lam_old[n+1] - 2*lam_old[n] + lam_old[n+ni-1])
                        + dt*k0*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*s*u_p*lam_old[n]/(1 + exp(r*(u_old[n] - u_p)));
                } else if (i == ni-1) {
                    // solve for energy
                    u[n] = u_old[n] + (nu1*dt/(dx*dx))*(u_old[n-(ni-1)] - 2*u_old[n] + u_old[n-1])
                        - (dt/(2.0*dx))*u_old[n]*(u_old[n-(ni-1)] - u_old[n-1])
                        + dt*k0*q*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*eps*u_old[n]*u_old[n];
                    // solve for mass fraction
                    lambda[n] = lam_old[n] + (nu2*dt/(dx*dx))*(lam_old[n-(ni-1)] - 2*lam_old[n] + lam_old[n-1])
                        + dt*k0*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*s*u_p*lam_old[n]/(1 + exp(r*(u_old[n] - u_p)));
                } else {
                    // solve for energy
                    u[n] = u_old[n] + (nu1*dt/(dx*dx))*(u_old[n+1] - 2*u_old[n] + u_old[n-1])
                        - (dt/(2.0*dx))*u_old[n]*(u_old[n+1] - u_old[n-1])
                        + dt*k0*q*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*eps*u_old[n]*u_old[n];
                    // solve for mass fraction
                    lambda[n] = lam_old[n] + (nu2*dt/(dx*dx))*(lam_old[n+1] - 2*lam_old[n] + lam_old[n-1])
                        + dt*k0*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*s*u_p*lam_old[n]/(1 + exp(r*(u_old[n] - u_p)));
                }
            }
        }

        // plot results every 100 time steps
		if (t%1000==0) { cout<<t<<endl; saveVTI(t/10,world,u);}
    }
    

	return 0;	// normal exit
}   // main end


// Gauss Seidel solver, receives references to a base Matrix type
vector<double> solveGS(Matrix &A, vector<double> &g, bool verbose) {

    int nr = A.getNr();     // number of rows

    vector<double> x(nr);   // solution vector

    /* solve matrix system */
    for (int it = 0; it < 10000; it++) { // solver iteration
        for (int r = 0; r < nr; r++) { // loop over rows
            
            double sum = A.dotRow(r,x) - A(r,r)*x[r];
            double x_star = (g[r] - sum) / A(r,r); // new estimate for x[r]

            x[r] += 1.4*(x_star-x[r]);    // SOR
        }

        // convergence check, only every 50 iterations
        if(it%50 == 0) {

            vector<double> R = A*x - g; // residue vector

            // compute average error
            double L2 = sqrt(mag2(R)/nr);

            if (verbose) {  // output iteration data if desired
                cout<<"solver iteration "<<it<<", L2 norm: "<<L2<<endl;
            }
            if (L2 < 1e-4) break; // break out of solver loop
        }
    }

    // return solution vector, hint to compiler to try moving instead of copying
    return std::move(x);
}

// writes out density field to a "results" folder (needs to be manually created!)
void saveVTI(int ts, const World &world, const vector<double> &den)  {
	// output data
	stringstream ss;
	ss<<"results/field_"<<setfill('0')<<setw(6)<<ts<<".vti";

	ofstream out(ss.str());            //open output file
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData WholeExtent=\"0 "<<world.nn[0]-1<<" 0 "<<world.nn[1]-1<<" 0 "<<world.nn[2]-1<<"\" "
	   <<"Origin=\""<<world.x0<<"\" Spacing=\""<<world.dh<<"\">\n";
	out<<"<Piece Extent=\"0 "<<world.nn[0]-1<<" 0 "<<world.nn[1]-1<<" 0 "<<world.nn[2]-1<<"\">\n";

	out<<"<PointData>\n";
	
	out<<"<DataArray Name=\"den (#/m^3)\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">\n";
	for (const double &val:den) out<<val<<" ";
	out<<"</DataArray>\n";

	out<<"</PointData>\n";

	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
}