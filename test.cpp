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
#include "Equations.h"
#include "Operators.h"

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
    world.nn = {51,26,1};           // number of cells
	
	int ni = world.nn[0];
	int nj = world.nn[1];
	int nu = ni*nj;
	
	// set temporal / spatial constants
	double dx = world.dh[0];    // x spacing
	double dy = world.dh[1];    // y spacing
    double dt = 1e-3;           // time step
	int nt = 5001;			// number of timesteps

    // set physical constants
    double K = 1;
    double q = 50;
    double Ti = 300;
    double gamma = 1.4;
    double R = 287;

    //************** END INITIALIZATION **************//



    //************** CREATE VECTORS AND MATRICES **************//

    // set solution vectors
    vector<double> rho(nu);     // density
    vector<double> P(nu);       // pressure
    vector<double> T(nu);       // temperature
    vector<double> u(nu);       // x velocity
    vector<double> v(nu);       // y velocity
    vector<double> E(nu);       // internal energy
    vector<double> Y(nu);       // reactant mass fraction

    // set matrices for crank nicolson
    // create a nu*nu banded matrix with ni-offset
    FiveBandMatrix Arho(nu,ni);      // solve for rho
    FiveBandMatrix Brho(nu,ni);      // solve for rho

    FiveBandMatrix Au(nu,ni);        // solve for u
    FiveBandMatrix Bu(nu,ni);        // solve for u

    FiveBandMatrix Av(nu,ni);        // solve for v
    FiveBandMatrix Bv(nu,ni);        // solve for v

    FiveBandMatrix AE(nu,ni);        // solve for E
    FiveBandMatrix BE(nu,ni);        // solve for E

    FiveBandMatrix AY(nu,ni);        // solve for Y
    FiveBandMatrix BY(nu,ni);        // solve for Y

    //************** SET MATRIX VALUES **************//

    // assign initial and boundary conditions to vectors
    setInitialConds(dt,ni,nj,rho,P,T,u,v,E,Y);

    // set matrix values
    setRho(dt,dx,dy,ni,nj,Arho,Brho,u,v);
    setUvel(dt,dx,dy,ni,nj,Au,Bu,rho,P,u,v);
    setVvel(dt,dx,dy,ni,nj,Av,Bv,rho,P,u,v);
    setE(dt,dx,dy,ni,nj,AE,BE,P,u,v);
    setY(dt,dx,dy,ni,nj,AY,BY,rho,u,v);

    //************** SOLVE CRANK NICOLSON **************//

    // initial pressure or energy (possibly needed, arbitrary rn)
    //P = EtoP(E,rho,u,v,Y,gamma,q);
    //E = PtoE(P,rho,u,v,Y,gamma,q);

    for (int i = 0; i < nt; i++) {  // Crank Nicolson time marching

        // Evaluate rho
        vector<double> b1 = Brho*rho;
        // solve Gauss Seidel
        rho = solveGS(Arho,b1);
        // reset matrix values
        setRho(dt,dx,dy,ni,nj,Arho,Brho,u,v);

        // Calculate P, T
        //P = EtoP(E,rho,u,v,Y,gamma,q);
        //T = iglTemp(P,rho,R);

        /********** NOTE *********/
        // Code currently set up to only solve rho
        // idea is that you would solve each diffeq and then update the matrices
        // A and B needed for crank nicolson ( An^k+1 = Bn^k + R)
        // Right now rho, u, v work/converge when solved by themselves, E and Y do not
        // matrix elements should be correct  for rho, u, v, unsure about E and Y
        // Possible idea is to set a constant flow field (u,v) so no need to solve them
        // Equations for rho, E, Y have similar forms, so could kind of generalize the 
        // solution for rho to a general "variable" such as energy
        // Then set initial conditions for rho or something as a moving shockwave
        // Set everything but rho as constant (now rho could be energy not density)
        // Add the w source/reaction term to the equation to simulate stuff happening
        // Add also an arbitrary "loss" term (i.e. to remove energy/simulate all other variables)
        // Essentially, idea is to solve the equation with Y in it while keeping most other
        // variables constant
        // we have the E = p/(gamma - 1) + ... equation and ideal gas law to relate certain properties
        // but if we only have to solve one differential equation that may be ideal
        // and whe can kind of throw everything else under the rug

        // Evaluate u
        vector<double> dPx = diffPx(dx,dy,ni,nj,P);
        vector<double> b2 = Bu*u - dPx;
        // solve Gauss Seidel
        //u = solveGS(Au,b2);
        // reset matrix values
        //setUvel(dt,dx,dy,ni,nj,Au,Bu,rho,P,u,v);

        // Evaluate v
        vector<double> dPy = diffPy(dx,dy,ni,nj,P);
        vector<double> b3 = Bv*v - dPy;
        // solve Gauss Seidel
        //v = solveGS(Av,b3);
        // reset matrix values
        //setVvel(dt,dx,dy,ni,nj,Av,Bv,rho,P,u,v);

        // Evaluate E
        vector<double> dPu = diffPUx(dx,dy,ni,nj,P,u);
        vector<double> dPv = diffPVy(dx,dy,ni,nj,P,v);
        vector<double> b4 = BE*E - dPu - dPv;
        // solve Gauss Seidel
        //E = solveGS(AE,b4);
        // reset matrix values
        //setE(dt,dx,dy,ni,nj,AE,BE,P,u,v);

        // Calculate P, T
        P = EtoP(E,rho,u,v,Y,gamma,q);
        T = iglTemp(P,rho,R);

        // Evaluate Y
        vector<double> w = wCalc(K,Ti,rho,Y,T);
        vector<double> b5 = BY*Y + w;
        // solve Gauss Seidel
        //Y = solveGS(AY,b5);
        // reset matrix values
        //setY(dt,dx,dy,ni,nj,AY,BY,rho,u,v);

        // plot results every 100 time steps
		if (i%100==0) { cout<<i<<endl; saveVTI(i,world,rho);}
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