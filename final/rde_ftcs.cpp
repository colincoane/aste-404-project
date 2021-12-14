/*
ASTE-404 Final Project
1D Simulation of a Rotating Detonation Rocket Engine (RDE)
Nondimensional simulation of detonation wave energy
that captures basic physics underlying an RDE while
avoiding the complexities and nonlinearities of reaction 
kinetics and 3-dimensional detonation phenomena
*/

// Header Files
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <iomanip>
#include "Vec.h"

// Define world structure for storing cell data
struct World { // world structure
  double3 x0;	// origin
  double3 dh;   // cell spacing
  int3 nn;  	// number of nodes
  int U(int i, int j, int k=0) {return k*nn[1]*nn[0]+j*nn[0]+i;}
};


// function prototypes
void saveVTI(int ts, const World &world, const std::vector<double> &den);


int main() {

    //********** USER DEFINED PARAMETERS **********//

    // Domain Parameters
    double L = 2.16*M_PI;   // Domain length
    int ni = 161;           // number of simulation nodes (x-dir)
	int nj = 26;            // plot height (y-dir) - for plotting only, no bearing on simulation

    double dt = 1e-4;           // time step
	int nt = 350001;			// number of timesteps in simulation

    // Output parameters

    // output snapshot of system to vti for animation
    int steps_per_snapshot = 1000;  // save snapshot in .vti every x timesteps

    // output entire time history of simulation to csv
    std::string out_file_name = "time_history.csv"; // time history output file name
    int steps_per_csv = 100;        // save time history in csv every x timesteps

    // set physical constants for simulation (nondimensional)
    double q = 1.0;             // propellant heat release
    double alpha = 0.3;         // activation energy
    double u_c = 1.1;           // ignition energy / temperature
    double s = 4.5;             // injector area
    double k0 = 1;              // reaction rate constant
    double epsilon = 0.8;       // energy loss coefficient
    double r = 5;               // injector activation     
    double u_p = 0.8;           // injector pressure

    // diffusion to force continuous solutions
    double nu1 = 0.0075;        // combustion diffusitivity
    double nu2 = 0.0075;        // injection diffusivity

    //********** END USER DEFINED PARAMETERS **********//


    //************** INITIALIZATION **************//

    // set world variables
    World world;
    world.x0 = {0.0,0.0,0.0};		// origin for plotting
    world.dh = {0.02,0.02,0.02};	// cell spacing (arbitrary)

    // reassign x-direction node spacing based on domain length
    world.dh[0] = L*1/(ni-1);

    world.nn = {ni,nj,1};           // number of cells
	int nu = ni*nj;                 // total nodes per timestep
	
	// set temporal & spatial constants
	double dx = world.dh[0];    // x-cell spacing
	double dy = world.dh[1];    // y-cell spacing


    //************** CREATE SOLUTION VECTORS **************//

    // set solution vectors
    std::vector<double> u(nu);           // internal energy
    std::vector<double> lambda(nu);      // propellant mass fraction

    // initial conditions, loop over solution vectors
    for (int j = 0; j < nj; j++) {
        for (int i = 0; i < ni; i++) {
            int n = j*ni + i;   // unknown indexing

            // initial energy spike to facilitate detonation
            u[n] = 1.5*pow((1/cosh(1*(i*dx-1))),20);    // sech()^20 function

            // constant initial combustion progress
            lambda[n] = 0.5;        // initially half burnt, half unburnt
        }
    }


    //************** INTEGRATE SYSTEM **************//

    std::ofstream out_file(out_file_name); // open file for output

    for (int t = 0; t < nt; t++) {  // FTCS time marching

        // create copy of current step properties
        std::vector<double> u_old = u;          // energy at current step
        std::vector<double> lam_old = lambda;   // mass fraction at current step

        // Apply FTCS to internal nodes and solve properties at next step
        for (int j = 0; j < nj; j++) {
            for (int i = 0; i < ni; i++) {
                int n = j*ni + i;   // unknown indexing

                // Apply periodic boundary conditions to simulate annulus
                if (i == 0) {           // left-most node circles back to right-most node
                    
                    // solve for energy
                    u[n] = u_old[n] + (nu1*dt/(dx*dx))*(u_old[n+1] - 2*u_old[n] + u_old[n+ni-1])
                        - (dt/(2.0*dx))*u_old[n]*(u_old[n+1] - u_old[n+ni-1])
                        + dt*k0*q*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*epsilon*u_old[n]*u_old[n];

                    // solve for mass fraction
                    lambda[n] = lam_old[n] + (nu2*dt/(dx*dx))*(lam_old[n+1] - 2*lam_old[n] + lam_old[n+ni-1])
                        + dt*k0*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*s*u_p*lam_old[n]/(1 + exp(r*(u_old[n] - u_p)));
                    
                } else if (i == ni-1) { // right-most node circles back to left-most node
                    
                    // solve for energy
                    u[n] = u_old[n] + (nu1*dt/(dx*dx))*(u_old[n-(ni-1)] - 2*u_old[n] + u_old[n-1])
                        - (dt/(2.0*dx))*u_old[n]*(u_old[n-(ni-1)] - u_old[n-1])
                        + dt*k0*q*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*epsilon*u_old[n]*u_old[n];

                    // solve for mass fraction
                    lambda[n] = lam_old[n] + (nu2*dt/(dx*dx))*(lam_old[n-(ni-1)] - 2*lam_old[n] + lam_old[n-1])
                        + dt*k0*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*s*u_p*lam_old[n]/(1 + exp(r*(u_old[n] - u_p)));
                    
                } else {                // standard internal node
                    
                    // solve for energy
                    u[n] = u_old[n] + (nu1*dt/(dx*dx))*(u_old[n+1] - 2*u_old[n] + u_old[n-1])
                        - (dt/(2.0*dx))*u_old[n]*(u_old[n+1] - u_old[n-1])
                        + dt*k0*q*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*epsilon*u_old[n]*u_old[n];

                    // solve for mass fraction
                    lambda[n] = lam_old[n] + (nu2*dt/(dx*dx))*(lam_old[n+1] - 2*lam_old[n] + lam_old[n-1])
                        + dt*k0*(1 - lam_old[n])*exp((u_old[n] - u_c) / alpha)
                        - dt*s*u_p*lam_old[n]/(1 + exp(r*(u_old[n] - u_p)));
                    
                }
            }
        }
        // end FTCS loop, output results if desired

        // Note: results output can be modified to output lambda if desired

        // plot snapshot every x time steps
		if (t%steps_per_snapshot==0) { std::cout<<t<<std::endl; saveVTI(t,world,u);}

        // output time history results to csv file every x time steps
        if (t%steps_per_csv == 0) {
            // since snapshot is repeated nj times, only take
            // from first ni values in the energy vector u
            for (int idx = 0; idx < ni; idx++) {
                if (idx < ni - 1) {
                    out_file<<u[idx]<<",";  // populate row in csv
                } else {
                    out_file<<u[idx]<<"\n"; // end row, next snapshot will fill next row
                }
            }
        }

        // end loop for single timestep
    }
    // end time integration
    
	return 0;	// normal exit
}   // end of main function

// writes snapshot of energy to a MANUALLY CREATED "results" folder
void saveVTI(int ts, const World &world, const std::vector<double> &den)  {
	// output data
	std::stringstream ss;
	ss<<"results/field_"<<std::setfill('0')<<std::setw(6)<<ts<<".vti";

	std::ofstream out(ss.str());            //open output file
	out<<"<?xml version=\"1.0\"?>\n";
	out<<"<VTKFile type=\"ImageData\">\n";
	out<<"<ImageData WholeExtent=\"0 "<<world.nn[0]-1<<" 0 "<<world.nn[1]-1<<" 0 "<<world.nn[2]-1<<"\" "
	   <<"Origin=\""<<world.x0<<"\" Spacing=\""<<world.dh<<"\">\n";
	out<<"<Piece Extent=\"0 "<<world.nn[0]-1<<" 0 "<<world.nn[1]-1<<" 0 "<<world.nn[2]-1<<"\">\n";

	out<<"<PointData>\n";
	
	out<<"<DataArray Name=\"u (-)\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">\n";
	for (const double &val:den) out<<val<<" ";
	out<<"</DataArray>\n";

	out<<"</PointData>\n";

	out<<"</Piece>\n";
	out<<"</ImageData>\n";
	out<<"</VTKFile>\n";
}