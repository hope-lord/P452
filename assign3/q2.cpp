/*
  Here I have assumed that h = m_e = 1. So energy of ground state = 1/8 and energy of excited state = 1/2.
  So for ground state the Schroedinger equation is: d^2y/dx^2 = -y*pi^2
  For 1st excited state the Schroedinger equation is: d^2y/dx^2 = -y*4*pi^2
 */

# include "package_cpp/UTILITY.h"
# include "package_cpp/ode.cpp"
# include "package_cpp/io.cpp"

/* const double pi = M_PI; */

VECTOR F(double x,VECTOR &y,void * energy){
	double E = *(double *)energy;
	VECTOR derivs(2);
	derivs[0] = y[1];
	derivs[1] = -y[0]*E*8*pi*pi;
	return derivs;
}

void EnergyState(double energy,double guess1,double guess2,string f1,string f2,string f3){
    double x0 = 0; // left position of wall
	double x1 = 1; // right position of wall
	double y0 = 0,y1 = 0; // Wave function vanishes at x = 0 and x = 1
	int N = 100; // Number of points in between 0,1
    // With first guess
	VECTOR init = {y0,guess1};
	MATRIX out = rk4(&F,init,x0,x1,&energy,N);
	FILE *file1 = fopen(f1.c_str(),"w");
	print_to_file(file1,out);
    // With second guess
	init[1] = guess2;
	out = rk4(&F,init,x0,x1,&energy,N);
	FILE *file2 = fopen(f2.c_str(),"w");
	print_to_file(file2,out);
    // With Shooting Method
	out = Dirichlet_Boundary(&F,y0,y1,x0,x1,guess1,guess2,&energy,1e-6,100);
	FILE *file3 = fopen(f3.c_str(),"w");
	print_to_file(file3,out);
}

int main(){
	double E0 = 1.0/8.0,E1 = 0.5;
	EnergyState(E0,3,3.5,"q2_ground_guess1.txt","q2_ground_guess2.txt","q2_ground_sol.txt");

	EnergyState(E1,6,6.5,"q2_1st_guess1.txt","q2_1st_guess2.txt","q2_1st_sol.txt");
	return 0;
}
