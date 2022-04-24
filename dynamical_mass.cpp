#include<iostream>
#include<vector>
#include<cmath>
using namespace std;


const double pi = M_PI;
const double D = 45;
const double w2 = 0.16;
const double m0 = 0;
const double eps = 1e-8;
const double eps2 = sqrt(-w2*log(1e-18));
const int N = 501;
const int Nz = 31;
const double MAX_ty = log10(1000);
const double MIN_ty = -4;
const double GuessA = 1.4;
const double GuessB = 0.5;


typedef vector<double> VECTOR;

// Calculates the roots of legendre polynomial and their respective weights
void gauleg(const double x1, const double x2, VECTOR &x,VECTOR &w){
    const int N=x.size();  // number of points
    const double EPS=1e-14; // tolerance
    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for(i=0;i<N;i++){
        
        z=cos(pi*(i+0.75)/(N+0.5));
        while(abs(z-z1)>EPS){  // Newton's method of root finding
            p1 =1;
            p2 = 0;
            for(j=0;j<N;j++){  // Legendre polynomial at z
                p3 = p2;
                p2 = p1;
                p1 = ((2*j+1)*z*p2-j*p3)/(j+1);
            }
            pp = N*(z*p1-p2)/(z*z-1); // derivative of Legendre polynomial at z
            z1 = z;
            z = z1-p1/pp;
        }
        x[i] = xm-xl*z;
        x[N-1-i] = xm+xl*z;
        w[i] = 2*xl/((1-z*z)*pp*pp);  // weight Calculation
        w[N-1-i] = w[i]; // since the weights are the symmetic about midpoint
    }
}

// The integral of z in the integration of A(x) 
double fza(double z,double tx,double ty){
	return 2*sqrt(1-z*z)*(-2*exp10(ty)/3+(1+exp10(ty-tx))*exp10((tx+ty)/2)*z-4*exp10(ty)*z*z/3)*
		exp((-exp10(tx)-exp10(ty)+2*z*exp10((tx+ty)/2))/w2)/pi/w2;
}

// The integral of z in the integration of B(x)
double fzb(double z,double tx,double ty){
	return 2*sqrt(1-z*z)*(exp10(tx)+exp10(ty)-exp10((ty+tx)/2)*2*z)*exp((-exp10(tx)-exp10(ty)+2*z*exp10((tx+ty)/2))/w2)/pi/w2;
}

// yA/(y^2A+B)
double integral_a(double a,double b,double tx,double ty){
	return D*exp10(2*ty)*log(10)*a/(exp10(ty)*a*a+b*b);// * int_fza(tx,ty);
}

// yB/(y^2A+B)
double integral_b(double a,double b,double tx,double ty){
	return D*exp10(2*ty)*log(10)*b/(exp10(ty)*a*a+b*b); //* int_fzb(tx,ty);
}

// integration of fza function
double integrate_fza(double tx,double ty,VECTOR & z,VECTOR & wz){
	double result = 0;
	int M = z.size();
	for(int i = 0;i<M;i++){
		result += wz[i]*fza(z[i],tx,ty);
	}
	return result;
}

//integration of fzb function
double integrate_fzb(double tx,double ty,VECTOR & z,VECTOR & wz){
	double result = 0;
	int M = z.size();
	for(int i = 0;i<M;i++){
		result += wz[i]*fzb(z[i],tx,ty);
	}
	return result;
}


int main(void){
	VECTOR wz(Nz),z(Nz);
	gauleg(-1,1,z,wz); // Caclulate the roots of legendre polynomial and their respective weights in -1 to 1
	VECTOR ty(N);VECTOR w(N);
	gauleg(MIN_ty,MAX_ty,ty,w); // Calculate the roots of legendre polynomial and their respective weights in MIN_ty to MAX_ty
	VECTOR &tx = ty;
	
	VECTOR A(N,GuessA); //Guess the intial function A(x)
	VECTOR B(N,GuessB); // Guess the initial function B(x)
	VECTOR A1(N);
	VECTOR B1(N);
	bool converged_A = false;
	bool converged_B = false;
	int iter = 1;
	double sumA,sumB;
	while(!converged_A or !converged_B ){
		converged_A = converged_B = true;
		cout << "Iteration number : "<<iter <<endl;
		iter++;
		for(int i=0;i<N;i++){
			sumA = sumB = 0;
			int j = i;
			// IGNORE THE POINTS IF THEY ARE TOO FAR CURRENT POINT
			while(j>-1 and (std::abs(exp10(ty[j]/2)-exp10(tx[i]/2)<eps2))){
				sumA += w[j]*integral_a(A[j],B[j],tx[i],ty[j])*integrate_fza(tx[i],ty[j],z,wz);
				sumB += w[j]*integral_b(A[j],B[j],tx[i],ty[j])*integrate_fzb(tx[i],ty[j],z,wz);
				j--;
			}
			j = i+1;
			while(j<N and (std::abs(exp10(ty[j]/2)-exp10(tx[i]/2)<eps2))){
				sumA += w[j]*integral_a(A[j],B[j],tx[i],ty[j])*integrate_fza(tx[i],ty[j],z,wz);
				sumB += w[j]*integral_b(A[j],B[j],tx[i],ty[j])*integrate_fzb(tx[i],ty[j],z,wz);
				j++;
			}
			A1[i] = sumA+1;
			B1[i] = sumB+m0;
			converged_A *= (std::abs(A1[i]-A[i])<eps);
			converged_B *= (std::abs(B1[i]-B[i])<eps);
		}	
		swap(A,A1);
		swap(B,B1);
	}
	 // Print the result to a file
    FILE* data_file = fopen("out_data.txt","w");
    fprintf(data_file,"# The values of D, w2, m0 are respectively:\n%g %g %g\n",D,w2,m0);
    fprintf(data_file,"# x A(x) B(x)\n");
    for(int i=0;i<N;i++)
        fprintf(data_file,"%lg %lg %lg\n",tx[i],A[i],B[i]);
    fclose(data_file);
	cout << "M(0) = " << B[0]/A[0] << endl;
	return 0;
}
