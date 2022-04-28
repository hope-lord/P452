# include <iostream>
# include <cmath>
using namespace std;

const double e = M_E;

double f(double x) {
	return exp(-x*x);
}

// Montecarlo integration with uniform distribution
double Uniform_Sample_Integration(int N){
	double sum = 0;
	double x;
	for (int i = 0; i < N; i++) {
		x = (double)rand() / RAND_MAX;
		sum += f(x); // w(x) = 1
	}
	return sum / N;
}

// Montecarlo integration with Exponential distribution
double Expoential_Sample_Integration(int N){
	double sum = 0;
	double x;
	for (int i = 0; i < N; i++) {
		x = -log(1 - (e-1)/e*(double)rand() / RAND_MAX); // generating Exponential random number form 0 to 1 from uniform random numbers 
		sum += (e-1)*f(x)/exp(-x)/e; // w(x) = e*exp(-x)/(e-1)
	}
	return sum / N;
}


int main(){
	srand(1598731);
	// Store N vs integration for both methods in a file
	FILE * file = fopen("q1_Integration.txt", "w");
	for (int N = 100; N <= 100000000; N *= 10) {
		fprintf(file, "%d %lf %lf\n", N, Uniform_Sample_Integration(N), Expoential_Sample_Integration(N));
	}
	fclose(file);
}
