# include <iostream>
# include <cmath>
using namespace std;

const double e = M_E;

double f(double x) {
	return exp(-x*x);
}

// Montecarlo integration with uniform distribution
void Uniform_Sample_Integration(int N, double &result,double &error  ){
	result = 0;
	double square = 0;
	double x;
	for (int i = 0; i < N; i++) {
		x = (double)rand() / RAND_MAX;
		result += f(x); // w(x) = 1
		square += f(x)*f(x);
	}
	result *= (1.0 / N);
	square *= (1.0 / N);
	error = sqrt(square - result*result)/sqrt(N);
}

// Montecarlo integration with Exponential distribution
void Expoential_Sample_Integration(int N,double &result,double &error){
	result = 0;
	double square = 0;
	double x;
	for (int i = 0; i < N; i++) {
		x = -log(1 - (e-1)/e*(double)rand() / RAND_MAX); // generating Exponential random number form 0 to 1 from uniform random numbers 
		result += (e-1)*f(x)/exp(-x)/e; // w(x) = e*exp(-x)/(e-1)
		square += (e - 1)*(e - 1)*f(x)*f(x)/exp(-2*x)/e/e;
	}
	result *= (1.0 / N);
	square *= (1.0 / N);
	error = sqrt(square - result*result)/sqrt(N);
}


int main(){
	srand(1598731);
	// Store N vs integration for both methods in a file
	FILE * file1 = fopen("q1_uniform.txt", "w");
	FILE * file2 = fopen("q1_expon.txt", "w");
	double result,error;

	for (int N = 100; N <= 100000000; N *= 10) {
		cout << "At N = "<<N<<":" << endl;
		Uniform_Sample_Integration(N, result, error);
		fprintf(file1, "%d %lf %lf\n", N, result,error);
		cout << "Uniform: " << result << " +/- " << error << endl;
		Expoential_Sample_Integration(N, result, error);
		cout << "Exponential: " << result << " +/- " << error << endl;
		fprintf(file2, "%d %lf %lf\n", N, result,error);
	}
	fclose(file1);
	fclose(file2);
}
/*
 At N = 100:
Uniform: 0.74866 +/- 0.0197903
Exponential: 0.744452 +/- 0.00554889
At N = 1000:
Uniform: 0.739751 +/- 0.00627652
Exponential: 0.745411 +/- 0.00175903
At N = 10000:
Uniform: 0.748328 +/- 0.0020065
Exponential: 0.746198 +/- 0.000553959
At N = 100000:
Uniform: 0.747267 +/- 0.000635817
Exponential: 0.746977 +/- 0.00017412
At N = 1000000:
Uniform: 0.747092 +/- 0.000200955
Exponential: 0.746832 +/- 5.49985e-05
At N = 10000000:
Uniform: 0.746947 +/- 6.35594e-05
Exponential: 0.746829 +/- 1.73977e-05
At N = 100000000:
Uniform: 0.746829 +/- 2.00995e-05
Exponential: 0.746825 +/- 5.50143e-06

 * */
