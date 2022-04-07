# include "package_cpp/UTILITY.h"

#include "package_cpp/matrix_utility.cpp"

# include "package_cpp/sparse_Linear.cpp"

# include "package_cpp/io.cpp"

const int N = 20;
const int m = 0.2;


double operate_A(int i,VECTOR & x)
{
	double sum = 0;
	if (i-2<0) sum += x[i-2+N];
	else sum += x[i-2];
	if (i+2>N-1) sum += x[i+2-N];
	else sum += x[i+2];
	return 0.5*sum  + (1+m*m)*x[i];
}

int main(){
	VECTOR x(N,1);
	FILE* file = fopen("q3_out.txt","w");
	MATRIX A;
	for (int i=0;i<N;i++){
		VECTOR b(N,0);
		b[i] = 1;
		VECTOR x(N,1);
		Conjugate_Gradient(&operate_A,x,b,1e-6);
		for (int j =0;j<N;j++) fprintf(file,"%lf ",x[j]);
		fprintf(file,"\n");
		A.push_back(x);
	}
	fclose(file);
	VECTOR c = matrix_mult(&operate_A, A[0]);
	display_vector(c);
	return 0;
}
