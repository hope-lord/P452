# include "package_cpp/UTILITY.h"

# include "package_cpp/Linear.cpp"

# include "package_cpp/io.cpp"

int main(void){
	int N = 6;
	FILE *file = fopen("q1.txt","r");
	// Define vector and matrix
	VECTOR b(N);
	MATRIX A(N,b);
	// read them from the file
	for(int i = 0;i<N;i++){
		for(int j = 0;j<N;j++){
			fscanf(file,"%lf",&A[i][j]);
		}
	}
	for(int i = 0;i<N;i++){
		fscanf(file,"%lf",&b[i]);
	}
	fclose(file);
	// copy matrix A and vector b since after beign used at Gauss Jordan it will be identity
	MATRIX A1(A);
	VECTOR b1(b);
	// Gauss Jordan Method
	cout <<"Solution by Gauss-Jordan Method"<<endl;
	gauss_jordan(A,b);
	display_vector(b);
	
	cout<<endl;
	//LU Method
	cout <<"Solution by LU Method"<<endl;
	solve(A1,b1);
	display_vector(b1);
	cout<<endl;
	return 0;
}

/* Solution by Gauss-Jordan Method */
/* [-1.76182 0.896228 4.05193 -1.61713 2.04191 0.151832] */
 
/* Solution by LU Method */
/* [-1.76182 0.896228 4.05193 -1.61713 2.04191 0.151832] */

