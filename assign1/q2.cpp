# include "package_cpp/UTILITY.h"

# include "package_cpp/Linear.cpp"

# include "package_cpp/io.cpp"


bool Gauss_Seidel_print(string file_name,MATRIX &A, VECTOR &x, VECTOR &b,const double eps);
bool Jacobi_print(string file_name,MATRIX &A, VECTOR &x, VECTOR &b,const double eps);


int main(void){
	int N = 6;
	double eps = 1e-6;
	FILE *file = fopen("q2.txt","r");
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
	// copy A and b
	MATRIX A1(A);
	VECTOR b1(b);
	// LU Method
	cout <<"Solution by LU Method"<<endl;
	VECTOR x = solve(A1,b1); // This function uses LU Decompostion
	display_vector(x);
	cout<<endl;
	// Jacobi Method
	cout <<"Solution by Jacobi Method"<<endl;
	VECTOR x0(N,1); // initial guess
	Jacobi(A, x0, b,eps);
	display_vector(x0);
	cout<<endl;

	// Inverse of A
	// /// / //
	// Using Gauss-Seidel Method
	cout <<"Inverse by Gauss-Seidel Method"<<endl;
	MATRIX Ainv_Gauss(N,VECTOR(N));
	for(int i = 0;i<N-1;i++){
		VECTOR m(N,0);
		VECTOR x1(N,1); //guesss vector
		m[i] = 1; // ith column of identity
		Gauss_Seidel(A, x1, m,eps);
		for(int j = 0;j<N;j++){
			Ainv_Gauss[j][i] = x1[j];
		}
    }
	// Do the last column inversion and plot residual vs iteration to the file
	VECTOR m(N,0);
	m[N-1] = 1; // last column which will be plotted 
	VECTOR x1(N,1);
	Gauss_Seidel_print("q2_Gauss_Seidal_iter.txt",A, x1, m,eps);
	for(int j = 0;j<N;j++){
		Ainv_Gauss[j][N-1] = x1[j];
	}
	display_matrix(Ainv_Gauss);
	cout<<endl;
	// Using Jacobi Method
	cout <<"Inverse by Jacobi Method"<<endl;
	MATRIX Ainv_Jacobi(N,VECTOR(N));
	for(int i = 0;i<N-1;i++){
		VECTOR m1(N,0);
		VECTOR x2(N,1);
		m1[i] = 1;
		Jacobi(A, x2, m1,eps);
		for(int j = 0;j<N;j++){
			Ainv_Jacobi[j][i] = x2[j];
		}
    }
	// Do the last column inversion and plot residual vs iteration to the file
	VECTOR m1(N,0);
	m1[N-1] = 1;
	VECTOR x2(N,1);
	Jacobi_print("q2_Jacobi_iter.txt",A, x2, m1,eps);
	for(int j = 0;j<N;j++){
		Ainv_Jacobi[j][N-1] = x2[j];
	}
	display_matrix(Ainv_Jacobi);
	cout<<endl;

	return 0;
}

// this function is same as Gauss_Seidel function with extra feature of printing residual vs iteration no to a file
bool Gauss_Seidel_print(string file_name,MATRIX &A, VECTOR &x, VECTOR &b,const double eps){
	FILE *file = fopen(&file_name[0],"w");
	int N = A.size();
	double err = eps*100; // there will be iteration as long as the error is not less than eps
		for(int i=0;i<5000;i++){
			VECTOR r = matrix_mult(A,x);
			for(int j=0;j<N;j++) r[j]-=b[j];
			double r_norm = dot(r,r);
			fprintf(file,"%d %lf\n",i,sqrt(r_norm));
		if(!__Gauss_Seidel__helper__(A, b, x, err)) return false;
			if(err<=eps) {
			/* cout<<"Converged at "<<i<<endl; */
				return true;
		}
	}
	fprintf(stderr,"Did not converge with 5000 iterations.\n");
	fclose(file);
	return false;
}

// this function is same as Jacobi function with extra feature of printing residual vs iteration no to a file
bool Jacobi_print(string file_name,MATRIX &A, VECTOR &x, VECTOR &b,const double eps){
	FILE *file = fopen(&file_name[0],"w");
	int N = A.size();
	double err = eps*100; // there will be iteration as long as the error is not less than eps
		for(int i=0;i<5000;i++){
			VECTOR r = matrix_mult(A,x);
			for(int j=0;j<N;j++) r[j]-=b[j];
			double r_norm = dot(r,r);
			fprintf(file,"%d %lf\n",i,sqrt(r_norm));
		if(!__Jacobi_helper__(A, b, x, err)) return false;
			if(err<=eps) {
			/* cout<<"Converged at "<<i<<endl; */
				return true;
		}
	}
	fprintf(stderr,"Did not converge with 5000 iterations.\n");
	fclose(file);
	return false;
}
/* Solution by LU Method */
/* [-0.333312 0.333349 1.00001 -0.666605 3.00433e-05 0.666677] */

/* Solution by Jacobi Method */
/* [-0.33331 0.33335 1.00001 -0.666602 3.14741e-05 0.666677] */

/* Inverse by Gauss-Seidel Method */
/* [[0.9351 0.8701 0.2597 0.2078 0.4156 0.1688] */
/* [0.2900 0.5801 0.1732 0.1385 0.2771 0.1126] */
/* [0.0866 0.1732 0.3203 0.0563 0.1126 0.1082] */
/* [0.2078 0.4156 0.1688 0.9351 0.8701 0.2597] */
/* [0.1385 0.2771 0.1126 0.2900 0.5801 0.1732] */
/* [0.0563 0.1126 0.1082 0.0866 0.1732 0.3203]] */

/* Inverse by Jacobi Method */
/* [[0.9351 0.8701 0.2597 0.2078 0.4156 0.1688] */
/* [0.2900 0.5801 0.1732 0.1385 0.2771 0.1126] */
/* [0.0866 0.1732 0.3203 0.0563 0.1126 0.1082] */
/* [0.2078 0.4156 0.1688 0.9351 0.8701 0.2597] */
/* [0.1385 0.2771 0.1126 0.2900 0.5801 0.1732] */
/* [0.0563 0.1126 0.1082 0.0866 0.1732 0.3203]] */



