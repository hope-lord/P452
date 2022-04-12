# include "package_cpp/UTILITY.h"

# include "package_cpp/Linear.cpp"
# include "package_cpp/data_stats.cpp"
/* # include "cond.cpp" */

double phi(double x,int j){
      if (j==0) return 1.0;
	  else if (j==1) return 2*x-1;
	  else if (j==2) return 8*x*x-8*x+1;
	  else 	return 32*x*x*x-48*x*x+18*x-1;
}

void Chebyshev_poly_fit(VECTOR &xx,VECTOR &yy,VECTOR & ww){
	const int n = 3;

    MATRIX A(n+1,VECTOR(n+1,0));
    const int N = xx.size();
    VECTOR x_vec(2*n+1);
    VECTOR yx_vec(n+1);
	double S_yy;
    // read data files
    for(int i=0;i<N;i++){
		if (ww[i]==0) ww[i] = 1;
        else ww[i] = 1/ww[i]/ww[i];
        /* sum_x += x;sum_y += y;sum_xy += x*y;sum_xx += x*x;sum_yy += y*y; */
		S_yy += yy[i]*yy[i]*ww[i];
        for(int j=0;j<n+1;j++){
            // x_vec[j] = pow(xx[i],j);
            yx_vec[j] += yy[i]*phi(xx[i],j)*ww[i];
            for(int k=j;k<n+1;k++){
                A[j][k] += ww[i]*phi(xx[i],k)*phi(xx[i],j);
                if(k!=j) A[k][j] = A[j][k];   
                // display_matrix(A); 
                // cout<<endl;
            }
            
        }
    }
    
	double S_xx = A[1][1];
	double S_xy = yx_vec[1];
    // solve linear equation
    A = inv(A);
    VECTOR coeff = matrix_mult(A,yx_vec);
    VECTOR coeff_err(n+1);
    for(int i=0;i<n+1;i++)
        coeff_err[i] = sqrt(A[i][i]);

    double chi2 = 0;
    for(int i=0;i<N;i++){
        double y_fit = 0;
		for(int l = 0;l<n+1;l++) y_fit += coeff[l]*phi(xx[i],l);
		/* polynomial(xx[i],coeff); */
        chi2 += (yy[i]-y_fit)*(yy[i]-y_fit)*ww[i];
    }
    double r = S_xy/S_xx/S_yy;

    /* fprintf(out,"Pearson's r = %lg\n",r); */
    fprintf(stdout,"𝜒2/ν = %lg/%d\n",chi2,N-n-1);
    fprintf(stdout,"Coefficients:\n");
    for(int i=0;i<n+1;i++)
        fprintf(stdout,"a%d = %lg \n",i,coeff[i]);

    // Covariance matrix
    fprintf(stdout,"Covariance matrix:\n");
    for(int i=0;i<n+1;i++){
        for(int j=0;j<n+1;j++){
            fprintf(stdout,"%lg ",A[i][j]);
        }
        fprintf(stdout,"\n");
}

}

int main(){
	// The Condition numbers are calculated using L-1 NORM 
	int N = 21;
	FILE * file = fopen("assign2fit.txt","r");
	VECTOR xx(N),yy(N),ww(N);
	for(int i=0;i<N;i++){
		fscanf(file,"%lg %lg",&xx[i],&yy[i]);
		ww[i] = 1;
	}
	cout << "For Simple Polynomial Fit :\n";
	cout << "Condition No = "<< 22030.866<<endl;
	// For normal polynomial fit
	poly_fit(3,xx,yy,ww);
	// For Chebyshev polynomial fit
	cout << endl;
	cout << "For Chebyshev Polynomial Fit :\n";
	cout << "Condition No = "<< 4.798<<endl;
	Chebyshev_poly_fit(xx,yy,ww);

	return 0;
}

/*OUTPUT::  */
/*
For Simple Polynomial Fit :
Condition No = 22030.9
𝜒2/ν = 0.0371505/17
Coefficients:
a0 = 0.574659
a1 = 4.72586
a2 = -11.1282
a3 = 7.66868
Covariance matrix:
0.544043 -3.96041 7.71692 -4.39174 
-3.96041 42.8691 -97.3558 60.1489 
7.71692 -97.3558 238.277 -154.096 
-4.39174 60.1489 -154.096 102.731 

For Chebyshev Polynomial Fit :
Condition No = 4.798
𝜒2/ν = 0.0371505/17
Coefficients:
a0 = 1.16097 
a1 = 0.393514 
a2 = 0.0468498 
a3 = 0.239646 
Covariance matrix:
0.055544 8.08007e-18 0.0297186 1.76042e-17 
8.08007e-18 0.143456 1.19027e-18 0.0369189 
0.0297186 1.19027e-18 0.111445 9.05639e-19 
1.76042e-17 0.0369189 9.05639e-19 0.100323 
*/

/*Since Chebyshev Polynomial Fitting has less condition number, The change is covariance matrix is less prone due to perturbation than simple polynomial fitting.*/
