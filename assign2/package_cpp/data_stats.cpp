// JACKNIFE METHOD
/* #include "package_cpp/UTILITY.h" */
/* #include "package_cpp/Linear.cpp" */
/* #include "package_cpp/io.cpp" */
double Jacknife(general_function &func,VECTOR &y){
    int N = y.size();
    VECTOR f(N);
    // double f_avg;    
    double y_avg = 0;
    double sigma2_jk_y = 0;
    for(int i=0;i<N;i++){ 
        y_avg += y[i]/N;
        // f_avg += func.function(y[i],func.paras)/N;    
    }
    for(int i=0;i<N;i++){
        y[i] =  (N*y_avg-y[i])/N;
        sigma2_jk_y += (y_avg-y[i])*(y_avg-y[i])/N;
        // f[i] = N*func.function(y[i],func.paras)/N;
    }
    return sigma2_jk_y;
}

double polynomial(double x,VECTOR &list){
    int n=list.size();
    double sum = 0;
    for(int i=0;i<n;i++) sum += list[i] * pow(x ,i);
    return sum;
}


// n order polynomial fit with N data points
void poly_fit(const int n, VECTOR &xx,VECTOR &yy,VECTOR & ww){
    MATRIX A(n+1,VECTOR(n+1,0));
    const int N = xx.size();
    VECTOR x_vec(2*n+1);
    VECTOR yx_vec(n+1);
	double S_yy;
    // read data files
    for(int i=0;i<N;i++){
		if (ww[i]==0) ww[i] = 1;
        else ww[i] = 1/ww[i]/ww[i];
    // double sum_x =0, sum_y = 0, sum_xy = 0, sum_xx = 0, sum_yy = 0;
        /* sum_x += x;sum_y += y;sum_xy += x*y;sum_xx += x*x;sum_yy += y*y; */
		S_yy += yy[i]*yy[i]*ww[i];
        for(int j=0;j<n+1;j++){
            // x_vec[j] = pow(xx[i],j);
            yx_vec[j] += yy[i]*pow(xx[i],j)*ww[i];
            for(int k=j;k<n+1;k++){
                A[j][k] += ww[i]*pow(xx[i],k+j);
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
        double y_fit = polynomial(xx[i],coeff);
        chi2 += (yy[i]-y_fit)*(yy[i]-y_fit)*ww[i];
    }
    double r = S_xy/S_xx/S_yy;

    /* fprintf(out,"Pearson's r = %lg\n",r); */
    fprintf(stdout,"𝜒2/ν = %lg/%d\n",chi2,N-n-1);
    fprintf(stdout,"Coefficients:\n");
    for(int i=0;i<n+1;i++)
        fprintf(stdout,"a%d = %lg\n",i,coeff[i]);

    // Covariance matrix
    fprintf(stdout,"Covariance matrix:\n");
    for(int i=0;i<n+1;i++){
        for(int j=0;j<n+1;j++){
            fprintf(stdout,"%lg ",A[i][j]);
        }
        fprintf(stdout,"\n");
    }
}
