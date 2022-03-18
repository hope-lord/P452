/*-------------------------------------Sparse Mult--------------------*/
VECTOR matrix_mult(double (*operate_A)(int,VECTOR&),VECTOR& x0){
	int N = x0.size();
	VECTOR mult(N);
	for(int i=0;i<N;i++){
		mult[i] = operate_A(i,x0);
	}
	return mult;
}

/* --------------------------------Conjugate Gradient Method------------------------------------------------ */

bool Conjugate_Gradient(double (*func)(int,VECTOR &),VECTOR &x0,VECTOR &b,double err=1e-8){
    size_t N=x0.size();
    VECTOR r(N);
    double alpha,beta,rr,dad,rr_new,r_norm;
    VECTOR Ax = matrix_mult(func,x0);
    for(size_t i=0;i<N;i++) r[i]=b[i]-Ax[i];
    VECTOR d(r);  //copy initial residual to d
    rr = dot(r,r);
    r_norm = sqrt(rr);
    if(r_norm<err) return true;
    size_t i= 0;
    while(i<N){
        VECTOR Ad = matrix_mult(func,d);
        dad = dot(d,Ad);
        alpha = rr/dad;
        //beta = rr/rr_prev;
        for(size_t j=0;j<N;j++) {
            x0[j] += alpha*d[j];
            r[j] -= alpha*Ad[j];
        }
        rr_new = dot(r,r);
        r_norm = sqrt(rr_new);
        if(r_norm<err) i=N;
        beta = rr_new/rr;
        for(size_t j=0;j<N;j++) d[j] = r[j] + beta*d[j];
        i++;
        rr = rr_new;
    }

    return (r_norm<err);//true;
}

/* --------------------Gauss Sedel Method ------------------------- */



bool __Gauss_Seidel__helper__(double (*diag)(int),double (*operate_A)(int,VECTOR&),VECTOR &b, VECTOR &x, double &err ){
	int N = b.size();
	err = 0;
	for(int i=0;i<N;i++) {
		if (diag(i)==0) { // diag function returns the i the diagonal of A
			fprintf(stderr,"The %dth diagonal is 0.\n",i);
			return false;
		}
		double a = x[i];
		x[i] = (b[i] - operate_A(i,x)+diag(i)*x[i])/diag(i); // x[i] = (b[i]-(L+U)x0 [i])/A_ii
		err += (x[i]-a)*(x[i]-a); // error addition
	}
	err = sqrt(err);
	return true;
}
bool Gauss_Seidel(double (*diag)(int),double (*operate_A)(int,VECTOR&),VECTOR &x, VECTOR &b,const double eps=1e-6){	
	double err = eps*100; // there will be iteration as long as the error is not less than eps
		for(int i=0;i<5000;i++){ // iterate N times
		if(!__Gauss_Seidel__helper__(diag,operate_A,b,x,err)) return false;
			if(err<=eps) {
			/* cout<<"Converged at "<<i<<endl; */
				return true;
		}
	}
	fprintf(stderr,"Did not converge with 5000 iterations.\n");
	return false;
}


/* ---------------------Jacobi Method ------------------------- */

bool __Jacobi_helper__(double (*diag)(int),double (*operate_A)(int,VECTOR&),VECTOR &b, VECTOR &x, double &err ){
	int N = b.size();
	VECTOR x0(x);
	x0 = matrix_mult(operate_A,x0);
	
	err = 0;
	for(int i=0;i<N;i++) {
		if (diag(i)==0) {
			fprintf(stderr,"The %dth diagonal is 0.\n",i);
			return false;
		}
		double a = x[i];
		x[i] = (b[i] - x0[i]+diag(i)*x[i])/diag(i);
		err += (x[i]-a)*(x[i]-a);
	}
	err = sqrt(err);
	return true;
}
bool Jacobi(double (*diag)(int),double (*operate_A)(int,VECTOR&),VECTOR &x, VECTOR &b,const double eps=1e-6){	
	double err = eps*100; // there will be iteration as long as the error is not less than eps
		for(int i=0;i<5000;i++){
		if(!__Jacobi_helper__(diag,operate_A,b,x,err)) return false;
			if(err<=eps) {
			cout<<"Converged at "<<i<<endl;
				return true;
		}
	}
	fprintf(stderr,"Did not converge with 5000 iterations.\n");
	return false;
}

