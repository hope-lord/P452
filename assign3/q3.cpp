// Solve Laplace equation with Dirichlet boundary conditions 
// For 2D system
# include "package_cpp/UTILITY.h"
/* # include "package_cpp/sparse_Linear.cpp" */
# include "package_cpp/matrix_utility.cpp"


class Solve_Laplace_2D{
	
	int nx, ny, n;
	double x0,x1,y0,y1,h;
    double (*BoundaryCondition)(double ,double ) = NULL;

	VECTOR phi,B;
	public:
    Solve_Laplace_2D(double(*f)(double,double),double x0,double x1,double y0,double y1,double h){
		this->x0 = x0;
		this->x1 = x1;
		this->y0 = y0;
		this->y1 = y1;
		this->h = h;
		this->nx = (x1-x0)/h -1;
		this->ny = (y1-y0)/h -1;
		this->n = nx*ny;
		/* cout << nx << " " << ny << endl; */
		this->BoundaryCondition = f;
		phi.resize(n);
		for(int i = 0; i < n; i++)
			phi[i] = random_double(-1,1);
		
		B.resize(n,0); // Boundry Condition 
		for(int i = 0; i < nx; i++)
		    B[i] -= BoundaryCondition(x0+i*h+h,y0);
		for(int i = 0; i < nx; i++)
		    B[i+nx*(ny-1)] -= BoundaryCondition(x0+i*h+h,y1);
		for(int i = 0; i < ny; i++)
			B[i*nx] -= BoundaryCondition(x0,y0+i*h+h);
		for(int i = 0; i < ny; i++)
			B[i*nx+nx-1] -= BoundaryCondition(x1,y0+i*h+h);
	}
	
	void solve(){
		Conjugate_Gradient(B);
	}
	
	void Possion_Function(double(*rho)(double,double)){
		for(int i = 0; i < nx; i++)
			for(int j = 0; j < ny; j++)
				B[i *nx+j] += h*h*rho(x0+i*h+h,y0+j*h+h);
	}

	double Laplace_Operator(int i, VECTOR& vec){
		int x = i%nx;
		int y = i/nx;

		double sum = -4*vec[i] ;
		if(x>0) sum += vec[i-1];
		/* else  sum += BoundaryCondition(x0,h*(y+1)+y0); */
		if(x<nx-1) sum += vec[i+1];
		/* else  sum += BoundaryCondition(x1,h*(y+1)+y0); */
		if(y>0) sum += vec[i-nx];
		/* else  sum += BoundaryCondition(h*(x+1)+x0,y0); */
		if(y<ny-1) sum += vec[i+nx];
		/* else  sum += BoundaryCondition(h*(x+1)+x0,y1); */
		return sum;
   }

	VECTOR get_phi(){
		return phi;
	}

	void print_to_file (string filename){
		FILE *file = fopen(filename.c_str(),"w");
		for ( int i = 0; i < ny+2; i++ )
			fprintf(file,"%lf %lf %lf\n",x0,y0+i*h,BoundaryCondition(x0,y0+i*h));
		for(int i=0;i<nx;i++){
			fprintf(file,"%lf %lf %lf\n",x0+i*h+h,y0,BoundaryCondition(x0+i*h+h,y0));
			for(int j=0;j<ny;j++){
				fprintf(file,"%lf %lf %lf\n",x0+i*h+h,y0+j*h+h,phi[j*nx+i]);
			}
			fprintf(file,"%lf %lf %lf\n",x0+i*h+h,y1,BoundaryCondition(x0+i*h+h,y1));
		}

		/* for ( int i = 0; i < nx+2; i++ ) */
		/* 	fprintf(file,"%lf %lf %lf\n",x0+i*h,y0,BoundaryCondition(x0+i*h,y0)); */
		/* for ( int i = 0; i < nx+2; i++ ) */
		/* 	fprintf(file,"%lf %lf %lf\n",x0+i*h,y1,BoundaryCondition(x0+i*h,y1)); */
		/* for ( int i = 0; i < ny+2; i++ ) */
			/* fprintf(file,"%lf %lf %lf\n",x0,y0+i*h,BoundaryCondition(x0,y0+i*h)); */
		for ( int i = 0; i < ny+2; i++ )
			fprintf(file,"%lf %lf %lf\n",x1,y0+i*h,BoundaryCondition(x1,y0+i*h));

		fclose(file);
	}

	VECTOR matrix_mult(VECTOR& vec){
		VECTOR mult(n);
		for(int i=0;i<n;i++){
		mult[i] =  Laplace_Operator(i,vec);
		}
		return mult;
	}

/* --------------------------------Conjugate Gradient Method------------------------------------------------ */

	bool Conjugate_Gradient(VECTOR &b,double err=1e-8){
		VECTOR r(n);
		double alpha,beta,rr,dad,rr_new,r_norm;
		VECTOR Ax = matrix_mult(phi);
		for(size_t i=0;i<n;i++) r[i]=b[i]-Ax[i];
		VECTOR d(r);  //copy initial residual to d
		rr = dot(r,r);
		r_norm = sqrt(rr);
		if(r_norm<err) return true;
		size_t i= 0;
		while(i<n and r_norm>err){
        VECTOR Ad = matrix_mult(d);
        dad = dot(d,Ad);
        alpha = rr/dad;
        //beta = rr/rr_prev;
        for(size_t j=0;j<n;j++) {
            phi[j] += alpha*d[j];
            r[j] -= alpha*Ad[j];
        }
        rr_new = dot(r,r);
        r_norm = sqrt(rr_new);
        if(r_norm<err) i=n;
        beta = rr_new/rr;
        for(size_t j=0;j<n;j++) d[j] = r[j] + beta*d[j];
        i++;
        rr = rr_new;
    	}
		return (r_norm<err);//true;
	}
	
	~Solve_Laplace_2D(){
		phi.clear();
		B.clear();
	}

};

double f(double x,double y){
	if (y==0) return 1;
	else return 0;
	/* return 0; */
}

double rho(double x,double y){
	/* return -1; */
	return 0;
}

int main(){
	Solve_Laplace_2D *solver= new Solve_Laplace_2D(&f,0,1,0,1,0.05);
	solver->Possion_Function(&rho);
	solver->solve();
	solver->print_to_file("q3_out.txt");
	delete solver;
	return 0;
}


