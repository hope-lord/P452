# include "package_cpp/UTILITY.h"

#include "package_cpp/matrix_utility.cpp"

# include "package_cpp/sparse_Linear.cpp"

# include "package_cpp/io.cpp"

bool Conjugate_Gradient_print(string file_name,double (*func)(int,VECTOR &),VECTOR &x0,VECTOR &b,double err=1e-8);


const int n = 20;
const int N = 400; // 20 X 20 lattice sites

const double m2 = 0.04;


double operate_A(int i,VECTOR & vec){
	int x = i%n;
	int y = i/n;
	double sum = 0;
	// Periodic boundry condition n = 0
	// // For x component
	if (x==0) sum += 0.5*vec[y*n+n-1];
	else sum+= 0.5*vec[y*n+x-1];

	if (x==n-1) sum += 0.5*vec[y*n];
	else sum += 0.5*vec[y*n+x+1];

	// // For y component
	if (y==0) sum += 0.5*vec[n*n-n+x];
	else sum += 0.5*vec[y*n-n+x];

	if(y==n-1) sum += 0.5*vec[x];
	else sum += 0.5*vec[y*n+n+x]; 

	return sum+ (m2-2)*vec[i];
}



int main(){
	srand(5019);
	FILE* file = fopen("q3_out.txt","w"); // print inverse matrix to this file
	VECTOR x0(N,1); // random guess field for Choose all values to be 1 for first column
	VECTOR x1(N,1); // random guess field for Choose all values to be 1 for second column
	VECTOR b0(N,0); 
	VECTOR b1(N,0);
	b0[0] = 1; // first column of identity matrix
	b1[1] = 1; // second column of identity matrix

	Conjugate_Gradient_print("q3_residual_col_0.txt",&operate_A,x0,b0,1e-6);
	Conjugate_Gradient_print("q3_residual_col_1.txt",&operate_A,x1,b1,1e-6);

    // print the first two columns of the inverse matrix to q3_out.txt file
	for (int j =0;j<N;j++) fprintf(file,"%lf %lf\n",x0[j],x1[j]);
	fprintf(file,"\n");
	fclose(file);
	return 0;
}


bool Conjugate_Gradient_print(string file_name,double (*func)(int,VECTOR &),VECTOR &x0,VECTOR &b,double err){
	FILE* file = fopen(&file_name[0],"w");
    size_t N=x0.size();
    VECTOR r(N);
    double alpha,beta,rr,dad,rr_new,r_norm;
    VECTOR Ax = matrix_mult(func,x0);
    for(size_t i=0;i<N;i++) r[i]=b[i]-Ax[i]; // residual vector
    VECTOR d(r);  //copy initial residual to d
    rr = dot(r,r);
    r_norm = sqrt(rr);
    if(r_norm<err) return true; // return incase guess is correct
    size_t i= 0;
    while(i<N and r_norm>err){
        VECTOR Ad = matrix_mult(func,d);
        dad = dot(d,Ad);
        alpha = rr/dad;
        for(size_t j=0;j<N;j++) {
            x0[j] += alpha*d[j]; // update x0
            r[j] -= alpha*Ad[j]; // update residual
        }
        rr_new = dot(r,r);
        r_norm = sqrt(rr_new); //update r_norm
		fprintf(file,"%ld %lg\n",i+1,r_norm);
        if(r_norm<err) i=N;
        beta = rr_new/rr;
        for(size_t j=0;j<N;j++) d[j] = r[j] + beta*d[j]; // update d
        i++;
        rr = rr_new; // update residual norm 
    }
    return (r_norm<err);//true;
}
