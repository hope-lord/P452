#include "matrix_utility.cpp"

void gauss_jordan(MATRIX &ar, 
					VECTOR &br=null_vector,
					MATRIX &cr=null_matrix){
	int LIMIT=ar.size();

	for(int i=0;i<LIMIT;i++){
		if (ar[i][i]==0){
			bool c=partial_pivot(ar,br,cr);
			if(!c) {cout<<"WARNING: DETERMINANT IS 0.\n";return;}
		}
		D_t pivot=ar[i][i];
		if (br!=null_vector) br[i]/=pivot;
		for(int j=0;j<LIMIT;j++){
			ar[i][j]/=pivot;
			if (cr!=null_matrix) cr[i][j]/=pivot;
		}
		for(int j=0;j<LIMIT;j++){
			if(i!=j and ar[j][i]!=0){
				D_t fact=ar[j][i];
				for(int t=0;t<LIMIT;t++){
					ar[j][t]-=fact*ar[i][t];
					if (cr!=null_matrix) cr[j][t]-=fact*cr[i][t];
				}
				if(br!=null_vector)br[j]-=fact*br[i];
			}
		}
	}
}

/*-------------------------------------------*/


int LU(MATRIX &ar,VECTOR &br=null_vector,MATRIX &cr=null_matrix){
	if (ar.size()!=ar[0].size()){cout<<"The operator matrix is not square one\n";return -1;}
	if(br!=null_vector) {
		if (br.size()!=ar.size()){cout<<"The Dimensions are not equal.\n";return -1;}
	}
	if(cr!=null_matrix) {
		if (cr.size()!=ar.size()){cout<<"The Dimensions are not equal.\n";return -1;}
		if(cr[0].size()!=ar.size()){cout<<"The Dimensions are not equal.\n";return -1;}
	}
	int LIMIT=ar.size();
	int row_rotations=0;
	int i=0;int s=0;bool t=false;
    
	while(i<LIMIT){
		if(!t) s=i;
		VECTOR new_row;
		new_row=ar[i];
		int j=0;
		while(j<LIMIT){
			if(i>j){
				D_t sum=0;
				for(int k=0;k<j;k++)sum+=ar[i][k]*ar[k][j];
				ar[i][j]=(ar[i][j]-sum)/ar[j][j];
				j++;
			}
			else if(i==j){
				D_t sum=0;
				for(int k=0;k<j;k++)sum+=ar[i][k]*ar[k][j];
				ar[i][i]-=sum;
				if (ar[i][i]!=0){
					j++;t=false;
				}
				else{
					j=LIMIT;
					t=true;
					ar[i]=new_row;
				}
			}
			else{
				D_t sum=0;
				for(int k=0;k<i;k++)sum+=ar[i][k]*ar[k][j];
				ar[i][j]-=sum;
				j++;
			}
		}
		if(!t) i++;
		else{
			s++;
			if(s==LIMIT) return -1;
			partial_pivot_swap(ar,i,s,br,cr);
			row_rotations++;
		}
	}
	return row_rotations;
}

/*--------------------------------------------*/

//br or cr is not requird for calculation of determinant but br will have to undergo partial_pivoting
D_t lu_and_det(MATRIX &ar,VECTOR &br=null_vector,MATRIX &cr=null_matrix){
	int LIMIT=ar.size();
	int p=LU(ar,br,cr);
	if(p==-1) return 0;
	D_t mult=1;
	for(int i=0;i<LIMIT;i++) mult*=ar[i][i];
	if(p%2==0) return mult;
	else return -mult;
}

/*-----------------------------------------------------------*/

void forward_substitution(MATRIX &ar,VECTOR &br){
	int LIMIT=ar.size();
    for(int i=0;i<LIMIT;i++){
        D_t sum=0;
        for(int j=0;j<i;j++) sum+=ar[i][j]*br[j];
        br[i]-=sum;
	}
}

/*---------------------------------------------------------------*/

void backward_substitution(MATRIX &ar,VECTOR &br){
	int LIMIT=ar.size();
    int i = LIMIT - 1;
    while (i > -1){
        D_t sum = 0;
        for(int j=i+1;j<LIMIT;j++) sum += ar[i][j] * br[j];
        br[i] -= sum;
        br[i] /= ar[i][i];
        i--;
	}
}

/*--------------------------------------------------------------------*/

VECTOR solve(MATRIX &ar,VECTOR &br){
	D_t s=lu_and_det(ar,br);
	if(s==0)cout<<"Determinant of the operator is zero. Hence no unique solution possible\n";
	else{
		forward_substitution(ar,br);
		backward_substitution(ar,br);
	}
	return br;
}

/*-----------------------------------------------------------------------*/

MATRIX inv_lu(MATRIX &ar){
	int LIMIT=ar.size();
    MATRIX inv;
	for(int i=0;i<LIMIT;i++) {
		VECTOR m;
		for(int j=0;j<LIMIT;j++){
			if(i==j) m.push_back(1);
			else m.push_back(0);
		}
		inv.push_back(m);
	}
	D_t s=lu_and_det(ar,null_vector,inv);
	if(s==0) cout<<"Determinant of the operator is zero. Hence no unique solution possible\n";
	else{
	    for(int i=0;i<LIMIT;i++){
			VECTOR column;
	    	for(int j=0;j<LIMIT;j++)column.push_back(inv[j][i]);
    		forward_substitution(ar,column);
            backward_substitution(ar,column);
		    for(int j=0;j<LIMIT;j++)inv[j][i]=column[j];
	    }
    }
	return inv;
}

/*--------------------------------------------------------------------*/
MATRIX inv(MATRIX &ar){
	int LIMIT=ar.size();
    MATRIX inv;
	for(int i=0;i<LIMIT;i++) {
		VECTOR m;
		for(int j=0;j<LIMIT;j++){
			if(i==j) m.push_back(1);
			else m.push_back(0);
		}
		inv.push_back(m);
	}
	gauss_jordan(ar,null_vector,inv);
	return inv;
}

/*---------------------------------------------------------------------*/



/*---------------------Gauss Seidel Method-----------------------------*/
bool __Gauss_Seidel__helper__(MATRIX &A,VECTOR &b, VECTOR &x, double &err ){
	int N = A.size();
	err = 0;
	for(int i=0;i<N;i++) {
		if (A[i][i]==0) {
			fprintf(stderr,"The %dth diagonal is 0.\n",i);
			return false;
		}
		double a = x[i];
		double sum =  0;
		for(int j=0;j<N;j++) {
			if(j!=i) sum += A[i][j]*x[j];
		}
		x[i] = (b[i] - sum)/A[i][i];
		err += (x[i]-a)*(x[i]-a);
	}
	err = sqrt(err);
	return true;
}
bool Gauss_Seidel(MATRIX &A,VECTOR &x, VECTOR &b,const double eps=1e-6){	
	double err = eps*100; // there will be iteration as long as the error is not less than eps
		for(int i=0;i<5000;i++){
		if(!__Gauss_Seidel__helper__(A, b, x, err)) return false;
			if(err<=eps) {
			/* cout<<"Converged at "<<i<<endl; */
				return true;
		}
	}
	fprintf(stderr,"Did not converge with 5000 iterations.\n");
	return false;
}

/*-------------------Jacobi Method------------------------------------*/
bool __Jacobi_helper__(MATRIX &A,VECTOR &b, VECTOR &x, double &err ){
	int N = A.size();
	VECTOR x0(x);
	x0 = matrix_mult(A,x0);
	
	err = 0;
	for(int i=0;i<N;i++) {
		if (A[i][i]==0) {
			fprintf(stderr,"The %dth diagonal is 0.\n",i);
			return false;
		}
		double a = x[i];
		x[i] = (b[i] - x0[i]+A[i][i]*x[i])/A[i][i];
		err += (x[i]-a)*(x[i]-a);
	}
	err = sqrt(err);
	return true;
}
bool Jacobi(MATRIX &A,VECTOR &x, VECTOR &b,const double eps=1e-6){	
	double err = eps*100; // there will be iteration as long as the error is not less than eps
		for(int i=0;i<5000;i++){
		if(!__Jacobi_helper__(A,b,x,err)) return false;
			if(err<=eps) {
			/* cout<<"Converged at "<<i<<endl; */
				return true;
		}
	}
	fprintf(stderr,"Did not converge with 5000 iterations.\n");
	return false;
}

/*---------------------------------Conjugate GRadient----------------*/
bool Conjugate_Gradient(MATRIX &A,VECTOR &x0,VECTOR &b,double err=1e-8){
    size_t N=x0.size();
    VECTOR r(N);
    double alpha,beta,rr,dad,rr_new,r_norm;
    VECTOR Ax = matrix_mult(A,x0);
    for(size_t i=0;i<N;i++) r[i]=b[i]-Ax[i];
    VECTOR d(r);  //copy initial residual to d
    rr = dot(r,r);
    r_norm = sqrt(rr);
    if(r_norm<err) return true;
    size_t i= 0;
    while(i<N){
        VECTOR Ad = matrix_mult(A,d);
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


