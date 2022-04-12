
void swap(D_t &a,D_t &b){
	D_t t;
	a=b;
	b=t;
}

/*---------------------------------------------------*/

//This partial pivot is useful in solving linear equations
void partial_pivot_swap(MATRIX &ar,int t,int s,VECTOR &br=null_vector,MATRIX &cr=null_matrix){
	int LIMIT=ar.size();
	if (t<0 or t>=LIMIT) return;
	if(s<0 or s>=LIMIT) return;
	VECTOR m;
	m=ar[t];
	ar[t]=ar[s];
	ar[s]=m;
	if(cr!=null_matrix){
		m=cr[t];cr[t]=cr[s];cr[s]=m;
	}
	if(br!=null_vector){
		swap(br[s],br[t]);
	}
}

/*---------------------------------------------------*/
bool partial_pivot(MATRIX &ar,VECTOR &br=null_vector,MATRIX &cr=null_matrix){
	int LIMIT=ar.size();
	for(int i=0;i<LIMIT;i++){
		if (ar[i][i]==0){
			bool isit=true;
			for(int t=0;t<LIMIT;t++){
				if(0!=ar[i][t])isit=false; 
			}
			if(isit) return (!isit);
			for(int j=i+1;i<LIMIT;i++){
				if(ar[j][i]!=0){
					for(int t=0;t<LIMIT;t++){
						swap(ar[j][t],ar[i][t]);
						if (cr!=null_matrix) swap(cr[j][t],cr[i][t]);
					}
					if (br!=null_vector)swap(br[j],br[i]);
					i=LIMIT;
				}
			}
		}
	}
	return true;
}

/*----------------------------------------------------------------*/

D_t dot(VECTOR &ar,VECTOR &br){
	int len=ar.size();
	if(len!=br.size()){cout<<"Dimensions are not equal.\n";return 0;}
	D_t sum=0;
	for(int i=0;i<len;i++) sum+=ar[i]*br[i];
	return sum;
}

/*-------------------------------------------------------------*/

complex<D_t> dot(vector<complex<D_t>> &ar,vector<complex<D_t>> &br){
	int len=ar.size();
	if(len!=br.size()){cout<<"Dimensions are not equal.\n";return 0;}
	complex<D_t> sum=0;
	for(int i=0;i<len;i++) sum+=ar[i]*br[i];
	return sum;
}

/*--------------------------------------------------------------*/

complex<D_t> inner_product(vector<complex<D_t>> &ar,vector<complex<D_t>> &br){
	vector<complex<D_t>> ddot;
	for(int i=0;i<ar.size();i++) ddot.push_back(conj(ar[i])) ;
	return dot(ddot,br);
}
/*---------------------------------------------------------------*/

MATRIX matrix_mult(MATRIX &a,MATRIX &b){
	MATRIX axb;
	for(int i=0;i<a.size();i++){
		VECTOR m;
		for(int j=0;j<b[0].size();j++){
			D_t sum=0;
			for(int k=0;k<b.size();k++)sum+=a[i][k]*b[k][j];
			m.push_back(sum);
		}
		axb.push_back(m);
	}
	return axb;
}

/*-------------------------------------------------------------*/

VECTOR matrix_mult(MATRIX &a,VECTOR &b){
	VECTOR axb;
	for(int i=0;i<a.size();i++) axb.push_back(dot(a[i],b));
	return axb;
}

/*-------------------------------------------------------------*/
/* VECTOR matrix_mult(double (*A)(int,VECTOR &),VECTOR &x){ */
/* 	VECTOR axb; */
/* 	for(int i=0;i<x.size();i++) axb.push_back(A(i,x)); */
/* 	return axb; */
/* } */
/*-------------------------------------------------------------*/

CVECTOR matrix_mult(CMATRIX &a,CVECTOR &b){
	CVECTOR axb;
	for(int i=0;i<a.size();i++) {
		C_t sum=0;
			for(int j=0;j<b.size();j++)sum+=a[i][j]*b[j];
		axb.push_back(sum);
	} // axb.push_back(dot(a[i],b));
	return axb;
}

/*-------------------------------------------------------------*/
CMATRIX matrix_mult(CMATRIX &a,CMATRIX &b){
	CMATRIX axb;
	for(int i=0;i<a.size();i++){
		CVECTOR m;
		for(int j=0;j<b[0].size();j++){
			C_t sum=0;
			for(int k=0;k<b.size();k++)sum+=a[i][k]*b[k][j];
			m.push_back(sum);
		}
		axb.push_back(m);
	}
	return axb;
}
