/* 
-> To find eigenvalues and eigenvectors
-> Jacobi and Power methods are here
 */

// INCLUDE "UTILITY.h" and "matrix_utilitly.cpp" before including this file


/*----------------------------JACOBI METHOD----------------------------------*/

// Find theta that correspends to make zero ij element
double theta(MATRIX &A,int i,int j){
    // tan(θ)=t satisfies t^2 - t*(A11-A22)/A12 -1 = 0
    double b = (A[i][i]-A[j][j])/A[i][j];
    return atan(b/2+sqrt(b*b+4)/2);
}
// Find the maximum off diagonal element to reduce it to zero
double find_max_offdiogonal(MATRIX &A,int &i,int &j){
    double max = 0;
    for(int k=0;k<A.size();k++){
        for(int l=k+1;l<A.size();l++){
            if(std::abs(A[k][l])>max){
                max = abs(A[k][l]);
                i = k;
                j = l;
            }
        }
    }
    return max;
}

// Jacobi method funciton
MATRIX Jacobi_eigenvalues(MATRIX &A,const double eps = 1e-9,const int max_iter = 1000){
    int N = A.size();
    int i=0,j=0;
    double a,b,c,s,theta_ij,temp0,temp1;
    double max = find_max_offdiogonal(A,i,j);
    int iter = 0;
    MATRIX ROTATION(N,VECTOR(N,0));  // The rotation matrix
    for(int i = 0;i<N;i++) ROTATION[i][i] = 1; // Since the initial rotation is identity
    while(max>eps && iter<max_iter){
        iter++;
        // iterating over columns
        double theta_ij = theta(A,i,j); // find θ
        c = cos(theta_ij);s = sin(theta_ij); // cos(θ) and sin(θ)
        temp0 = A[i][i];
        temp1 = A[j][j];
        // Update the matrix after θ rotation
        A[i][i] = A[i][i]*c*c+A[j][j]*s*s-2*A[i][j]*c*s;
        A[j][j] = temp0*s*s+A[j][j]*c*c+2*A[i][j]*c*s;
        A[i][j] = A[j][i] =0;
        for(int k=0;k<N;k++){
            if(k!=i && k!=j){
                temp0 = A[k][i];
                A[k][i] = c*temp0 - s*A[k][j];
                A[k][j] = s*temp0+ c*A[k][j];
                A[i][k] = A[k][i];
                A[j][k] = A[k][j];
            }
            //update rotation matrix by multiplying new rotation matrix
            temp0 = ROTATION[k][i];
            temp1 = ROTATION[k][j];
            ROTATION[k][i] = c*temp0-s*temp1;
            ROTATION[k][j] = s*temp0+c*temp1;
        }
        // find new maximum off-diagonal element
        max = find_max_offdiogonal(A,i,j);
    }
    cout<<"Iterations: "<<iter<<", Max = "<<max<<endl;
    return ROTATION;
}


/*----------------------------POWER METHOD-----------------------------------*/

VECTOR power_method(MATRIX &A,int no_eigens=1){
    int N = A.size();
    VECTOR x0(N),x1(N);
    double l1,l2,err ;
    for(int no=0;no<no_eigens;no++){
        err = 1; 
        if (no>0){
            for(int i = 0;i<N;i++){
                for(int j=0;j<N;j++) A[i][j]-= l1*x0[i]*x0[j];
            }
        }
        for(int i=0;i<N;i++) x0[i] = 1;//random_double(0,1);
        x1 = matrix_mult(A,x0);
        l1 = dot(x0,x1)/dot(x0,x0);
        double norm_x1 = sqrt(dot(x1,x1));
        for(int i =0;i<N;i++) x0[i] = x1[i]/norm_x1;
        int iter = 0;
        while (err>1e-11){
            iter++;
            x1 = matrix_mult(A,x0);
            l2 = dot(x0,x1);
            err = std::abs(l1-l2);
            double norm_x1 = sqrt(dot(x1,x1));
            for(int i =0;i<N;i++) x0[i] = x1[i]/norm_x1;
            l1 = l2;
        }
        cout<<"Largest Eigenvalue = "<<l1<<endl;
		display_vector(x0);
    }
    return x0;
}
