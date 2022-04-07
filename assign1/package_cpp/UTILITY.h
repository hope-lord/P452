#include<iostream>
#include<vector>
#include<cmath>
#include<complex>
#include<cstring>
#include<random>
#include<time.h>
using namespace std;

typedef double D_t;
typedef complex<D_t> C_t;


const D_t pi = M_PI;
const C_t I = C_t(0, 1);



typedef vector<D_t> VECTOR;
typedef vector<VECTOR> MATRIX; 
typedef vector<size_t> INDEX;
typedef vector<C_t> CVECTOR;
typedef vector<CVECTOR> CMATRIX;


double random_double(double min,double max){
    return min + (max-min)*rand()/RAND_MAX;
}


static VECTOR null_vector;
static MATRIX null_matrix;

typedef struct {
    D_t(*function)(D_t,void * );
    void * paras;
}general_function;




