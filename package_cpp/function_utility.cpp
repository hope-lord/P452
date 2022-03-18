#include<iostream>
#include<vector>
#include<complex>
using namespace std;

static VECTOR null_vector1;
static vector<complex<D_t>> null_vector2;

/*--------------------------------------------------*/

D_t polynomial(D_t x,VECTOR &list){
    int n=list.size();
    D_t sum = 0;
    for(int i=0;i<n;i++) sum += list[i] * pow(x ,(n - 1 - i));
    return sum;
}

/*-------------------------------------------*/

complex<D_t> polynomial(complex<D_t> x,vector<complex<D_t>> &list){
    int n=list.size();
    complex<D_t> sum = 0;
    for(int i=0;i<n;i++) sum += list[i] * pow(x ,(n - 1 - i));
    return sum;
}

/*-------------------------------------------------------------------------------*/

D_t derivative(D_t (*f)(D_t),D_t x,D_t h=0.00001){
    return (f(x+h)-f(x-h))/(2*h);
}

/*----------------------------------------------------------------------------*/

D_t derivative(D_t (*f)(D_t,VECTOR&),D_t x,VECTOR &list,D_t h=0.00001){
    return (f(x+h,list)-f(x-h,list))/(2*h);
}

/*-------------------------------------------------------------------------------*/
complex<D_t> derivative(complex<D_t> (*f)(complex<D_t>),complex<D_t> x,D_t h=0.00001){
    return (f(x+h)-f(x-h))/(2*h);
}

/*-------------------------------------------------------------------------------*/

complex<D_t> derivative(complex<D_t> (*f)(complex<D_t>,vector<complex<D_t>>&),complex<D_t> x,vector<complex<D_t>> &list,D_t h=0.00001){
    return (f(x+h,list)-f(x-h,list))/(2*h);
}

/*----------------------------------------------------------------------------*/

D_t d_derivative(D_t (*f)(D_t),D_t x,D_t h=0.00001){
    return (f(x+h)+f(x-h)-2*f(x))/(h*h);
}

/*----------------------------------------------------------------------------*/

D_t d_derivative(D_t (*f)(D_t,VECTOR&),D_t x,VECTOR &list,D_t h=0.00001){
    return (f(x+h,list)+f(x-h,list)-2*f(x,list))/(h*h);
}

/*-------------------------------------------------------------------------------*/

complex<D_t> d_derivative(complex<D_t> (*f)(complex<D_t>),complex<D_t> x,D_t h=0.00001){
    return (f(x+h)+f(x-h)-f(x)-f(x))/(h*h);
}

/*----------------------------------------------------------------------------*/

complex<D_t> d_derivative(complex<D_t> (*f)(complex<D_t>,vector<complex<D_t>>&),complex<D_t> x,vector<complex<D_t>> &list,D_t h=0.00001){
    return (f(x+h,list)+f(x-h,list)-f(x,list)-f(x,list))/(h*h);
}

/*-------------------------------------------------------------------------------*/


