#include"UTILITY.h"


D_t midpoint(general_function *f,D_t a,D_t b,D_t f2p_max,D_t prec=0.001){
    if (a==b) return 0;
    float c=1;
    if (a>b){
        D_t t=a;a=b;b=t;
        c=-1;
    }
    int N=pow((pow(b-a,3) *std::abs(f2p_max)/24/prec),0.5);
    N=ceil(N);
    if (N==0) N=1;  //if linear function is given
    D_t h=(b-a)/N;
    D_t sum=0;
    for(int i=0;i<N;i++){
        D_t mid_x=a+i*h+h/2;
        sum+=f->function(mid_x,f->paras); //as weight=1 for each x_i;
    }
    return c*h*sum;
}
/*------------------------------------------------------------------------------------------------*/

D_t trapezoidal(general_function*f,D_t a,D_t b,D_t f2p_max,D_t prec=0.001){
    if (a==b) return 0;
    float c=1;
    if (a>b){
        D_t t=a;a=b;b=t;
        c=-1;
    }
    int N=pow((pow(b-a,3) *std::abs(f2p_max)/12/prec),0.5);
    N=ceil(N);
    if (N==0) N=1;  //if linear function is given
    D_t h=(b-a)/N;
    D_t sum=0;
    for(int i=0;i<=N;i++){
        int w=2;
        if ((i==0) or (i==N)) w=1; //weight is 1 for first and last
        D_t xi=a+i*h;
        sum+=w*f->function(xi,f->paras);
    }
    return c*h*sum/2;
}
/*-------------------------------------------------------------------------*/

D_t simpson(general_function*f,D_t a,D_t b,D_t f4p_max,D_t prec=0.001){
    if (a==b) return 0;
    float c=1;
    if (a>b){
        D_t t=a;a=b;b=t;
        c=-1;
    }
    int N=pow((pow(b-a,5) *std::abs(f4p_max)/180/prec),0.25);
    N=ceil(N);
    if (N==0) N=1;  //if linear function is given
    D_t h=(b-a)/N;
    D_t sum=0;
    for(int i=0;i<=N;i++){
        int w;
        if (i==0 or i==N) w=1; //weight is 1 for first and last
        else if (i%2==1) w=4; //for odd weight=4
        else w=2; //for even weight=2
        D_t xi=a+i*h;
        sum+=w*f->function(xi,f->paras);
    }
    return c*h*sum/3;
}
/*--------------------------------------------------------------------------------------------------------------------------------------------*/

/*VECTOR monte_carlo(D_t(*f)(D_t),D_t a,D_t b,size_t N=1e6) {
    D_t sum=0;
    D_t sum2=0;
    D_t inv_pdf=b-a; //inverse pdf
    if (inv_pdf==0){ 
        VECTOR c{0,0};    
        return c;
    }
    std::default_random_engine generator;
    std::uniform_real_distribution<D_t> distribution(a,b);
    for(int i=0;i<N;i++){
        D_t value=f(distribution(generator));
        sum+=value;
        sum2 += pow(value,2);
    }
    VECTOR c{inv_pdf*sum/N, pow((sum2/N/inv_pdf - pow((sum/N/inv_pdf),2)),0.5)};    
    return c;
}

/*-------------------------------------------------------------------------*/
/*VECTOR monte_carlo(D_t(*f)(D_t,VECTOR&),D_t a,D_t b,VECTOR &args,long long int N=1e6) {
    D_t sum=0;
    D_t sum2=0;
    D_t inv_pdf=b-a; //inverse pdf
    if (inv_pdf==0){ 
        VECTOR c{0,0};    
        return c;
    }
    std::default_random_engine generator;
    std::uniform_real_distribution<D_t> distribution(a,b);
    for(int i=0;i<N;i++){
        D_t value=f(distribution(generator),args);
        sum+=value;
        sum2 += pow(value,2);
    }
    VECTOR c{inv_pdf*sum/N, pow((sum2/N/inv_pdf - pow((sum/N/inv_pdf),2)),0.5)};    
    return c;
}
*/
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/

VECTOR monte_carlo(general_monte_function*f,VECTOR &a,VECTOR &b,size_t N=1e6) {
    D_t sum=0;
    D_t sum2=0;
    D_t inv_pdf=1; //inverse pdf
    for(size_t i=0;i<a.size();i++) 
        inv_pdf*=(b[i]-a[i]);//for multivariate function inverse pdf is multiplication of range differences
    
    if (inv_pdf==0){ 
        VECTOR c{0,0};    
        return c;
    }
    std::default_random_engine generator;
    vector<std::uniform_real_distribution<D_t>> distribution;
    
    for(size_t i=0;i<a.size();i++){
        std::uniform_real_distribution<D_t> c(a[i],b[i]);
        distribution.push_back(c);
    }

    for(size_t i=0;i<N;i++){
        VECTOR x;
        for(size_t i=0;i<a.size();i++)
            x.push_back(distribution[i](generator));
        
        D_t value=f->function(x,f->paras);
        sum+=value;
        sum2 += pow(value,2);
    }
    VECTOR c{inv_pdf*sum/N, pow((sum2/N/inv_pdf - pow((sum/N/inv_pdf),2)),0.5)};    
    return c;
}

/*-------------------------------------------------------------------------*/
/*VECTOR monte_carlo(D_t(*f)(VECTOR&,VECTOR&),VECTOR &a,VECTOR &b,VECTOR &args,long long int N=1e6) {
    D_t sum=0;
    D_t sum2=0;
    D_t inv_pdf=1; //inverse pdf
    for(int i=0;i<a.size();i++) 
        inv_pdf*=(b[i]-a[i]);//for multivariate function inverse pdf is multiplication of range differences
    
    if (inv_pdf==0){ 
        VECTOR c{0,0};    
        return c;
    }
    std::default_random_engine generator;
    vector<std::uniform_real_distribution<D_t>> distribution;
    
    for(int i=0;i<a.size();i++){
        std::uniform_real_distribution<D_t> c(a[i],b[i]);
        distribution.push_back(c);
    }

    for(int i=0;i<N;i++){
        VECTOR x;
        for(int i=0;i<a.size();i++)
            x.push_back(distribution[i](generator));
        
        D_t value=f(x,args);
        sum+=value;
        sum2 += pow(value,2);
    }
    VECTOR c{inv_pdf*sum/N, pow((sum2/N/inv_pdf - pow((sum/N/inv_pdf),2)),0.5)};    
    return c;
}
*/
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/



/*-------------------------------------------------------------------------*/