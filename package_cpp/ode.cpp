


MATRIX rk2(general_deff_function*f,VECTOR y0,D_t x0,D_t xstop,size_t N=1000){
    D_t h = (xstop - x0) / N;
    VECTOR result;
    result.push_back(y0);
    VECTOR x;
    x.push_back(x0);
    while (std::abs(x0-xstop)>=std::abs(h)){
        int n=result.size();
        D_t d=f(result[-1],x0);
        D_t k1_by_2=h/2*d;
        D_t yn=result[n-1]+k1_by_2;
        d=f(yn,x0+h/2);
        D_t k2=h*d;
        yn=result[n-1]+k2;
        result.push_back(yn);
        x0+=h;
        x.push_back(x0);
    }
    MATRIX c{x,result};
    return c;
}
/*-----------------------------------------------------------------------*/
MATRIX rk2(D_t(*f)(D_t,D_t,VECTOR&),D_t y0,D_t x0,D_t xstop,VECTOR &args,long int N=1000){
    D_t h = (xstop - x0) / N;
    VECTOR result;
    result.push_back(y0);
    VECTOR x;
    x.push_back(x0);
    while (std::abs(x0-xstop)>=std::abs(h)){
        int n=result.size();
        D_t d=f(result[-1],x0,args);
        D_t k1_by_2=h/2*d;
        D_t yn=result[n-1]+k1_by_2;
        d=f(yn,x0+h/2,args);
        D_t k2=h*d;
        yn=result[n-1]+k2;
        result.push_back(yn);
        x0+=h;
        x.push_back(x0);
    }
    MATRIX c{x,result};
    return c;
}
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

void one_to_many(VECTOR& y0,MATRIX& y){
    int n=y0.size();
    for(int i=0;i<n;i++){
        y[i].push_back(y0[i]);
    }
}

VECTOR mul(D_t h,VECTOR &a){
    VECTOR b;
    for(int i=0;i<a.size();i++)
        b.push_back(a[i]*h);
    return b;
}

VECTOR add(VECTOR &a,VECTOR &b){
    VECTOR c;
    for(int i=0;i<a.size();i++)
        c.push_back(a[i]+b[i]);
    return c;
}


MATRIX rk2(VECTOR(*f)(VECTOR &,D_t),VECTOR &y0,D_t x0,D_t xstop,long int N=1000){
    D_t h = (xstop - x0) / N;
    int n=y0.size();
    MATRIX final_result;
    VECTOR e;
    for(int i=0;i<n;i++)
        final_result.push_back(e);
    one_to_many(y0,final_result);
    VECTOR result;
    for(int i=0;i<n;i++)
        result.push_back(y0[i]);
    VECTOR x;
    x.push_back(x0);
    while (std::abs(x0-xstop)>=std::abs(h)){
        VECTOR d=f(result,x0);
        VECTOR k1_by_2=mul(h/2,d);
        VECTOR yn=add(result,k1_by_2);
        d=f(yn,x0+h/2);
        VECTOR k2=mul(h,d);
        yn=add(result,k2);
        result=yn;
        one_to_many(yn,final_result);
        x0+=h;
        x.push_back(x0);
    }
    final_result.insert(final_result.begin(),x);
    return final_result;
}

/*-----------------------------------------------------------------------*/
MATRIX rk2(VECTOR(*f)(VECTOR &,D_t,VECTOR &),VECTOR &y0,D_t x0,D_t xstop,VECTOR &args,long int N=1000){
    D_t h = (xstop - x0) / N;
    int n=y0.size();
    MATRIX final_result;
    VECTOR e;
    for(int i=0;i<n;i++)
        final_result.push_back(e);
    one_to_many(y0,final_result);
    VECTOR result;
    for(int i=0;i<n;i++)
        result.push_back(y0[i]);
    VECTOR x;
    x.push_back(x0);
    while (std::abs(x0-xstop)>=std::abs(h)){
        VECTOR d=f(result,x0,args);
        VECTOR k1_by_2=mul(h/2,d);
        VECTOR yn=add(result,k1_by_2);
        d=f(yn,x0+h/2,args);
        VECTOR k2=mul(h,d);
        yn=add(result,k2);
        result=yn;
        one_to_many(yn,final_result);
        x0+=h;
        x.push_back(x0);
    }
    final_result.insert(final_result.begin(),x);
    return final_result;
}

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

MATRIX rk4(D_t(*f)(D_t,D_t),D_t y0,D_t x0,D_t xstop,long int N=1000){
    D_t h=(xstop-x0)/N;
    VECTOR result{y0};
    VECTOR x{x0};
    while (std::abs(x0 - xstop) >= std::abs(h)){
        int n=result.size();
        //finding k1 and k1/2
        D_t d = f(result[n-1], x0);
        D_t k1=h* d;
        D_t k1_by_2 = h / 2 * d;
        D_t k1_by_6=h / 6* d;
        //finding k2 and k2/2
        D_t yn = result[n-1]+ k1_by_2;
        d = f(yn, x0 + h / 2);
        D_t k2 = h* d;
        D_t k2_by_2=h/2*d;
        D_t k2_by_3=h/3*d;
        //finding k3
        yn = result[n-1]+ k2_by_2;
        d = f( yn, x0 + h / 2);
        D_t k3 = h* d;
        D_t k3_by_3= h/3* d;
        //finding k4
        yn = result[n-1]+ k3;
        d = f (yn, x0+h);
        D_t k4 = h* d;
        D_t k4_by_6=h/6*d;
        //finding the correct yn
        yn=result[n-1]+k1_by_6+k2_by_3+k3_by_3+k4_by_6;
        result.push_back(yn);
        x0 += h;
        x.push_back(x0);
    }
    MATRIX c{x,result};
    return c;
}
/*-----------------------------------------------------------------------*/
MATRIX rk4(D_t(*f)(D_t,D_t,VECTOR&),D_t y0,D_t x0,D_t xstop,VECTOR& args,long int N=1000){
    D_t h=(xstop-x0)/N;
    VECTOR result{y0};
    VECTOR x{x0};
    while (std::abs(x0 - xstop) >= std::abs(h)){
        int n=result.size();
        //finding k1 and k1/2
        D_t d = f(result[n-1], x0,args);
        D_t k1=h* d;
        D_t k1_by_2 = h / 2 * d;
        D_t k1_by_6=h / 6* d;
        //finding k2 and k2/2
        D_t yn = result[n-1]+ k1_by_2;
        d = f(yn, x0 + h / 2,args);
        D_t k2 = h* d;
        D_t k2_by_2=h/2*d;
        D_t k2_by_3=h/3*d;
        //finding k3
        yn = result[n-1]+ k2_by_2;
        d = f( yn, x0 + h / 2,args);
        D_t k3 = h* d;
        D_t k3_by_3= h/3* d;
        //finding k4
        yn = result[n-1]+ k3;
        d = f (yn, x0+h,args);
        D_t k4 = h* d;
        D_t k4_by_6=h/6*d;
        //finding the correct yn
        yn=result[n-1]+k1_by_6+k2_by_3+k3_by_3+k4_by_6;
        result.push_back(yn);
        x0 += h;
        x.push_back(x0);
    }
    MATRIX c{x,result};
    return c;
}

/*-----------------------------------------------------------------------*/

MATRIX rk4(VECTOR(*f)(VECTOR&,D_t),VECTOR& y0,D_t x0,D_t xstop,long int N=1000){
    D_t h=(xstop-x0)/N;
    int n=y0.size();
    MATRIX final_result;
    VECTOR e;
    for(int i=0;i<n;i++)
        final_result.push_back(e);
    one_to_many(y0,final_result);
    VECTOR result;
    for(int i=0;i<n;i++)
        result.push_back(y0[i]);
    VECTOR x;
    x.push_back(x0);
    while (std::abs(x0 - xstop) >= std::abs(h)){
        //finding k1 and k1/2
        VECTOR d = f(result, x0);
        VECTOR k1=mul(h,d);
        VECTOR k1_by_2 = mul(h / 2 , d);
        VECTOR k1_by_6=mul(h / 6, d);
        //finding k2 and k2/2
        VECTOR yn = add(result, k1_by_2);
        d = f(yn, x0 + h / 2);
        VECTOR k2 = mul(h, d);
        VECTOR k2_by_2=mul(h/2,d);
        VECTOR k2_by_3=mul(h/3,d);
        //finding k3
        yn = add(result, k2_by_2);
        d = f( yn, x0 + h / 2);
        VECTOR k3 = mul(h, d);
        VECTOR k3_by_3= mul(h/3, d);
        //finding k4
        yn = add(result, k3);
        d = f (yn, x0+h);
        VECTOR k4 = mul(h, d);
        VECTOR k4_by_6=mul(h/6,d);
        //finding the correct yn
        yn=add(result,k1_by_6);
        yn=add(yn,k2_by_3);
        yn=add(yn,k3_by_3);
        yn=add(yn,k4_by_6);
        result=yn;
        one_to_many(yn,final_result);
        x0 += h;
        x.push_back(x0);
    }
    final_result.insert(final_result.begin(),x);
    return final_result;
}



MATRIX rk4(VECTOR(*f)(VECTOR&,D_t,VECTOR&),VECTOR& y0,D_t x0,D_t xstop,VECTOR& args,long int N=1000){
    D_t h=(xstop-x0)/N;
    int n=y0.size();
    MATRIX final_result;
    VECTOR e;
    for(int i=0;i<n;i++)
        final_result.push_back(e);
    one_to_many(y0,final_result);
    VECTOR result;
    for(int i=0;i<n;i++)
        result.push_back(y0[i]);
    VECTOR x;
    x.push_back(x0);
    while (std::abs(x0 - xstop) >= std::abs(h)){
        //finding k1 and k1/2
        VECTOR d = f(result, x0,args);
        VECTOR k1=mul(h,d);
        VECTOR k1_by_2 = mul(h / 2 , d);
        VECTOR k1_by_6=mul(h / 6, d);
        //finding k2 and k2/2
        VECTOR yn = add(result, k1_by_2);
        d = f(yn, x0 + h / 2,args);
        VECTOR k2 = mul(h, d);
        VECTOR k2_by_2=mul(h/2,d);
        VECTOR k2_by_3=mul(h/3,d);
        //finding k3
        yn = add(result, k2_by_2);
        d = f( yn, x0 + h / 2,args);
        VECTOR k3 = mul(h, d);
        VECTOR k3_by_3= mul(h/3, d);
        //finding k4
        yn = add(result, k3);
        d = f (yn, x0+h,args);
        VECTOR k4 = mul(h, d);
        VECTOR k4_by_6=mul(h/6,d);
        //finding the correct yn
        yn=add(result,k1_by_6);
        yn=add(yn,k2_by_3);
        yn=add(yn,k3_by_3);
        yn=add(yn,k4_by_6);
        result=yn;
        one_to_many(yn,final_result);
        x0 += h;
        x.push_back(x0);
    }
    final_result.insert(final_result.begin(),x);
    return final_result;
}


/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/