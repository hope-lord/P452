

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

/*-----------------------------------------------------------------------*/
MATRIX rk2(VECTOR(*f)(double, VECTOR &,void*),VECTOR &y0,double x0,double xstop,void*args,long int N=1000){
    double h = (xstop - x0) / N;
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
        VECTOR d=f(x0,result,args);
        VECTOR k1_by_2=mul(h/2,d);
        VECTOR yn=add(result,k1_by_2);
        d=f(x0+h/2,yn,args);
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
MATRIX rk4(VECTOR(*f)(double,VECTOR&,void*),VECTOR& y0,double x0,double xstop,void*args,long int N=1000){
    double h=(xstop-x0)/N;
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
	int i = 0;
    while (i<N){
		i++;
        //finding k1 and k1/2
        VECTOR d = f(x0,result,args);
        VECTOR k1=mul(h,d);
        VECTOR k1_by_2 = mul(h / 2 , d);
        VECTOR k1_by_6=mul(h / 6, d);
        //finding k2 and k2/2
        VECTOR yn = add(result, k1_by_2);
        d = f(x0+h/2,yn,args);
        VECTOR k2 = mul(h, d);
        VECTOR k2_by_2=mul(h/2,d);
        VECTOR k2_by_3=mul(h/3,d);
        //finding k3
        yn = add(result, k2_by_2);
        d = f( x0+h/2,yn, args);
        VECTOR k3 = mul(h, d);
        VECTOR k3_by_3= mul(h/3, d);
        //finding k4
        yn = add(result, k3);
        d = f (x0+h,yn, args);
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
MATRIX Dirichlet_Boundary(VECTOR (*f)(double ,VECTOR&,void*),double ya,double yb,double a,double b,double choice1,double choice2,void*args,const double eps = 1e-6,int N = 1000){
  VECTOR init(2) ;
  MATRIX out;
  double c1,c2;
  // Find the second boundary at choice1
  init[0]=ya;init[1]=choice1;
  out = rk4(f,init,a,b,args,N);
  c1 = out[1][N];
  if (abs(c1-yb)<eps) return out; // return if the guess is good enough
  // Find the first boundary at choice2
  init[1] = choice2;
  out = rk4(f,init,a,b,args,N);
  c2 = out[1][N];
  if (abs(c2-yb)<eps) return out; // return if the guess is good enough
  // if both the boundry value for both choices are in one side then return
  if ((c1-yb)*(c2-yb)>0 ) {
	  printf("The choices are not valid.\n");
	  return out;
  }
  // interpolation to find the correct choice
  double choice = (yb-c1)*(choice2-choice1)/(c2-c1)+choice1;
  // Use the choice to find the boundry
  init[1] = choice;
  out = rk4(f,init,a,b,args,N);
  double c = out[1][N];
  if (abs(c-yb)<eps) return out; // return if the guess is good enough
  // Then use recursion to replace one of choice and find the correct choice
  if ((c1-yb)*(c-yb)>0) return Dirichlet_Boundary(f,ya,yb,a,b,choice,choice2,args,eps,N);
  else return Dirichlet_Boundary(f,ya,yb,a,b,choice1,choice,args,eps,N);
}
/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/
