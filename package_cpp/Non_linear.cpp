

VECTOR make_bracket(D_t(*f)(D_t),D_t a,D_t b,D_t factor){
    VECTOR ar;
    int i=0;
    if (a>b){D_t t=a;a=b;b=t;}
    while((i<15) and (f(a)*f(b)>0)){
        i++;
        if (std::abs(f(a)) < std::abs(f(b))) a -= (b - a) * factor;
        else b += (b - a) * factor;
    }
    if (f(a)*f(b)>0){
        cout<<"Give proper a & b.\n";
        return ar;
    }
    else if(f(a)*f(b)==0){
        if (f(a)==0) {
            ar[0]=a;ar[1]=a;
            return ar;
        }
        else{
            ar[0]=b;ar[1]=b;
            return ar;
        }
    }
    else{
        ar[0]=a;ar[1]=b;
        return ar;
    }
}
/*-----------------------------------------*/

VECTOR make_bracket(D_t(*f)(D_t,VECTOR&),D_t a,D_t b,D_t factor,VECTOR &args){
    VECTOR ar;
    int i=0;
    if (a>b){D_t t=a;a=b;b=t;}
    while((i<15) and (f(a,args)*f(b,args)>0)){
        i++;
        if (std::abs(f(a,args)) < std::abs(f(b,args))) a -= (b - a) * factor;
        else b += (b - a) * factor;
    }
    if (f(a,args)*f(b,args)>0){
        cout<<"Give proper a & b.\n";
        return ar;
    }
    else if(f(a,args)*f(b,args)==0){
        if (f(a,args)==0) {
            ar[0]=a;ar[1]=a;
            return ar;
        }
        else{
            ar[0]=b;ar[1]=b;
            return ar;
        }
    }
    else{
        ar[0]=a;ar[1]=b;
        return ar;
    }
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
D_t bisection(D_t (*func)(D_t),D_t a,D_t b,D_t factor,D_t epsi){
    VECTOR ar=make_bracket(func,a,b,factor);
    //cout<<*ar<<", "<<*(ar+1)<<endl;
    if (ar==null_vector1) return 0;  //bracketing was not possible
    else if(ar[0]==ar[1]) return ar[1];
    else{
        a=ar[0];
        b=ar[1];
        D_t c=(a+b)/2;int i=0;
        while(i<200 and std::abs(b-a)>epsi and func(c)!=0){
            i++;
            if (func(c)*func(a)>0) a=c;
            else b=c;
            c=(a+b)/2;
        }
        if (i==50) cout<<"The limit has reached.\n";
        return c;
    }
}

/*--------------------------------------------------------------------------*/

D_t bisection(D_t (*func)(D_t,VECTOR&),D_t a,D_t b,D_t factor,D_t epsi,VECTOR& args){
    VECTOR ar=make_bracket(func,a,b,factor,args);
    //cout<<*ar<<", "<<*(ar+1)<<endl;
    if (ar==null_vector1) return 0;  //bracketing was not possible
    else if(ar[0]==ar[1]) return ar[1];
    else{
        a=ar[0];
        b=ar[1];
        D_t c=(a+b)/2;int i=0;
        while(i<200 and std::abs(b-a)>epsi and func(c,args)!=0){
            i++;
            if (func(c,args)*func(a,args)>0) a=c;
            else b=c;
            c=(a+b)/2;
        }
        if (i==50) cout<<"The limit has reached.\n";
        return c;
    }
}

/*--------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
D_t regula_falsi(D_t (*func)(D_t),D_t a,D_t b,D_t factor,D_t epsi){
    VECTOR ar=make_bracket(func,a,b,factor);
    if (ar==null_vector1) return 0; //bracketing was not possible
    else if (ar[0]==ar[1]) return ar[1]; //a or b hits the root
    else{
        D_t c= (func(b)*a-b*func(a))/(func(b)-func(a));int i=0;
        //min(abs(c-a),abs(b-a))<epsi to stop loop
        while (i < 200 and std::abs(c-a)>epsi and std::abs(c-b)>epsi and func(c) != 0){
            i++;
            if (func(c)*func(a)>0) a=c;
            else b=c;
            c = (func(b) * a - b * func(a)) / (func(b) - func(a));
        }
        if (i==200) cout<<"The limit has reached.\n";
        return c;
    }
}
/*----------------------------------------------------------------------------*/
D_t regula_falsi(D_t (*func)(D_t,VECTOR&),D_t a,D_t b,D_t factor,D_t epsi,VECTOR &args){
    VECTOR ar=make_bracket(func,a,b,factor,args);
    if (ar==null_vector1) return 0; //bracketing was not possible
    else if (ar[0]==ar[1]) return ar[1]; //a or b hits the root
    else{
        D_t c= (func(b,args)*a-b*func(a,args))/(func(b,args)-func(a,args));int i=0;
        //min(abs(c-a),abs(b-a))<epsi to stop loop
        while (i < 200 and std::abs(c-a)>epsi and std::abs(c-b)>epsi and func(c,args) != 0){
            i++;
            if (func(c,args)*func(a,args)>0) a=c;
            else b=c;
            c = (func(b,args) * a - b * func(a,args)) / (func(b,args) - func(a,args));
        }
        if (i==200) cout<<"The limit has reached.\n";
        return c;
    }
}

/*-----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

D_t newton_raphson(D_t (*f)(D_t),D_t x0,D_t epsi){
    if (f(x0)==0) return x0;
    D_t fp=derivative(f,x0);
    if (fp==0){
        cout<<"Re-enter proper x_0.\n";
        return 0;
    }
    D_t x1=x0-f(x0)/fp;
    int i=0;
    while (std::abs(x1-x0)>epsi and f(x1)!=0 and i<200){
        i++;
        x0=x1;
        fp=derivative(f, x0);
        if (fp == 0){
            cout<<"Re-enter proper x_0.\n";
            return 0;
        }
        x1 = x0 - f(x0) / fp;
    }
    return x1;
}
/*---------------------------------------------------------------*/

D_t newton_raphson(D_t (*f)(D_t,VECTOR&),D_t x0,D_t epsi,VECTOR &args){
    if (f(x0,args)==0) return x0;
    D_t fp=derivative(f,x0,args);
    if (fp==0){
        cout<<"Re-enter proper x_0.\n";
        return 0;
    }
    D_t x1=x0-f(x0,args)/fp;
    int i=0;
    while (std::abs(x1-x0)>epsi and f(x1,args)!=0 and i<200){
        i++;
        x0=x1;
        fp=derivative(f, x0,args);
        if (fp == 0){
            cout<<"Re-enter proper x_0.\n";
            return 0;
        }
        x1 = x0 - f(x0,args) / fp;
    }
    return x1;
}

/*-----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/

vector<complex<D_t>> laguerre(VECTOR &the_list,D_t x0,D_t epsi){
    vector<complex<D_t>> list;
    complex<D_t> a0=x0;
    int t=the_list.size();
    for(int i=0;i<t;i++) list.push_back(the_list[i]);
    vector<complex<D_t>> root;
    while (root.size()<t){
        int nn=list.size();
        complex<D_t> n=nn;
        if (std::abs(polynomial(a0,list))==0) root.push_back(a0);
        else{
            if (list.size()==2){
                    root.push_back(-list[1]/list[0]);
                    return root;
            }
            complex<D_t> G=derivative(polynomial,a0,list)/polynomial(a0,list);
            complex<D_t> H=pow(G,2)-d_derivative(polynomial,a0,list)/polynomial(a0,list);
            complex<D_t> a,ll=pow(n*(n*H-pow(G,2))-(n*H-pow(G,2)),0.5);

            if (std::abs(G+ll)>std::abs(G-ll)) a=n/(G+ll);
            else if (std::abs(G+ll)==0){
                cout<<"Re- enter proper x_0.\n";
                return root;
            }
            else a=n/(G-ll);
            a0-=a;
            int j=0;
            while (std::abs(a)>epsi and std::abs(polynomial(a0,list))!=0 and j<50){
                j++;
                G = derivative(polynomial,a0,list) / polynomial(a0,list);
                H = pow(G,2) - d_derivative(polynomial,a0,list) / polynomial(a0,list);
                ll=pow(n*(n*H-pow(G,2))-(n*H-pow(G,2)),0.5);
                if (std::abs(G+ll)>std::abs(G-ll)) a=n/(G+ll);
                else if (std::abs(G+ll)==0){
                cout<<"Re- enter proper x_0.\n";
                return root;
            }
            else a=n/(G-ll);
                a0-=a;
            }
            if (j==50) cout<<"The limit has reached.\n";
            root.push_back(a0);
        }
        if (root.size() != 0){
            vector<complex<D_t>> new_list;
            new_list.push_back(list[0]);
            complex<D_t> last_root = root.back();
            for (int i=1;i<list.size()-1;i++)
                new_list.push_back(last_root * new_list.back() + list[i]);
            list = new_list;
        }
    }
    return root;

   
}