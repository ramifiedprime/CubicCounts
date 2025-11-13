#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long gcd(long a, long b){//TODO
    long r;
    while(b!=0){
        r=a%b;
        a=b;
        b=r;
    }
    return labs(a);
}

_Bool member(long x, long* list){
    long n = list[0];
    long lo = 1, hi = n;
    while (lo <= hi) {
        long mid = lo + (hi - lo) / 2;
        if (list[mid] == x) return 1;
        if (list[mid] < x) lo = mid + 1;
        else hi = mid - 1;
    }
    return 0;
}
////////////////////////////////////////////////////////////////////////////////
// Lazy implementation of Erastosthenes sieve
// Input:
//      - int X: integer bound
// Output:
//      - int* primes: array with primes[0] being the number of primes up to X
//                     and primes[i] being the ith prime for i>=1
long* primes_up_to(long X, long sqrtX){
    // printf("Initialising p...\n");
    _Bool* comp = (_Bool*)calloc(X+1, sizeof(_Bool));//comp[n]=true iff n composite
    long i,j;
    long pi=0;
    for(i=2; i<= X; i++){
        if(!comp[i]){
            pi++;
            if(i<=sqrtX){
                for(j=2; j<= X/i; j++){
                    comp[i*j] = 1;
                }
            }
        }
    }
    // printf("\tcomp initialised\n")
    long* primes = (long*)calloc(pi+1,sizeof(long));//extra entry for storing pi
    primes[0]=pi;
    j=1;
    for(i=2; i<= X; i++){
        if(!comp[i]){
            primes[j] = i;
            j++;
        }
    }
    // printf("\tprimes initialised")
    free(comp);
    // printf("...done")
    return primes; //first entry is number of primes, rest are primes
}

////////////////////////////////////////////////////////////////////////////////
// implementation of finding index in sub-algorithm 5.1 (init)
// Input:
//      - int* q: array of primes
//      - int tM: M from algorithm
// Output:
//      - int i: index of correct prime in list
long init_index_find(long* p, long X, long M){
    long r[2] = {1,p[0]};
    long i = (r[0]+r[1]+1)/2;
    long tM = 3*M;
    double num=1;
    while(i!=r[1]){
        num = X/(p[i]*log(p[i]))*(1+(1.0/(2*log(p[i]))));
        if(num <= tM){
            r[1]=i;
        }else{
            r[0]=i;
        }
        i = (r[0]+r[1]+1)/2;
    }
    if(num > tM){
        printf("Error: Memory constraint too small for bound.\n");
        exit(1);
    }//At this stage we've completed the first pass
    if(p[i]<=53){
        double J=0;
        double JJ=(((double)tM)/X - (1/log(59))*(1+(1/(2*log(59)))));
        i = 16;
        while(J<=JJ && i>3){
            J+=1.0/(p[i]*p[i]);
            i--;
        }
    }
    return i;
}

//gets ordered list of all integers up to X which are divisible by pp[i] for some i>=index 
long* get_list(long* pp, long index, long X){
    long i,k,x,num=0;
    _Bool* blist = (_Bool*)calloc(X+2, sizeof(_Bool));
    for(i=index; i<=pp[0]; i++){
        for(x=pp[i]; x<X; x+=pp[i]){
            if(!blist[x]){num++; blist[x]=1;}
        }
    }
    // printf("\tsuccessfully passed blist construction\n");
    long* list = (long*)calloc(num+2, sizeof(long));
    i=1;
    for(x=6; x<=X+1; x+=6){
        for(k=-1;k<=1 && x+k<=X;k+=2){
            if(blist[x+k]){list[i]=x+k;i++;}
        }
    }
    list[0]=i-1;
    // printf("\tsuccessfully passed list construction\n");
    free(blist);
    return list;
}

long* get_sqfull(long* pp, long X){
    long i,n;
    long* sqfull = (long*)calloc((X+1),sizeof(long));
    sqfull[0]=X;
    sqfull[1]=0;
    for(i=3; i<=pp[0] && pp[i]<=X; i++){
        for(n=pp[i]; n<=X; n+=pp[i]){
            sqfull[n]=1;
        }
    }
    return sqfull;
}

// u a divisor of 72 divisible by 3, t coprime to 72, K=Q(sqrt(-ut/3))
long get_quad_disc(long t, long u){
    t=t*(u/3);
    if(u%4==0){
        t=t/4;
    }
    if(t%4!=3){
        return 4*t;
    }
    return t;
}



// Currently returns t or 0
long test(long a, long b, long c, long d, long D, long f, long* sqfull, long sqrtX, long* list, long* pp, long index){
    long a9,d9,t,i,u;
    double sqrt3X=floor(1.7320508075688772*sqrtX);
    if ((D%27)==0 && (f%3!=0)){return 0;} // Not in V3 and not (1^3)
    else if(f%3==0){// checking if in U3 given that it's 1^3
        a9 = a%9;
        d9 = d%9;
        if(a9==0 || d9==0){return 0;}
        else if(((a9%3)==0)&&((d9%3)==0)){return 0;}
        else if((((a9-d9)%3)==0) && ((a9-b+c-d9)%9)==0){return 0;}
        else if((((a9+d9)%3)==0) && ((a9+b+c+d9)%9)==0){return 0;}
    }
    if(f>sqfull[0]){printf("Overflow on sqfull: sqfull[0]=%ld, f=%ld\n", sqfull[0],f);}
    if(sqfull[f]){return 0;}
    t = labs(D/(f*f));
    u = gcd(t,72);
    t = t/u;
    if(gcd(t,f)!=1){return 0;}
    if (t<=sqrt3X && sqfull[t]){return 0;} 
    else if(member(t,list)){return 0;}
    else{
        for(i=2;i<=index-1;i++){
            if(t%pp[i]==0){return 0;}
        }
    }
    return get_quad_disc(t, u);
}

long is_complex_field(long a, long b, long c, long d, long P, long Q, long R, long D, long f,long* sqfull, long sqrtX, long* list, long* pp, long index){
        if((D <= 0) || (D%16==0) || (D%16==4 && (P%2!=0 || R%2!=0))){
            return 0;
        }
        // printf("\tpassed U2 check, now running test\n");
        return test(a,b,c,d,D,f,sqfull,sqrtX,list,pp,index);
}

double U(long a, long b){
    if(a>=2.0*b/3){return ((double)b*b)/(3.0*a);}
    else{return b-3*((double)a)/4.0;}
}

int CCFCCF(long X, FILE *fptr, _Bool verbose){
    long a,b,c,d,D,f,P,Q,R;
    long i=0, tX=3*X, B=(long)sqrt(X);
    if(verbose){printf("Initialising...\n");}
    long* p = primes_up_to(X,B);
    if(verbose){printf("...p initialised\n");}
    long* pp = (long*)malloc((p[0]+1)*sizeof(long));
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    if(verbose){printf("...pp initialised\n");}
    long index = init_index_find(p, X, X);
    if(verbose){printf("...index found\n");}
    long* sqfull = get_sqfull(pp,sqrt(3*X));
    if(verbose){printf("...sqfull initialised\n");}
    long* list = get_list(pp, index, X);
    if(verbose){printf("...list initialised\n");}
    if(verbose){printf("...done.\n");}
    // do looping
    double A_bd,B_bd,C_bd,D_ubd;
    long D_lbd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    long quadbd;
    // b=0 looping
    if(verbose){printf("Working through b=0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        C_bd=pow(((double)X)/(4*a), 1.0/3);
        for(c=1; c<=C_bd; c++){
            D_lbd=0;
            if(c<=a){D_lbd = (long)floor(sqrt(a*(a-c)));}
            for(d=D_lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                P= -3*a*c;
                Q= -9*a*d;
                R= c*c;
                D= Q*Q-4*P*R;
                if(D>tX || D<=0){continue;} // IF DISC
                f=gcd(P,gcd(Q,R));
                if(is_complex_field(a,0,c,d,P,Q,R,D,f,sqfull,B,list,pp,index)){
                    fprintf(fptr,"[%ld,%d,%ld,%ld]\n",a,0,c,d);
                }
            }
        }
    }
    if(verbose){printf("...done.\nWorking through b>0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
        for(b=1; b<=B_bd; b++){ 
            C_bd=U(a,b)+pow(((double)X)/(4.0*a), 1.0/3);
            for(c=(1-b); c<=C_bd;c++){
                D_lbd =(long)floor(-(((double)a-b)*(a-b+c))/a);
                D_ubd = (((double)a+b)*(a+b+c))/a;
                quadbd = a*(a-c);
                for(d=D_lbd+1; d<D_ubd; d++){
                    if(d*(d-b) < quadbd){continue;}//Improvement poss here, inc function
                    P= b*b-3*a*c;
                    Q= b*c-9*a*d;
                    R= c*c-3*b*d;
                    D= Q*Q-4*P*R;
                    if (D>tX || D<=0){continue;} //IF DISC
                    f=gcd(P,gcd(Q,R));
                    if(is_complex_field(a,b,c,d,P,Q,R,D,f,sqfull,B,list,pp,index)){
                        fprintf(fptr,"[%ld,%ld,%ld,%ld]\n",a,b,c,d);
                    }
                }
            }
        }
    }
    if(verbose){printf("...done.\n");}
    free(p);
    free(pp);
    free(sqfull);
    free(list);
    return 0;
}

int CCFCCFPRP(long B, FILE *fptr, _Bool verbose){
    long a,b,c,d,D,f,P,Q,R,check;
    long i=0, X=B*B, tX=3*X;
    if(verbose){printf("Initialising...\n");}
    long* p = primes_up_to(3*B,sqrt(3*B));
    if(verbose){printf("...p initialised\n");}
    long* pp = (long*)malloc((p[0]+1)*sizeof(long));
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    if(verbose){printf("...pp initialised\n");}
    long index = init_index_find(p, 3*B, 3*B);
    if(verbose){printf("...dumping p...");}
    free(p);
    if(verbose){printf("done.\n");}
    if(verbose){printf("...index found\n");}
    long* sqfull = get_sqfull(pp,3*B);
    if(verbose){printf("...sqfull initialised\n");}
    long* list = get_list(pp, index, 3*B);
    if(verbose){printf("...list initialised\n");}
    if(verbose){printf("...done.\n");}
    // do looping
    double A_bd,B_bd,C_bd,D_ubd;
    long D_lbd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    long quadbd;
    // b=0 looping
    if(verbose){printf("Working through b=0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        if(verbose){printf("...a=%ld\n", a);}
        C_bd=pow(((double)X)/(4*a), 1.0/3);
        for(c=1; c<=C_bd; c++){
            // if(verbose){printf("...a=%ld, c = %ld\n", a, c);}
            P= -3*a*c;
            R= c*c;
            D_lbd=0;
            if(c<=a){D_lbd = (long)floor(sqrt(a*(a-c)));}
            for(d=D_lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                Q= -9*a*d;
                D= Q*Q-4*P*R;
                if(D<=0){continue;}
                if(D>tX){break;} // check disc, note as b=0 then  once D gets big it only gets bigger as d increases for fixed a,c
                f=gcd(P,gcd(Q,R));
                if(D>3*B*f){continue;} // IF PRP
                check=is_complex_field(a,0,c,d,P,Q,R,D,f,sqfull,B,list,pp,index);
                if(check){
                    fprintf(fptr,"%ld,%d,%ld,%ld,%ld\n",a,0,c,d,check);
                }
            }
        }
    }
    if(verbose){printf("...done.\nWorking through b>0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
        for(b=1; b<=B_bd; b++){ 
            if(verbose){printf("...a=%ld, b=%ld\n", a, b);}
            C_bd=U(a,b)+pow(((double)X)/(4.0*a), 1.0/3);
            for(c=(1-b); c<=C_bd;c++){
                P= b*b-3*a*c;
                D_lbd =(long)floor(-(((double)a-b)*(a-b+c))/a);
                D_ubd = (((double)a+b)*(a+b+c))/a;
                quadbd = a*(a-c);
                // printf("a=%ld, b=%ld, c=%ld, %ld<=d<=%lf\n", a, b, c, D_lbd, D_ubd);
                for(d=D_lbd+1; d<D_ubd; d++){
                    if(d*(d-b) < quadbd){continue;}//Improvement poss here, inc function
                    Q= b*c-9*a*d;
                    R= c*c-3*b*d;
                    D= Q*Q-4*P*R;
                    if(D>tX || D<=0){continue;} // IF DISC
                    f=gcd(P,gcd(Q,R));
                    if(D>3*B*f){continue;} //IF PRP
                    check=is_complex_field(a,0,c,d,P,Q,R,D,f,sqfull,B,list,pp,index);
                    if(check){
                    fprintf(fptr,"%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check);
                    }
                }
            }
        }
    }
    if(verbose){printf("...done.\n");}
    free(pp);
    free(sqfull);
    free(list);
    return 0;
}





//TODO: WRITE DISTRIBUTED VERSION for doing the n1 of n2 loads (break based on c variable?)
int CCFCCFPRP_distributed(long B, long n1, long n2, FILE *fptr, _Bool verbose){
    long a,b,c,d,D,f,P,Q,R,check,gcdPR,lowsplit,highsplit;
    long i=0, X=B*B, tX=3*X, lowsplit=((double)n1-1)/n2, highsplit=((double)n1)/n2;
    if(verbose){printf("Initialising...\n");}
    long* p = primes_up_to(3*B,sqrt(3*B));
    if(verbose){printf("...p initialised\n");}
    long* pp = (long*)malloc((p[0]+1)*sizeof(long));
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    if(verbose){printf("...pp initialised\n");}
    long index = init_index_find(p, 3*B, 3*B);
    if(verbose){printf("...dumping p...");}
    free(p);
    if(verbose){printf("done.\n");}
    if(verbose){printf("...index found\n");}
    long* sqfull = get_sqfull(pp,3*B);
    if(verbose){printf("...sqfull initialised\n");}
    long* list = get_list(pp, index, 3*B);
    if(verbose){printf("...list initialised\n");}
    if(verbose){printf("...done.\n");}
    // do looping
    double A_bd,B_bd,C_ubd,C_len,D_ubd;
    long C_lbd,D_lbd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    long quadbd;
    // b=0 looping
    if(verbose){printf("Working through b=0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        if(verbose){printf("...a=%ld\n", a);}
        C_bd=pow(((double)X)/(4*a), 1.0/3);
        for(c=1; c<=C_bd; c++){
            // if(verbose){printf("...a=%ld, c = %ld\n", a, c);}
            P= -3*a*c;
            R= c*c;
            D_lbd=0;
            gcdPR=gcd(P,R);
            if(c<=a){D_lbd = (long)floor(sqrt(a*(a-c)));}
            for(d=D_lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                Q= -9*a*d;
                D= Q*Q-4*Pb*R;
                if(D<=0){continue;}
                if(D>tX){break;} // check disc, note as b=0 then  once D gets big it only gets bigger as d increases for fixed a,c
                f=gcd(Q,gcdPR);
                if(D>3*B*f){continue;} // IF PRP
                check=is_complex_field(a,0,c,d,P,Q,R,D,f,sqfull,B,list,pp,index);
                if(check){
                    fprintf(fptr,"%ld,%d,%ld,%ld,%ld\n",a,0,c,d,check);
                }
            }
        }
    }
    if(verbose){printf("...done.\nWorking through b>0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
        for(b=1; b<=B_bd; b++){ 
            if(verbose){printf("...a=%ld, b=%ld\n", a, b);}
            C_len=U(a,b)+pow(((double)X)/(4.0*a), 1.0/3)-(1-b);
            if(lowsplit==0){
                C_lbd=1;
            }else{
                C_lbd=floor(1-b+lowsplit*C_len)+1;
            }            
            C_ubd=1-b+highsplit*C_len;
            for(c=(1-b); c<=C_bd;c++){
                P= b*b-3*a*c;
                D_lbd =(long)floor(-(((double)a-b)*(a-b+c))/a);
                D_ubd = (((double)a+b)*(a+b+c))/a;
                quadbd = a*(a-c);
                // printf("a=%ld, b=%ld, c=%ld, %ld<=d<=%lf\n", a, b, c, D_lbd, D_ubd);
                for(d=D_lbd+1; d<D_ubd; d++){
                    if(d*(d-b) < quadbd){continue;}//Improvement poss here, inc function
                    Q= b*c-9*a*d;
                    R= c*c-3*b*d;
                    D= Q*Q-4*P*R;
                    if(D>tX || D<=0){continue;} // IF DISC
                    f=gcd(P,gcd(Q,R));
                    if(D>3*B*f){continue;} //IF PRP
                    check=is_complex_field(a,b,c,d,P,Q,R,D,f,sqfull,B,list,pp,index);
                    if(check){
                    fprintf(fptr,"%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check);
                    }
                }
            }
        }
    }
    if(verbose){printf("...done.\n");}
    free(pp);
    free(sqfull);
    free(list);
    return 0;
}

int testing_suite(long a, long b, long c, long d, long B){
    long i=0;
    _Bool verbose=1;
    long P,Q,R,D,f;
    P= b*b-3*a*c;
    Q= b*c-9*a*d;
    R= c*c-3*b*d;
    D= Q*Q-4*P*R;
    f= gcd(P,gcd(Q,R));

    if(verbose){printf("Initialising...\n");}
    long* p = primes_up_to(3*B,sqrt(3*B));
    if(verbose){printf("...p initialised\n");}
    long* pp = (long*)malloc((p[0]+1)*sizeof(long));
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    if(verbose){printf("...pp initialised\n");}
    long index = init_index_find(p, 3*B, 3*B);
    if(verbose){printf("...index found\n");}
    long* sqfull = get_sqfull(pp,3*B);
    if(verbose){printf("...sqfull initialised\n");}
    long* list = get_list(pp, index, 3*B);
    if(verbose){printf("...list initialised\n");}
    if(verbose){printf("...done.\n");}
    //BELOW FOR LIVE TESTING
    printf("%ld\n", is_complex_field(a,b,c,d,P,Q,R,D,f,sqfull,B,list,pp,index));
    printf("%ld\n", pp[p[0]]);
    return 0;
}



int main(int argc, char** argv){
    FILE *fptr=fopen(argv[3], "w");
    // if(atoi(argv[1])==0){
    //     CCFCCF(pow(2,atol(argv[2])), fptr, 1);
    // }else if(atoi(argv[1])==1){
    //     CCFCCFPRP(pow(2,atol(argv[2])), fptr, 1);
    // }
    fclose(fptr);
}

