/**
 * @file CCFCCF.c
 * @brief Finds all cubic fields with bounded product of ramified primes.
 *
 * @details Implementation of Belabas' algorithm in the article _a fast algorithm to compute cubic fields_
 * with some minor modifications to count according to product of ramified primes.  Any references to 
 * algorithm numbers correspond to those in the published version of loc. cit..
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * @brief Greatest common divisor of two integers.
 * 
 * @param a First integer.
 * @param b Second integer.
 * @return greatest common divisor of a and b.
 */
long gcd(long a, long b){
    long r;
    while(b!=0){
        r=a%b;
        a=b;
        b=r;
    }
    return labs(a);
}

/**
 * @brief Obtains all primes up to a bound 
 * 
 * Implements Erastosthenes sieve, no cleverness applied, to return all primes
 * up to given bound.
 * 
 * @param X bound for the set of primes.
 * @return primes array primes with primes[0] being the number of primes up to X
 *         and primes[i] being the ith prime for i>=1.
 */
long* primes_up_to(long X){
    // printf("Initialising p...\n");
    long sqrtX = sqrt(X);
    //Precision error in sqrt capture:
    while (sqrtX * sqrtX > X) sqrtX--;      // round up correction
    while ((sqrtX+1) * (sqrtX+1) <= X) sqrtX++;  // round down correction
    _Bool* comp = (_Bool*)calloc(X+1, sizeof(_Bool));//comp[n]=true iff n composite
    long i,j;
    long pi=0;
    for(i=2; i<= X; i++){
        if(!comp[i]){
            pi++;
            if(i<=sqrtX){
                for(j=i; j<= X/i; j++){
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


/**
 * @brief Gets array `sqfull` from sub-algorithm 5.1 (`init`) 
 * 
 * Gets an array `sqfull` of length X+1 such that `sqfull[n]` is false
 * if and only if `n` is squarefree (away from 2 and 3).
 * 
 * @param pp squares of primes at most X
 * @param X bound for array
 * @return array `sqfull` as required for `init`
 */
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

/**
 * @brief recovers the quadratic discriminant associated to a complex cubic field
 * 
 * Given a reduced cubic form F associated to a cubic field K by Delone--Fadeev
 * correspondence, let D be the discriminant of the Hessian and f its content.  Recall that D=-3*Disc(F) and f^2|D.
 * Then we use D/f^2 to determine the associated quadratic discriminant.  t is the prime-to-6 part of this fraction, which
 * is necessarily the squarefree prime-to-6 part of the quadratic discriminant.  u is the gcd of this fraction and 72, which
 * precisely determines the 2 and 3 parts of the discriminant by standard arguments. 
 * 
 * Indeed, precisely $K=\mathbb{Q}(sqrt{-ut/3})$
 * 
 * @param t prime-to-6 squarefree part of D
 * @param u Greatest common divisor of D/f^2 and 3^5x2^4.
 * @return Associated quadratic discriminant.
 */
long get_quad_disc(long t, long u){
    t=t*3*u;
    while(t%9==0){
        t=t/9;
    }
    while(t%4==0){
        t=t/4;
    }
    if(t%4!=3){
        return 4*t;
    }
    return t;
}



/**
 * @brief tests whether a reduced cubic form corresponds to a field, implements `test` (sub-algorithm 5.2) with a modification for working in the product of ramified primes setting:  
 *  We call get_sqfull(pp,3*B), so that the range of the squarefull numbers in sqfull is up to 3B. Our input t will satisfy
 *     t=D/(f^2u)<=D/f^2<=3Bf/f^2<= 3B,
 * we can always just do the lookup in the sqfull table.
 * 
 * Assumes that the input is in U_2 already (i.e. satisfies 2-adic constraint)
 * Tests belonging to U_3 by congruences, then checks minmality at other primes
 * by checking if D/f^2 is squarefree and coprime to f away from 6.
 * 
 * returns 0 at any failure, returns the associated quadratic discriminant if true.
 * 
 * @param a First coefficient of cubic form.
 * @param b Second coefficient of cubic form.
 * @param c Third coefficient of cubic form.
 * @param d Fourth coefficient of cubic form.
 * @param D Discriminant of Hessian(a,b,c,d).
 * @param f Content of Hessian(a,b,c,d).
 * @param sqfull array of length X+1 such that `sqfull[n]` is false if and only if it's squarefree away from 6.
 * @return quadratic disriminant of the associated field, or 0 if it fails to correspond.
 */
long test(long a, long b, long c, long d, long D, long f, long* sqfull){
    long a9,d9,t,u;
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
    u = gcd(t,3888);
    t = t/u;
    if(gcd(t,f)!=1){return 0;}
    if (sqfull[t]){return 0;}
    return get_quad_disc(t, u);
}

/**
 * @brief Implements `is_complex_field` (sub-algorithm 5.6)
 * 
 * Tests belonging to U_2 and negative discriminant, then just runs test.
 * 
 * @param a First coefficient of cubic form.
 * @param b Second coefficient of cubic form.
 * @param c Third coefficient of cubic form.
 * @param d Fourth coefficient of cubic form.
 * @param P First coefficient of Hessian.
 * @param Q Second coefficient of Hessian.
 * @param R Third coefficient of Hessian.
 * @param D Discriminant of Hessian(a,b,c,d).
 * @param f Content of Hessian(a,b,c,d).
 * @param sqfull array of length X+1 such that `sqfull[n]` is false if and only if it's squarefree away from 6.
 * @return test(a,b,c,d,D,f,sqfull).
 */
long is_complex_field(long a, long b, long c, long d, long P, long Q, long R, long D, long f,long* sqfull){
        if((D <= 0) || (D%16==0) || (D%16==4 && (P%2!=0 || R%2!=0))){
            return 0;
        }
        // printf("\tpassed U2 check, now running test\n");
        return test(a,b,c,d,D,f,sqfull);
}

/**
 * @brief U function from lemma 4.4
 * 
 * @param a First coefficient of cubic form
 * @param b Second coefficient of cubic form.
 * @return U(a,b)
 */
double U(long a, long b){
    if(a>=2.0*b/3){return ((double)b*b)/(3.0*a);}
    else{return b-3*((double)a)/4.0;}
}

/**
 * @brief Recovers cubic fields of bounded product of ramified primes (adapts algorithm 5.7).
 * 
 * Implements Belabas' algorithm, but makes cuts for product of ramified primes
 * ordering.  The results are printed to a file, where results take the form of rows:
 *      a,b,c,d,|Disc(L)|,D/-3f\n
 * where ax^3+bx^2y+cxy^2+dy^3 is a reduced binary cubic form corresponding 
 * to a complex cubic field, D is the discriminant of the Hessian, f is its content, and L is the
 * associated quadratic resolvent.  Note that D/-3f is (away from 2) the product of ramified primes.  
 * This function obtains all such entries with product of ramified primes at most X.
 * 
 * @param B product of ramified primes bound for cubic fields.
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCFPRP(long B, FILE *fptr, int verbose){
    long a,b,c,d,f,P,Q,R,check,gcdPR;
    long i=0;
    // 54 is to account for wild ramification at 2 and 3
    // int128 is to account for B around 30 bits
    __int128 D;
    long D_cut, D_base, X=54*B*B, tX=3*X; 
    if(verbose>=2){printf("Initialising...\n");}
    long* p = primes_up_to(3*B);
    if(verbose>=2){printf("...p initialised\n");}
    long* pp = (long*)malloc((p[0]+1)*sizeof(long));
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    if(verbose>=2){printf("...pp initialised\n");}
    if(verbose>=2){printf("...dumping p...");}
    free(p);
    if(verbose>=2){printf("done.\n");}
    long* sqfull = get_sqfull(pp,3*B);
    if(verbose>=2){printf("...sqfull initialised\n");}
    if(verbose>=2){printf("...done.\n");}
    // do looping
    double A_bd,B_bd,C_bd,D_ubd;
    long D_lbd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    long quadbd;
    // b=0 looping
    if(verbose>=2){printf("Working through b=0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        if(verbose>=2){printf("...a=%ld\n", a);}
        C_bd=pow(((double)X)/(4*a), 1.0/3);
        for(c=1; c<=C_bd; c++){
            // if(verbose){printf("...a=%ld, c = %ld\n", a, c);}
            P= -3*a*c;
            R= c*c;
            D_base= -4*P*R;
            D_lbd=0;
            gcdPR=gcd(P,R);
            D_cut=3*B*gcdPR;
            if(c<=a){D_lbd = (long)floor(sqrt(a*(a-c)));}
            for(d=D_lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                Q= -9*a*d;
                D= (__int128)Q*Q+(__int128)D_base;
                if(D>D_cut){break;} // check disc, note as b=0 then  once D gets big it only gets bigger as d increases for fixed a,c
                f=gcd(gcdPR, Q);
                if(D>3*B*f){continue;} // IF PRP
                check=is_complex_field(a,0,c,d,P,Q,R,(long)D,f,sqfull);
                if(check){
                    fprintf(fptr,"%ld,%d,%ld,%ld,%ld,%ld\n",a,0,c,d,check,(long)D/(3*f)); // output is a,b,c,d,quadratic
                }
            }
        }
    }
    if(verbose>=2){printf("...done.\nWorking through b>0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
        for(b=1; b<=B_bd; b++){ 
            if(verbose>=2){printf("...a=%ld, b=%ld\n", a, b);}
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
                    D= (__int128)Q*Q-(__int128)4*P*R;
                    if(D>tX || D<=0){continue;} // IF DISC
                    f=gcd(P,gcd(Q,R));
                    if(D>3*B*f){continue;} //IF PRP
                    check=is_complex_field(a,b,c,d,P,Q,R,(long)D,f,sqfull);
                    if(check){
                    fprintf(fptr,"%ld,%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check,(long)D/(3*f));
                    }
                }
            }
        }
    }
    if(verbose>=1){printf("...done.\n");}
    free(pp);
    free(sqfull);
    return 0;
}



/**
 * @brief distributable version of CCFCCFPRP
 * 
 * Implements a sub-search of CCFCCF(B, *fptr, verbose), useful to allow running many searches in parallel.
 * Splits the range for the coefficient c into subranges of relative size 1/n2, and runs the n1-th subrange.
 * 
 * @param B product of ramified primes bound for cubic fields.
 * @param n1 integer between 1 and n2
 * @param n2 integer
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCFPRP_distributed(long B, long n1, long n2, FILE *fptr, int verbose){
    long a,b,c,d,f,P,Q,R,check,gcdPR;
    long i=0;
    __int128 D;
    long D_cut, D_base, X=54*B*B, tX=3*X;
    double lowsplit=((double)n1-1)/n2;
    double highsplit=((double)n1)/n2;
    if(verbose>=2){printf("Initialising...\n");}
    long* p = primes_up_to(3*B);
    if(verbose>=2){printf("...p initialised\n");}
    long* pp = (long*)malloc((p[0]+1)*sizeof(long));
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    if(verbose>=2){printf("...pp initialised\n");}
    if(verbose>=2){printf("...dumping p...");}
    free(p);
    if(verbose>=2){printf("done.\n");}
    long* sqfull = get_sqfull(pp,3*B);
    if(verbose>=2){printf("...sqfull initialised\n");}
    if(verbose>=2){printf("...done.\n");}
    // do looping
    double A_bd,B_bd,C_ubd,C_len,D_ubd;
    long C_lbd,D_lbd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    long quadbd;
    // b=0 looping
    if(verbose>=2){printf("Working through b=0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        if(verbose>=2){printf("...a=%ld\n", a);}
        C_len=pow(((double)X)/(4*a), 1.0/3);
        C_lbd=floor(lowsplit*C_len)+1;
        C_ubd=highsplit*C_len;
        for(c=C_lbd; c<=C_ubd; c++){
            // if(verbose){printf("...a=%ld, c = %ld\n", a, c);}
            P= -3*a*c;
            R= c*c;
            D_base= -4*P*R;
            D_lbd=0;
            gcdPR=gcd(P,R);
            D_cut=3*B*gcdPR;
            if(c<=a){D_lbd = (long)floor(sqrt(a*(a-c)));}
            for(d=D_lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                Q= -9*a*d;
                D= (__int128)D_base + (__int128)Q*Q;
                if(D>D_cut){break;} // check disc, note as b=0 then  once D gets big it only gets bigger as d increases for fixed a,c
                f=gcd(Q,gcdPR);
                if(D>3*B*f){continue;} // IF PRP
                check=is_complex_field(a,0,c,d,P,Q,R,(long)D,f,sqfull);
                if(check){
                    fprintf(fptr,"%ld,%d,%ld,%ld,%ld,%ld\n",a,0,c,d,check,(long)D/(3*f));
                }
            }
        }
    }
    if(verbose>=2){printf("...done.\nWorking through b>0 cases...\n");}
    for(a=1;a<=A_bd;a++){
        B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
        for(b=1; b<=B_bd; b++){ 
            if(verbose>=2){printf("...a=%ld, b=%ld\n", a, b);}
            C_len=U(a,b)+pow(((double)X)/(4.0*a), 1.0/3)+b;
            C_lbd=floor(lowsplit*C_len-b)+1;
            C_ubd=highsplit*C_len-b;
            for(c=C_lbd; c<=C_ubd;c++){
                P= b*b-3*a*c;
                D_lbd =(long)floor(-(((double)a-b)*(a-b+c))/a);
                D_ubd = (((double)a+b)*(a+b+c))/a;
                quadbd = a*(a-c);
                // printf("a=%ld, b=%ld, c=%ld, %ld<=d<=%lf\n", a, b, c, D_lbd, D_ubd);
                for(d=D_lbd+1; d<D_ubd; d++){
                    if(d*(d-b) < quadbd){continue;}// Improvement poss here, inc function
                    Q= b*c-9*a*d;
                    R= c*c-3*b*d;
                    D= (__int128)Q*Q-(__int128)4*P*R;
                    if(D>tX || D<=0){continue;} // IF DISC
                    f=gcd(P,gcd(Q,R));
                    if(D>3*B*f){continue;} //IF PRP
                    check=is_complex_field(a,b,c,d,P,Q,R,(long)D,f,sqfull);
                    if(check){
                    fprintf(fptr,"%ld,%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check,(long)D/(3*f));
                    }
                }
            }
        }
    }
    if(verbose>=1){printf("...done with run %ld of %ld.\n", n1, n2);}
    free(pp);
    free(sqfull);
    return 0;
}
