/**
 * @file CCFCCFPRP.c
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
static inline long lmax(long a, long b){ return a > b ? a : b; }
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
 * @brief Takes in integer and returns the prime-to-6-part as well as data of whether it was divisible by 2 or 3
 * 
 * @param f=2^a3^bm for some m coprime to 6
 * @return m, f_isdivby2, f_isdivby3
 */
long primeto6part(long f){
    if(f==0){return f;}
    while(f%2==0){f/=2;}
    while(f%3==0){f/=3;}
    return f;
}

/**
 * @brief Takes in disriminant of the hessian, which is the absolute value of the discriminant of the proposed cubic, and the content of the Hessian f.  If the form represents a cubic field then this returns the product of ramified primes, else it just returns the product of ramified primes with extra exponents in places
 * 
 * @param D discriminant of Hessian
 * @param f_primeto6 prime to 6 part of content of the Hessian
 * @return product of ramified primes (or approximation if invalid form)
 */
__int128 get_prospective_prp(__int128 D, long f_primeto6){
    if(D==0){return 0;}
    __int128 prp=(D/3); // initialised as abs disc of field
    prp/=f_primeto6;// correct prp for primes >3 if corresp to field
    while(prp%9==0){prp/=3;} // fix p=3
    while(prp%4==0){prp/=2;} // fix p=2
    return prp;
}


/**
 * @brief Obtains all primes up to a bound 
 * 
 * Implements Erastosthenes sieve, no cleverness applied, to return all squares of primes where the primes are bounded by X.
 * 
 * @param X bound for the set of primes.
 * @return prime_squares array prime_squares with prime_squares[0] being the number of primes up to X
 *         and prime_squares[i] being the square of the ith prime for i>=1 who is less than X.
 */
long* prime_squares_up_to(long X){
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
    long* prime_squares = (long*)calloc(pi+1,sizeof(long));//extra entry for storing pi
    prime_squares[0]=pi;
    j=1;
    for(i=2; i<= X; i++){
        if(!comp[i]){
            prime_squares[j] = i*i;
            j++;
        }
    }
    free(comp);
    return prime_squares; //first entry is number of primes, rest are prime squares
}


/**
 * @brief Gets array `sqfull` from sub-algorithm 5.1 (`init`) 
 * 
 * Gets an array `sqfull` of length X+1 such that `sqfull[n]` is false
 * if and only if `n` is squarefree (away from 2 and 3).
 * 
 * @param X bound for array
 * @return array `sqfull` as required for `init`
 */
long* get_sqfull(long X){
    long i,n;
    long* pp = prime_squares_up_to(X);
    long* sqfull = (long*)calloc((X+1),sizeof(long));
    sqfull[0]=X;
    sqfull[1]=0;
    for(i=3; i<=pp[0] && pp[i]<=X; i++){
        for(n=pp[i]; n<=X; n+=pp[i]){
            sqfull[n]=1;
        }
    }
    free(pp);
    return sqfull;
}

/**
 * @brief recovers the quadratic discriminant associated to a complex cubic field given positive u which is squarefree away from 6 and is an integer such that squareroot of -u gives the resolvent quadratic.
 * 
 * @param u: integer
 * @return absolute value of quadratic discriminant.
 */
long get_quad_disc(long u){
    while(u%9 == 0){u/=9;} // get squarefree above 3
    while(u%4 == 0){u/=4;} // get squarefree above 2
    if (u%4 == 3){return u;}
    return 4*u;
}



/**
 * @brief Implements `is_complex_field` and `test` (sub-algorithms 5.6 and 5.3)
 * 
 * Tests belonging to U_2 and negative discriminant, then just runs test.
 * 
 * @param a First coefficient of cubic form.
 * @param b Second coefficient of cubic form.
 * @param c Third coefficient of cubic form.
 * @param d Fourth coefficient of cubic form.
 * @param P First coefficient of Hessian.
 * @param R Third coefficient of Hessian.
 * @param D Discriminant of Hessian(a,b,c,d).
 * @param f Content of Hessian(a,b,c,d).
 * @param sqfull array of length X+1 such that `sqfull[n]` is false if and only if it's squarefree away from 6.
 * @return test(a,b,c,d,D,f,sqfull).
 */
long is_complex_field(long a, long b, long c, long d, long P, long R, __int128 D, int f_isdivby3, long f_primeto6, long* sqfull){
    long t,u;
    // Local conditions at p>3:  these are quick lookups
    if(f_primeto6>sqfull[0]){fprintf(stderr, "Overflow on sqfull: sqfull[0]=%ld, f_primeto6=%ld\n", sqfull[0],f_primeto6);}
    if(sqfull[f_primeto6]){return 0;}
    // TODO: this will probably overflow?
    u = D/(3*f_primeto6*f_primeto6); // u is the thing we take sqrt of for resolvent quadratic, it will be the cubic disc minus primes in the content larger than 3.
    t = primeto6part(u);
    if(gcd(t,f_primeto6)!=1){return 0;}
    if(t>sqfull[0]){fprintf(stderr, "Overflow on sqfull: sqfull[0]=%ld, t=%ld\n", sqfull[0],t);}
    if (sqfull[t]){return 0;}

    // Local condition at 2
    if((D <= 0) || (D%16==0) || (D%16==4 && (P%2!=0 || R%2!=0))){
        return 0;
    }

    // Local condition at 3;
    long a9,d9;
    if ((u%9)==0 && (!f_isdivby3)){return 0;} // Not in V3 and not (1^3)
    else if(f_isdivby3){// checking if in U3 given that it's 1^3
        a9 = a%9;
        d9 = d%9;
        if(a9==0 || d9==0){return 0;}
        else if(((a9%3)==0)&&((d9%3)==0)){return 0;}
        else if((((a9-d9)%3)==0) && ((a9-b+c-d9)%9)==0){return 0;}
        else if((((a9+d9)%3)==0) && ((a9+b+c+d9)%9)==0){return 0;}
    }
    // Checks all passed, it's a field, return the quadratic resolvent disc.  Since the checks have passed, u is squarefree away from 
    return get_quad_disc(u);
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
 * @brief Recovers cubic fields of bounded product of ramified primes (adapts algorithm 5.7) with fixed first two coefficients (a,b) where b is nonzero.
 * 
 * Implements Belabas' algorithm, but makes cuts for product of ramified primes
 * ordering.  The results are printed to a file, where results take the form of rows:
 *      a,b,c,d,|Disc(L)|,prp\n
 * where ax^3+bx^2y+cxy^2+dy^3 is a reduced binary cubic form corresponding 
 * to a complex cubic field, D is the discriminant of the Hessian, f is its content, and L is the
 * associated quadratic resolvent, and prp is the product of ramified primes.
 * This function obtains all such entries with fixed (a,b) and product of ramified primes at most B.
 * 
 * @param a first coefficient of the binary cubic form
 * @param b second coefficient of the binary cubic form
 * @param sqfull long list where sqfull[n] is false iff n is sqfree away from 6
 * @param X bound for discriminant
 * @param tX 3X 
 * @param B product of ramified primes bound for cubic fields.
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCFPRP_fixed_ab(long a, long b, long* sqfull, long X, long tX, long B, FILE *fptr, int verbose){
    int f_isdivby3;
    long c, d, f, f_primeto6, P, Q, R, check, gcdPR;
    // 54 is to account for wild ramification at 2 and 3
    // int128 is to account for B around 30 bits
    __int128 D, prp, D_cut;
    long D_base, quadbd;
    double C_bd, d_lin_ubd, temp, tempdisc, d_disc_lo, d_disc_hi;
    double alpha_d,beta_d,gamma_d;
    long d_lin_lbd, d_lbd, quad_hi, quad_lo, lo, hi, hi1, lo1;
    if(b==0){
        C_bd=pow(((double)X)/(4*a), 1.0/3);
        for(c=1; c<=C_bd; c++){
            // if(verbose){printf("...a=%ld, c = %ld\n", a, c);}
            P= -3*a*c;
            R= c*c;
            D_base= -4*P*R;
            d_lbd=0;
            gcdPR=gcd(P,R);
            D_cut=18*B*gcdPR;
            if(c<=a){d_lbd = (long)floor(sqrt(a*(a-c)));}
            for(d=d_lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                Q= -9*a*d;
                D= (__int128)Q*Q+(__int128)D_base;
                if(D>D_cut){break;} // check disc, note as b=0 then once D gets big it only gets bigger as d increases for fixed a,c
                f=gcd(gcdPR, Q);
                f_isdivby3 = (f%3 == 0);
                f_primeto6=primeto6part(f);
                prp=get_prospective_prp(D,f_primeto6);
                if(prp>B){continue;} // PRP bound
                check=is_complex_field(a,0,c,d,P,R,D,f_isdivby3,f_primeto6,sqfull);
                if(check){
                    fprintf(fptr,"%ld,%d,%ld,%ld,%ld,%ld\n",a,0,c,d,check,(long)prp); // output is a,b,c,d,quadratic disc, prp
                }
            }
        }
    }
    else{
        C_bd=U(a,b)+pow(((double)X)/(4.0*a), 1.0/3);
        for(c=(1-b); c<=C_bd; c++){
            P = b*b - 3*a*c;

            // lemma 4.2 linear bounds
            d_lin_lbd = (long)floor(-(((double)a-b)*(a-b+c))/a);
            d_lin_ubd = (((double)a+b)*(a+b+c))/a;

            // disc(d) = alpha*d^2 + beta*d + gamma <= tX
            // where alpha=81a^2, beta=12b^3-54abc, gamma=12ac^3-3b^2c^2.
            // disc(d) <= tX for d in [d_disc_lo, d_disc_hi]
            alpha_d  = 81.0*(double)a*a;
            beta_d   = 12.0*(double)b*b*b - 54.0*(double)a*b*c;
            gamma_d  = 12.0*(double)a*c*c*c - 3.0*(double)b*b*c*c;
            tempdisc = beta_d*beta_d - 4.0*alpha_d*(gamma_d - (double)tX);
            if(tempdisc < 0){ continue; } // D(d)>tX for all d as parabola doesn't meet axis
            temp = sqrt(tempdisc);
            d_disc_hi = (-beta_d + temp) / (2.0*alpha_d);
            d_disc_lo = (-beta_d - temp) / (2.0*alpha_d);

            // intersect linear and discriminant bounds
            lo = lmax((long)ceil(d_disc_lo) - 1, d_lin_lbd + 1);
            hi = (long)fmin(floor(d_disc_hi) + 1.0, floor(d_lin_ubd - 1e-9));
            if(lo > hi){ continue; }

            // Quadratic condition: d^2 - b*d - a*(a-c) >= 0
            // holds outwith the roots:
            //   d <= root_lo  OR  d >= root_hi
            // where root_lo/hi = (b -/+ sqrt(b^2+4a(a-c)))/2.
            // If form_disc = b^2+4a(a-c) < 0: no real roots, quadratic always positive,
            // condition always satisfied — single interval [lo, hi].
            // If form_disc >= 0: two valid sub-intervals, run two separate loops. 
            quadbd = a*(a-c);
            tempdisc = (double)b*b + 4.0*(double)quadbd;

            if(tempdisc < 0){
                // quadratic condition trivial: d is in [lo, hi]
                for(d=lo; d<=hi; d++){
                    Q = b*c - 9*a*d;
                    R = c*c - 3*b*d;
                    D = (__int128)Q*Q - (__int128)4*P*R;
                    if(D <= 0 || D > tX){ continue; } // safety 
                    f = gcd(P, gcd(Q,R));
                    f_isdivby3 = (f%3 == 0);
                    f_primeto6=primeto6part(f);
                    prp=get_prospective_prp(D,f_primeto6);
                    if(prp>B){continue;} // IF PRP
                    check = is_complex_field(a,b,c,d,P,R,D,f_isdivby3,f_primeto6,sqfull);
                    if(check){
                        fprintf(fptr,"%ld,%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check,(long)prp);
                    }
                }
            } else {
                // quadratic condition is nontrivial
                temp   = sqrt(tempdisc);
                quad_hi = (long)floor((b - temp) / 2.0);
                quad_lo = (long)ceil((b + temp) / 2.0);
                // lower interval: [lo, min(hi, quad_hi1)]
                hi1 = (quad_hi < hi) ? quad_hi : hi;
                for(d=lo; d<=hi1; d++){
                    if(d*(d-b) < quadbd){ continue; } // safety
                    Q = b*c - 9*a*d;
                    R = c*c - 3*b*d;
                    D = (__int128)Q*Q - (__int128)4*P*R;
                    if(D <= 0 || D > tX){ continue; } // safety
                    f = gcd(P, gcd(Q,R));
                    f_isdivby3 = (f%3 == 0);
                    f_primeto6=primeto6part(f);
                    prp=get_prospective_prp(D,f_primeto6);
                    if(prp>B){continue;}
                    check = is_complex_field(a,b,c,d,P,R,D,f_isdivby3,f_primeto6,sqfull);
                    if(check){
                        fprintf(fptr,"%ld,%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check,(long)prp);
                    }
                }

                // upper interval: [max(lo, quad_lo), hi]
                lo1 = (quad_lo > lo) ? quad_lo : lo;
                for(d=lo1; d<=hi; d++){
                    if(d*(d-b) < quadbd){ continue; } //safety net
                    Q = b*c - 9*a*d;
                    R = c*c - 3*b*d;
                    D = (__int128)Q*Q - (__int128)4*P*R;
                    if(D <= 0 || D > tX){ continue; } // safety net
                    f = gcd(P, gcd(Q,R));
                    f_isdivby3 = (f%3 == 0);
                    f_primeto6=primeto6part(f);
                    prp=get_prospective_prp(D,f_primeto6);
                    if(prp>B){continue;} // IF PRP
                    check = is_complex_field(a,b,c,d,P,R,D,f_isdivby3,f_primeto6,sqfull);
                    if(check){
                        fprintf(fptr,"%ld,%ld,%ld,%ld,%ld,%ld\n",a,b,c,d,check,(long)prp);
                    }
                }
            }
        }
    }
    return 0;
}






/**
 * @brief Recovers cubic fields of bounded product of ramified primes (adapts algorithm 5.7) with fixed first coefficient a.
 * 
 * Implements Belabas' algorithm, but makes cuts for product of ramified primes
 * ordering.  The results are printed to a file, where results take the form of rows:
 *      a,b,c,d,|Disc(L)|,D/-3f\n
 * where ax^3+bx^2y+cxy^2+dy^3 is a reduced binary cubic form corresponding 
 * to a complex cubic field, D is the discriminant of the Hessian, f is its content, and L is the
 * associated quadratic resolvent.  Note that D/-3f is (away from 2) the product of ramified primes.  
 * This function obtains all such entries with fixed a and product of ramified primes at most B.
 * 
 * @param a first coefficient of the binary cubic form
 * @param sqfull long list where sqfull[n] is false iff n is sqfree away from 6
 * @param X bound for discriminant
 * @param tX 3X 
 * @param B product of ramified primes bound for cubic fields.
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCFPRP_fixeda(long a, long* sqfull, long X, long tX, long B, FILE *fptr, int verbose){
    long b;
    double B_bd;
    if(verbose>=2){printf("...a=%ld, b=0\n", a);}
    B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
    for(b=0; b<=B_bd; b++){
        if(verbose>=2){printf("...a=%ld, b=%ld\n", a, b);}
        CCFCCFPRP_fixed_ab(a,b,sqfull, X, tX, B,fptr,verbose);
    }
    return 0;
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
 * This runs CCFCCFPRP_fixeda as the subroutine which does most work.
 * 
 * @param B product of ramified primes bound for cubic fields.
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCFPRP(long B, FILE *fptr, int verbose){
    if(verbose>=2){printf("Initialising...\n");}
    long* sqfull = get_sqfull(13*B);
    if(verbose>=2){printf("...sqfull initialised\n");}
    // do looping
    long X=54*B*B, tX=3*X, a;
    double A_bd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    for(a=1;a<=A_bd;a++){
        CCFCCFPRP_fixeda(a,sqfull, X, tX, B,fptr,verbose);
    }
    if(verbose>=1){printf("...done.\n");}
    free(sqfull);
    return 0;
}



/**
 * @brief distributable version of CCFCCFPRP
 * 
 * Implements a sub-search of CCFCCF(B, *fptr, verbose), useful to allow running many searches in parallel.
 * Parallelisation is on (a,b) paits, with a round-robin approach (so the n1th out of every n2).
 * 
 * @param B product of ramified primes bound for cubic fields.
 * @param n1 integer between 1 and n2
 * @param n2 integer
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */

int CCFCCFPRP_distributed(long B, long n1, long n2, FILE *fptr, int verbose){
    if(verbose>=2){printf("Initialising...\n");}
    long* sqfull = get_sqfull(13*B);
    if(verbose>=2){printf("...sqfull initialised\n");}
    
    // do looping
    long X=54*B*B, tX=3*X, a, b, idx=0;
    double A_bd,B_bd;
    A_bd=pow((16.0*X)/27, 1.0/4);
    for(a=1;a<=A_bd;a++){
        B_bd=(3.0*(double)a)/2 + sqrt(sqrt(((double)X)/3) - (3.0*a*a)/4);
        for(b=0; b<=B_bd; b++){
            if(verbose>=2){printf("...a=%ld, b=%ld\n", a, b);}
            if(idx%n2==n1-1){
                CCFCCFPRP_fixed_ab(a,b,sqfull, X, tX, B,fptr,verbose);
            }
            idx++;
        }
    }
    if(verbose>=1){printf("...done.\n");}
    free(sqfull);
    return 0;
}