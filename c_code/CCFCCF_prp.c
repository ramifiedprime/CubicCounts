#include <stdio.h>
#include <stdlib.h>
#include <math.h>
////////////////////////////////////////////////////////////////////////////////
// math.h was being a pain so I just wrote it
// Input:
//      - int b: base
//      - int e: exponent
// Output:
//      - int ret: b^e
int power(int b, int e){
    int p=1;
    int i=0;
    while(i<e){
        p=p*b;
        i++;
    }
    return p;
}
////////////////////////////////////////////////////////////////////////////////
// Lazy implementation of Erastosthenes sieve
// Input:
//      - int X: integer bound
// Output:
//      - int* primes: array with primes[0] being the number of primes up to X
//                     and primes[i] being the ith prime for i>=1
int* primes_up_to(int X, int sqrtX){
    _Bool* comp = (_Bool*)calloc(X+1, sizeof(_Bool));//comp[n]=true iff n composite
    int i,j;
    int pi=0;
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
    int* primes = (int*)calloc(pi+1,sizeof(int));//extra entry for storing pi
    primes[0]=pi;
    j=1;
    for(i=2; i<= X; i++){
        if(!comp[i]){
            primes[j] = i;
            j++;
        }
    }
    free(comp);
    return primes; //first entry is number of primes, rest are primes
}

////////////////////////////////////////////////////////////////////////////////
// implementation of finding index in sub-algorithm 5.1 (init)
// Input:
//      - int* q: array of primes
//      - int tM: M from algorithm
// Output:
//      - int i: index of correct prime in list
int init_index_find(int* p, int X, int M){
    int r[2] = {1,p[0]};
    int i = (r[0]+r[1]+1)/2;
    int tM = 3*M;
    float num;
    while(i!=r[1]){
        num = X/(p[i]*log(p[i]))*(1+(1.0/(2*log(p[i]))));
        printf("%d, %f, %d, %d, %d, %d\n", i, num, tM, p[0], r[0], r[1]);
        if((X/log(p[i]))*(1+(0.5*log(p[i]))) <= tM){
            r[1]=i;
        }else{
            r[0]=i;
        }
        i = (r[0]+r[1]+1)/2;
    }
    if((X/log(p[i]))*(1+(0.5*log(p[i]))) > tM){
        printf("Error: Memory constraint too small for bound.\n");
        exit(1);
    }//At this stage we've completed the first pass
    if(p[i]<=53){
        float J;
        i = 16;
        while(J<=(tM/X - (1/log(59))*(1+(1/(2*log(59))))) && i>3){
            J+=1.0/(p[i]*p[i]);
            i--;
        }
    }
    return i;
}


////////////////////////////////////////////////////////////////////////////////
// init functions from Belabas, sub-algorithm 5.1
// Input:
//      - 
// Output:
//      - 
// int init(int X, int sqrtX, int M){
//     int* p = primes_up_to(X, sqrtX);
//     int* pp = (int*)calloc(p[0],sizeof(int))
//     int i;
//     for(i=1;i<=p[0]; i++){
//         pp[i]=p[i]*p[i];
//     }
//     // TEST 1
//     i=(q+1)/2;

// }





int main(void){
    int X=power(2,20);
    int sqrtX=power(2,10);
    int M=power(2,20);
    int* p = primes_up_to(X,sqrtX);
    printf("%d primes found\n", p[0]);
    int index = init_index_find(p, X, M);
    return 0;
}


// int write_data(int a, int b, int c, int d, int D, int f){
//     FILE* fptr;
//     fptr = fopen("test.dat", "w");
//     fprintf(fptr, "%d,%d,%d,%d,%d,%d\n", a,b,c,d,D,f);
//     fclose(fptr);
//     return 0;
// }
