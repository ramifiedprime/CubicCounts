#include <stdio.h>
#include <stdlib.h>
#include <math.h>

////////////////////////////////////////////////////////////////////////////////
// Lazy implementation of Erastosthenes sieve
// Input:
//      - int X: integer bound
// Output:
//      - int* primes: array with primes[0] being the number of primes up to X
//                     and primes[i] being the ith prime for i>=1
int* primes_up_to(int X){
    _Bool* comp = (_Bool*)malloc(X*sizeof(_Bool));//comp[n]=true iff n composite
    int i,j;
    int pi = 0;
    for(i=2; i<= X; i++){
        if(!comp[i]){
            pi++;
            for(j=2; j<= X/i; j++){
                comp[i*j] = 1;
            }
        }
    }
    int* primes = (int*)malloc((pi+1)*sizeof(int));//extra entry for storing pi
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



int write_data(int a, int b, int c, int d, int D, int f){
    FILE* fptr;
    fptr = fopen("test.dat", "w");
    fprintf(fptr, "%d,%d,%d,%d,%d,%d\n", a,b,c,d,D,f);
    fclose(fptr);
    return 0;
}



int main(void){
    int X=pow(2,30);
    int* primes = primes_up_to(X);
    printf("%d primes found\n", primes[0]);
    // for(int i=1; i<=primes[0]; i++){
    //     printf("%d ", primes[i]);
    // }
    printf("\n");
    free(primes);
    return 0;
}