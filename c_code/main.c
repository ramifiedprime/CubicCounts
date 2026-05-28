/**
 * @file main.c
 * @brief Entry point for CCFCCFPRP, also supports distributed modes.
 *
 * Recommended compilation:
 * @code
 * gcc main.c CCFCCFPRP.c -lm -o main.exe
 * @endcode
 * 
 * Standard usage:
 * @code
 *   ./main.exe <exponent> <outfile>
 * @endcode  
 * e.g. ./main.exe 20 output.dat
 * runs CCFCCFPRP with B=2^exponent, writes to outfile.
 * 
 * Distributed usage:
 * @code
 *   ./main.exe <exponent> <n1> <n2> <outfile>
 * @endcode
 * e.g. ./main.exe 25 3 100 output.dat
 * runs CCFCCFPRP_distributed with B=2^exponent, writes to outfile.
 * Note:  this will run the CCRCCRPRP algorithm but only for leading coefficient a=n1 mod n2
 */
#include"CCFCCFPRP.h"
#include<math.h>
#include<stdlib.h>
/**
 * @brief main runs CCFCCFPRP in standard or distributed mode.

 * Useage is one of the two types below:
 * @code
 *   ./main.exe <exponent> <outfile>
 *   ./main.exe <exponent> <n1> <n2> <outfile>
 * @endcode
 * 
 * @param exponent the bound for the number of bits of the product of ramified primes, i.e. B=2^exponent for CCFCCFPRP.
 * @param n1 integer between 1 and n2
 * @param n2 integer
 * @param outfile name of output file.
 * @return 0
 */
int main(int argc, char** argv){
    if(argc == 3){ // standard mode
        FILE *fptr=fopen(argv[2], "w");
        long exponent = atol(argv[3]);
        printf("CCFCCFPRP output: B=2^%ld\n", exponent);
        CCFCCFPRP(pow(2,exponent), fptr, 0);
        fclose(fptr);
    }
    else if(argc == 5){// distributed mode
        FILE *fptr=fopen(argv[4], "w");
        long n1 = atol(argv[1]);
        long n2 = atol(argv[2]);
        long exponent = atol(argv[3]);
        if(n1 > n2 || n1 < 1){fprintf(stderr, "Error: requires 1 <= n1 <= n2"); return 1;}
        fprintf(fptr, "CCFCCFPRP distributed output: B=2^%ld, task %ld of %ld\n", exponent, n1, n2);
        CCFCCFPRP_distributed(pow(2,exponent), n1, n2, fptr, 0);
        fclose(fptr);
    }
    else{
        fprintf(stderr, "Error: incorrect number of inputs.");
        return 1;
    }
    return 0;
}
