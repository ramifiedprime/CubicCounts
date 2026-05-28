/**
 * @file main.c
 * @brief Entry point for CCFCCFPRP, also supports distributed modes.
 *
 * Recommended compilation:
 * @code
 * gcc -Wall main.c CCFCCFPRP.c -lm -o main.exe
 * @endcode
 * 
 * Standard usage:
 * @code
 *   ./main.exe <B> <outfile>
 *   ./main.exe -b <exponent> <outfile>
 * @endcode  
 * runs CCFCCFPRP with B=2^exponent, writes to outfile.
 * 
 * Distributed usage:
 * @code
 *   ./main.exe <B> <n1> <n2> <outfile>
 *   ./main.exe -b <exponent> <n1> <n2> <outfile>
 * @endcode
 * e.g. ./main.exe 25 3 100 output.dat
 * runs CCFCCFPRP_distributed with B=2^exponent, writes to outfile.
 */
#include"CCFCCFPRP.h"
#include<math.h>
#include<stdlib.h>
#include<unistd.h>

/**
 * @brief main runs CCFCCFPRP in standard or distributed mode.
 * 
 * @return 0
 */
long get_B(int bits_flag, char* input){
    if(bits_flag){return (long)pow(2,atol(input));}
    else{return atol(input);}
}

int main(int argc, char** argv){
    int bits_flag = 0;
    int opt;
    while((opt = getopt(argc, argv, "b")) != -1){
        if(opt == 'b'){bits_flag=1;}
        else{fprintf(stderr, "Error: unknown flag.\nUsage: %s [-b] <B_or_exponent> <n1> <n2> <outfile>\n", argv[0]); return 1;}
    }
    int n_inputs = argc-optind;
    if(n_inputs == 2){// standard mode
        if(bits_flag){printf("CCFCCFPRP: B=2^%s writing to %s\n", argv[optind], argv[optind+1]);}
        else{printf("CCFCCFPRP: B=%s writing to %s\n", argv[optind], argv[optind+1]);}
        FILE *fptr=fopen(argv[optind+1], "w");
        if(!fptr){fprintf(stderr, "Error: could not open %s\n", argv[optind+1]); return 1;}
        CCFCCFPRP(get_B(bits_flag, argv[optind]), fptr, 0);
        fclose(fptr);
    }
    else if(n_inputs == 4){// distributed mode
        long n1 = atol(argv[optind+1]);
        long n2 = atol(argv[optind+2]);
        if(n1 > n2 || n1 < 1){fprintf(stderr, "Error: requires 1 <= n1 <= n2\n"); return 1;}
        if(bits_flag){printf("CCFCCFPRP distributed: B=2^%s, task %ld of %ld, writing to %s\n", argv[optind], n1, n2, argv[optind+3]);}
        else{printf("CCFCCFPRP distributed: B=%s, task %ld of %ld, writing to %s\n", argv[optind], n1, n2, argv[optind+3]);}
        FILE *fptr=fopen(argv[optind+3], "w");
        if(!fptr){fprintf(stderr, "Error: could not open %s\n", argv[optind+3]); return 1;}
        CCFCCFPRP_distributed(get_B(bits_flag, argv[optind]), n1, n2, fptr, 0);
        fclose(fptr);
    }
    else{
        fprintf(stderr, "Error: incorrect number of inputs.\n");
        return 1;
    }
    return 0;
}
