#ifndef CCFCCFPRP_H
#define CCFCCFPRP_H
#include<stdio.h>

/**
 * @brief Recovers cubic fields of bounded product of ramified primes (adapts algorithm 5.7).
 * 
 * Implements Belabas' algorithm, but makes cuts for product of ramified primes
 * ordering.  Similar to CCFCCF this prints results to a file.  The printed
 * results take the form of lines:
 *      a,b,c,d,Disc(L),D/f^2\n
 * where ax^3+bx^2y+cxy^2+dy^3 is a reduced binary cubic form corresponding 
 * to a complex cubic field, D is the discriminant of the Hessian, f is its content, and L is the
 * associated quadratic resolvent (essentially product of ramified primes, away from 2).  
 * This function obtains all such entries with product of ramified primes at most X.
 * 
 * @param B product of ramified primes bound for cubic fields.
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCFPRP(long B, FILE *fptr, int verbose);
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
int CCFCCFPRP_distributed(long B, long n1, long n2, FILE *fptr, int verbose);

#endif //CCFCCFPRP_H