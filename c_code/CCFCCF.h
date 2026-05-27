#ifndef CCFCCF_H
#define CCFCCF_H
#include<stdio.h>


/**
 * @brief Recovers cubic fields of bounded discriminant (algorithm 5.7).
 * 
 * Implements Belabas' algorithm, printing results to a file.  The printed
 * results take the form of lines:
 *      a,b,c,d\n
 * where ax^3+bx^2y+cxy^2+dy^3 is a reduced binary cubic form corresponding 
 * to a complex cubic field.  Obtains all such forms of discriminant at most X.
 * 
 * @param X Discriminant bound for cubic fields.
 * @param *fptr file to output results to.
 * @param verbose verbosity flag for bug-tests.
 * @return 0
 */
int CCFCCF(long X, FILE *fptr, int verbose);
#endif //CCFCCF_H