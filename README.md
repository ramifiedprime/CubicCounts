# Cubic Counts
Repository for some code to count cubic fields ordered by product of ramified primes.  Essentially just implementing the algorithm in Belabas' article [1], which counts by discriminant, and modifying in appropriate places to get the required count.  Only interested in complex cubics for now, since that's where my practical interest comes in 

## Files
### magma\_code
- *CCFCCF.mg*: an implementation of Belabas' CCFCCF in magma, seems to agree with output of pari-gp equivalent so hopefully few bugs.
- *prototype.mg*: an initial attempt at implementing product of ramified primes in magma, probably full of bugs, decided midway that I was as well just doing it directly in C.

###c\_code
*CCFCCF_prp.c*: contains a C implementation of: 
- _CCFCCF_: Belabas' CCFCCF;
- _CCFCCFPRP_: (CCFCCF Product of Ramified Primes) Belabas but adapted for product of ramified primes; 
- _CCFCCFPRP\_distributed_ distributed variant for parallelising computation in CCFCCFPRP.

uses math libry, so compile with e.g. gcc -Wall CCFCCF.c -lm

## Todo
- [ ] wrap CCRCCR\_prp.c into a library rather than a .c
- [ ] comment code carefully
- [ ] make efficiency edits for product of ramified prime counts
- [ ] \(optional) fix prototype.mg so that it at least works
- [ ] \(optional) write CRFCRF version (for real cubics)


[1] K. Belabas, _A fast algorithm to compute cubic fields_, Math. Comp. *66* (1997), no. 219, 1213-1237


