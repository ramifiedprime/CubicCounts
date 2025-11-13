# Cubic Counts
Repository for some code to count cubic fields ordered by product of ramified primes.  Essentially just implementing the algorithm in Belabas' article [1], which counts by discriminant, and modifying in appropriate places to get the required count.  Only interested in complex cubics for now, since that's where my practical interest comes in 

##Files
- *CCFCCF.mg*: an implementation of Belabas' CCFCCF in magma, seems to agree with output of pari-gp equivalent so hopefully few bugs.
- *c_code/CCFCCF_prp.c*: contains a C implementation of Belabas' CCFCCF, as well as functions CCFCCFPRP (_CCFCCF Product of Ramified Primes_) and CCFCCFPRP\_distributed (distributed variant for parallelising computation).

## Todo
- [ ] wrap CCRCCR\_prp.c into a library rather than a .c
- [ ] comment code carefully
- [ ] make efficiency edits for product of ramified prime counts
- [ ] \(optional) write CRFCRF version (for real cubics)


[1] K. Belabas, _A fast algorithm to compute cubic fields_, Math. Comp. *66* (1997), no. 219, 1213-1237


