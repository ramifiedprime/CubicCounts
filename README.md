# Cubic Counts
Repository for some code to count cubic fields ordered by product of ramified primes.  Essentially just implementing the algorithm in Belabas' article [1], which counts by discriminant, and modifying in appropriate places to get the required count.  Only interested in complex cubics for now, since that's where my practical interest is.

**WARNING:** This is all still very much in development, and output is not guaranteed to be correct.

## Usage
The main useful thing in this repository is the code written in C.  To compile the code (do this once per pull/clone) from the main directory, just run the command below.
```bash
make
```
This will produce the packaged executable `CCFCCFPRP.exe`.  The syntax for using this executable is then the below.
```bash
./CCFCCFPRP.exe [-b] <B_or_exponent> <outfile>
./CCFCCFPRP.exe [-b] <B_or_exponent> <n1> <n2> <outfile>
```
if the `-b` flag is set then one should `<B_or_exponent>` to be the maximum number of bits for the product of ramified primes, else this should just be the bound for the product of ramified primes.  The second line is for parallelisation: it speficies that we only check the cubic forms with leading coefficient $a\equiv n_1\mod n_2$.

## Files
- *runme.sh*: compilation script to compile the CCFCCFPRP user program in c\_code/main.c. 

### c\_code
This is the main directory of useful code.
- *CCFCCF.c*: contains a C implementation of: 
    - _CCFCCF_: Belabas' CCFCCF.
- *CCFCCFPRP.c*: contains a C implementation of:
    - _CCFCCFPRP_: (CCFCCF Product of Ramified Primes) Belabas but adapted for product of ramified primes. 
    - _CCFCCFPRP\_distributed_: distributed variant for parallelising computation in CCFCCFPRP, distributes on the leading coefficient _a_.
- *main.c*: A main runner program which packages up CCFCCFPRP and its distributed equivalent for a user-friendly executable.
### magma\_code
**WARNING:  The magma code is very incomplete and untested, use at significant risk.**
- *CCFCCF.mg*: an implementation of Belabas' CCFCCF in magma, seems to agree with output of pari-gp equivalent so hopefully few bugs.
- *prototype.mg*: an initial attempt at implementing product of ramified primes in magma, probably full of bugs, decided midway that I was as well just doing it directly in C.

## Todo
- [x] wrap CCFCCF.c into a library rather than a .c
- [x] wrap CCFCCFPRP.c into a library rather than a .c
- [x] comment code carefully
- [x] make efficiency edits for product of ramified prime countsi
- [ ] Clean up code to remove tX and X etc, so things are clearer
- [ ] \(optional) fix prototype.mg so that it at least works
- [ ] \(optional) write CRFCRF version (for real cubics)


[1] K. Belabas, _A fast algorithm to compute cubic fields_, Math. Comp. *66* (1997), no. 219, 1213-1237


