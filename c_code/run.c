#include"CCFCCFPRP.h"
#include<math.h>
#include<stdlib.h>

int main(int argc, char** argv){
    FILE *fptr=fopen(argv[2], "w");
    CCFCCFPRP(pow(2,atol(argv[1])), fptr, 0);
    // CCFCCFPRP_distributed(pow(2,atol(argv[1])), atol(argv[2]), atol(argv[3]), fptr, 0);
    fclose(fptr);
    return 0;
}
