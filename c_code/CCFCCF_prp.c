#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int gcd(int a, int b){//TODO
    int r;
    while(b!=0){
        r=a%b;
        a=b;
        b=r;
    }
    return abs(a);
}

_Bool member(int x, int* list){
    int n = list[0];
    int lo = 1, hi = n;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (list[mid] == x) return 1;
        if (list[mid] < x) lo = mid + 1;
        else hi = mid - 1;
    }
    return 0;
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
    float num=1;
    while(i!=r[1]){
        num = X/(p[i]*log(p[i]))*(1+(1.0/(2*log(p[i]))));
        if(num <= tM){
            r[1]=i;
        }else{
            r[0]=i;
        }
        i = (r[0]+r[1]+1)/2;
    }
    if(num > tM){
        printf("Error: Memory constraint too small for bound.\n");
        exit(1);
    }//At this stage we've completed the first pass
    if(p[i]<=53){
        float J=0;
        float JJ=(((float)tM)/X - (1/log(59))*(1+(1/(2*log(59)))));
        i = 16;
        while(J<=JJ && i>3){
            J+=1.0/(p[i]*p[i]);
            i--;
        }
    }
    return i;
}

int* get_list(int* pp, int index, int X){
    int k,r,x,i,num=0;
    _Bool* blist = (_Bool*)calloc(X+2, sizeof(_Bool));
    for(k=4; k<=X/6; k++){
        for(r=-1; r<=1; r+=2){
            x=6*k+r;
            // printf("x=%d\n",x);
            for(i=index; i<=pp[0] && pp[i]<=x; i++){
                if(x%pp[i] == 0){
                    blist[x] = 1;
                    num++;
                    break;
                }
            }
        }
    }
    int* list = (int*)calloc(num+1, sizeof(int));
    list[0]=num;
    r = 1;
    for(k=0; k<=X+1; k++){
        if(blist[k]){
            list[r]=k;
            r++;
        }
    }
    free(blist);
    return list;
}

int* get_sqfull(int* pp, int X){
    int i,n;
    int* sqfull = (int*)malloc((sqrt(3*X)+1)*sizeof(int));
    for(n=2; n<=sqrt(3*X); n++){
        sqfull[n]=0;
        for(i=3; i<=pp[0]; i++){
            if(pp[i]<=n){
                break;
            }
            if(n%pp[i]==0){
                sqfull[n]=1;
            }
        }
    }
    return sqfull;
}



_Bool test(int a, int b, int c, int d, int D, int f, int* sqfull, int sqrtX, int* list, int* pp, int index){
    int a9,d9,t,i;
    float sqrt3X=floor(1.7320508075688772*sqrtX);
    if ((D%27)==0 && (f%3!=0)){return 0;} // Not in V3 and not (1^3)
    else if(f%3==0){// checking if in U3 given that it's 1^3
        a9 = a%9;
        d9 = d%9;
        if(a9==0 || d9==0){return 0;}
        else if(((a9%3)==0)&&((d9%3)==0)){return 0;}
        else if((((a9-d9)%3)==0) && ((a9-b+c-d9)%9)==0){return 0;}
        else if((((a9+d9)%3)==0) && ((a9+b+c+d9)%9)==0){return 0;}
    }

    if(sqfull[f]){return 0;}
    t = abs(D/(f*f));
    t = (t/gcd(t,72));
    if(gcd(t,f)!=1){return 0;}
    if (t<=sqrt3X && sqfull[t]){return 0;} 
    else if(member(t,list)){return 0;}
    else{
        for(i=2;i<=index-1;i++){
            if(t%pp[i]==0){return 0;}
        }
    }
    return 1;
}

_Bool is_complex_field(int a, int b, int c, int d, int P, int Q, int R, int D, int f,int* sqfull, int sqrtX, int* list, int* pp, int index){
        if((D <= 0) || (D%16==0) || (D%16==4 && (P%2==1 || R%2==1))){
            return 0;
        }
        return test(a,b,c,d,D,f,sqfull,sqrtX,list,pp,index);
}

float U(int a, int b){
    if(a>=2.0*b/3){return ((float)b*b)/(3.0*a);}
    else{return b-3*((float)a)/4.0;}
}


int main(void){
    FILE *fptr=fopen("output.dat","w");
    int i,r=6;
    int a,b,c,d,D,f,P,Q,R;
    int B=pow(2,r);
    int X=B*B;
    // make p
    int* p = primes_up_to(X,B);
    // make pp
    int pp[p[0]+1];
    pp[0]=p[0];
    for(i=1; i<=p[0]; i++){
        pp[i]=p[i]*p[i];
    }
    // get index, list, sqfull
    int index = init_index_find(p, X, X);
    int* sqfull = get_sqfull(pp,X);
    int* list = get_list(pp, index, X);
    // do looping
    float A_bd=pow((16.0*X)/27, 1.0/4);
    int lbd;
    int quadbd;
    // b=0 looping
    for(a=1;a<=A_bd;a++){
        for(c=1; c<=pow(((float)X)/(4*a), 1.0/3); c++){
            lbd=0;
            if(c<=a){lbd = (int)floor(sqrt(a*(a-c)));}
            for(d=lbd+1; d<=(a+c)-1; d++){ //Lemma 4.2 (12) for LB, (13) for UB 
                P= -3*a*c;
                Q= -9*a*d;
                R= c*c;
                D= Q*Q-4*P*R;
                f=gcd(P,gcd(Q,R));
                // if(D>3*B*f){break;}
                // printf("[%d,%d,%d,%d]\n",a,0,c,d);
                if(is_complex_field(a,0,c,d,P,Q,R,D,f,sqfull,B,list,pp,index)){
                    fprintf(fptr,"[%d,%d,%d,%d]\n",a,0,c,d);
                }
                // printf("[%d,%d,%d,%d]\n",a,0,c,d);
            }
        }
    }
    for(a=1;a<=A_bd;a++){
        for(b=1; b<=(3.0*(float)a)/2 + sqrt(sqrt(((float)X)/3) - (3.0*a*a)/4); b++){    
            for(c=(1-b); c<=U(a,b)+pow(((float)X)/(4.0*a), 1.0/3);c++){
                lbd =(int)floor(-(((float)a-b)*(a-b+c))/a);
                quadbd = a*(a-c);
                for(d=lbd+1; d<(((float)a+b)*(a+b+c))/a; d++){
                    if(d*(d-b) < quadbd){continue;}//Improvement poss here, inc function
                    P= -3*a*c;
                    Q= -9*a*d;
                    R= c*c;
                    D= Q*Q-4*P*R;
                    f=gcd(P,gcd(Q,R));
                    // if(D>3*X*f){break;}
                    if(is_complex_field(a,b,c,d,P,Q,R,D,f,sqfull,B,list,pp,index)){
                        fprintf(fptr,"[%d,%d,%d,%d]\n",a,b,c,d);
                    }
                }
            }
        }
    }
    // printf("%d\n", test(14, 1, 0, -7, 777924, 882, sqfull, sqrtX, list, pp, index));
    // for(int z=25; z<=100;z++){
    //     printf("sqfull[%d]=%d\n", z, sqfull[z]);//FOUND IT
    // }
    fclose(fptr);
    free(p);free(sqfull);free(list);
    return 0;
}

// failure: [5,0,3452,862]


// int write_data(int a, int b, int c, int d, int D, int f){
//     FILE* fptr;
//     fptr = fopen("test.dat", "w");
//     fprintf(fptr, "%d,%d,%d,%d,%d,%d\n", a,b,c,d,D,f);
//     fclose(fptr);
//     return 0;
// }
