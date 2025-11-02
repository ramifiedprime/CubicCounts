function init(X,M)
    p:=PrimesUpTo(Floor(Sqrt(X)));
    pp:=[x^2 : x in p];

    function inittest1(q,MMM)
        i:=Ceiling(#q/2);
        if #q eq 1 then
            if (X/Log(q[i]))*(1+(1/2*Log(q[i]))) gt MMM then 
                print("Error: memory constraint too small for bound");
                assert false;
            else
                return 1;
            end if;
        elif (X/Log(q[i]))*(1+(1/2*Log(q[i]))) le MMM then
            return inittest1(q[1..i],MMM);
        else 
            return i+inittest1(q[i+1..#q],MMM);
        end if;
    end function;

    function inittest2(p, MMM)
        i:=16;
        J:=p[i]^(-2);
        B:=MMM/X - (1/Log(59))*(1+1/(2*Log(59)));
        while (J le B) and (i gt 2) do
            J:=J+p[i-1]^2;
            i:=i-1;
        end while;
        return i+1;
    end function;

    index:=inittest1(p,3*M);
    if p[index] le 53 then
        index:=inittest2(p,3*M);
    end if;

    list:=[];
    for k in [4..Floor(X/6)] do
        for r in [-1,1] do
            x:=6*k+r;
            for qq in pp[index..#pp] do
                if (x mod qq) eq 0 then 
                    Append(~list, x);
                    break;
                elif qq gt x then break;
                end if;
            end for;
        end for;
    end for;

    sqfull:=[false:n in [1..Floor(Sqrt(3*X))]];
    for n in [1..#sqfull] do
        for qq in pp[3..#pp] do
            if (n mod qq) eq 0 then
                sqfull[n] := true;
            elif qq gt x then break;
            end if;
        end for;
    end for;

    return p,pp,index,list,sqfull;
end function;

X:=2^(20); M:=2^(18); stX:=Floor(Sqrt(3*X));
// p,pp,index,list,sqfull := init(2^(20),2^(18));
printf("\nInitialising...");
p,pp,index,list,sqfull := init(X,M);
printf("\nInitialisation complete.");

// Assumes input form is in U2
function test(f,a,b,c,d,D)
    if ((D mod 27) eq 0) and (f mod 3 ne 0) then return false; // Not in V3 and not (1^3)
    elif f mod 3 eq 0 then // checking if in U3 given that it's 1^3
        a9 := a mod 9;
        d9 := d mod 9;
        if a9 eq 0 or d9 eq 0 then return false;
        elif a9 mod 3 eq 0 and d9 mod 3 eq 0 then return false;
        elif ((a9-d9) mod 3) eq 0 and ((a9-b+c-d9) mod 9) eq 0 then return false; // something is wrong with this cond
        elif ((a9+d9) mod 3) eq 0 and ((a9+b+c+d9) mod 9) eq 0 then return false;
        end if;
    end if;
    if sqfull[f] then return false;
    end if;
    t := Integers()!Abs(D/(f^2));
    t := Integers()!(t/GCD(t,72));
    if GCD(t,f) ne 1 then return false; end if;
    if t le stX and sqfull[t] then return false; 
    elif t in list then return false;
    else
        for qq in pp[2..index-1] do
            if (t mod qq) eq 0 then 
                return false;
            end if;
        end for;
    end if;
    return true;
end function;

function is_complex_field(a,b,c,d,P,Q,R,D)
    if (D le 0) or (D mod 16 eq 0) or (D mod 16 eq 4 and (IsOdd(P) or IsOdd(R))) then
        return false; // second cdn disallows [1,0,3,2]
    end if;
    f:=GCD([P,Q,R]);
    return test(f,a,b,c,d,D);
end function;

function U(a,b)
    if a ge 2*b/3 then return b^2/(3*a);
    else return b-3*a/4;
    end if;
end function;



fields:=[];
printf("\nRunning loops with b=0...");
for a in [1..Floor((16*X/27)^(1/4))],
    c in [1..Floor((X/(4*a))^(1/3))] do
    lbd:=0;
    if c le a then
        lbd:=Floor(Sqrt(a*(a-c)));
    end if;
    for d in [lbd+1..(a+c)-1] //Lemma 4.2 (12) for LB, (13) for UB 
        do 
        P:=-3*a*c;
        Q:= -9*a*d;
        R:= c^2;
        D:=Q^2-4*P*R;
        if D gt 3*X then continue; end if;
        if is_complex_field(a, 0, c, d, P, Q, R, D) then
            Append(~fields, [a,0,c,d]);
        end if;
    end for;
end for;
printf("\nLoops complete.");

printf("\nRunning loops with b>0...");
for a in [1..Floor((16*X/27)^(1/4))],
    b in [1..Floor((3*a/2) + Sqrt((X/3)^(1/2)-(3*a^2/4)))],
    c in [(1-b)..Floor(U(a,b)+(X/(4*a))^(1/3))],
    d in [Floor((-(a-b)^2-a*c+b*c)/a)+1..Ceiling(((a+b)*(a+b+c))/a)-1] //Lemma 4.2 (13) for both 
    do
    // if d*(d-b) le a*(a-c) then continue; end if;
    // if d^2-a^2+a*c-d*b le 0 then printf("\n(12) violated");
    // elif -(a-b)^2-a*c ge a*d-b*c then print("\n(13) LB violated");
    // elif a*d-b*c ge (a+b)^2+a*c  then print("\n(13) UB violated");
    // elif Q^2-4*P*R le 0 then printf("\nDisc nonneg!");
    // end if; 
    if d*(d-b)-a*(a-c) le 0 then continue; end if;
    P:=b^2-3*a*c;
    Q:= b*c-9*a*d;
    R:= c^2-3*b*d;
    D:=Q^2-4*P*R;
    if D gt 3*X then continue; end if;
    if is_complex_field(a, b, c, d, P, Q, R, D) then
        Append(~fields, [a,b,c,d]);
    end if;
end for;

//Still not finding everything, must be a type.