// Construct and test Lyndon families of polynomials

intrinsic QNumber(n :: RngIntElt) -> RngUPolElt
{Return the q-analogue [n]_q=1+q+...+q^(n-1)}
    R<q> := PolynomialRing(Integers());
    if n le 0 then
        return Zero(R);
    end if;
    return &+[q^i : i in [0..n-1]];
end intrinsic;

intrinsic QFactorial(n :: RngIntElt) -> RngUPolElt
{Return the q-factorial [n]_q!=[1]_q...[n]_q}
    R<q> := PolynomialRing(Integers());
    if n le 0 then
        return Zero(R);
    end if;
    qnums := [&+[q^i : i in [0..m-1]] : m in [1..n]];
    return &*qnums;
end intrinsic;

intrinsic QMultinomial(a :: [RngIntElt]) -> RngUPolElt
{Return the q-multinomial [a_1,a_2,...]_q!=[a_1+a_2+...]_q/([a_1]_q[a_2]_q...)}
    R<q> := PolynomialRing(Integers());
    // If composition contains negative values return 0
    if &or[x lt 0 : x in a] then
        return Zero(R);
    elif #a eq 0 or &+a eq 0 then
        return One(R);
    end if;
    qnums := [&+[q^i : i in [0..m-1]] : m in [1..&+a]];
    return ExactQuotient(&*qnums[1..&+a], &*[&*qnums[1..x] : x in a | x gt 0]);
end intrinsic;

intrinsic QBinomial(n :: RngIntElt, k :: RngIntElt) -> RngUPolElt
{Return the q-binomial [(n,k)]_q!=([n]_q...[n-k+1]_q)/([k]_q...[1]_q)}
    return QMultinomial([k,n-k]);
end intrinsic;

intrinsic TestLyndonFamily(f :: UserProgram, n :: RngIntElt) -> BoolElt
{Given a polynomial family f_n (as f(n)=f_n), test the Lyndon family condition by substituting roots of unity}
    requirege n, 1;
    for d in Divisors(n) do
        require Type(f(d)) eq RngUPolElt: "Outputs of f must be integer polynomials";
    end for;
    A<z> := CyclotomicField(n);
    return &and[Evaluate(f(n), z^d) eq Evaluate(f(d), One(A)) : d in Divisors(n)];
end intrinsic;

intrinsic IsLyndonFamily(F::[RngUPolElt]) -> BoolElt
{Given a sequence of polynomialss, test the Lyndon family condition by substituting roots of unity}
    for n in [1..#F] do
        // Substitute powers of the nth root of unity
        A<z> := CyclotomicField(n);
        for d in Divisors(n) do
            if not Evaluate(F[n], z^d) eq Evaluate(F[d], One(A)) then
                return false;
            end if;
        end for;
    end for;
    return true;
end intrinsic;

intrinsic TestLyndonFamily(f :: UserProgram, rk :: RngIntElt, s :: [RngIntElt]) -> BoolElt
{Given a polynomial family indexed by s and given rank, test the Lyndon family condition by substituting roots of unity}
    require rk ge 1: "Rank must be positive";
    require &and[a ge 0 : a in s] : "Elements of s must be nonnegative";
    s_divisors := &meet[Set(Divisors(a)) : a in s | a gt 0];
    require s_divisors subset Divisors(rk) : "Divisors of s must be divisors of rk";
    for d in s_divisors do
        require Type(f([ExactQuotient(a,d) : a in s])) eq RngUPolElt: "Outputs of f must be integer polynomials";
    end for;
    A<z> := CyclotomicField(rk);
    // Evaluate for all divisors of rk, check Lyndon condition
    for d in Divisors(rk) do
        if d in s_divisors then
            if not Evaluate(f(s), z^ExactQuotient(rk,d)) eq Evaluate(f([ExactQuotient(a, d) : a in s]), One(A)) then
                return false;
            end if;
        else
            if not Evaluate(f(s), z^ExactQuotient(rk,d)) eq Zero(A) then
                return false;
            end if;
        end if;
    end for;
    // If all checks are passed, return true
    return true;
end intrinsic;

intrinsic ReducePolynomial(f::RngUPolElt, n::RngIntElt) -> RngUPolElt
{Reduce a polynomial mod q^n, but keep it in the integer ring}
    q := Parent(f).1;
    return &+[&+[Coefficient(f,m*n+j) : m in [0..Round(Degree(f)/n)+1]]*q^j : j in [0..n-1]];
end intrinsic;

intrinsic ReducePolynomialFamily(F::[RngUPolElt]) -> SeqEnum[RngUPolElt]
{Reduce a sequence of polynomials, where F[n] is reduced mod q^n-1}
    return [ReducePolynomial(F[n], n) : n in [1..#F]];
end intrinsic;

intrinsic LyndonSizesFromPolynomials(F::[RngUPolElt]) -> SeqEnum[RngIntElt]
{Return the evaluations of the polynomial family at q = 1}
    return [Evaluate(f,1) : f in F];
end intrinsic;