freeze;
// Construct q-analogues of binomial coefficients, and reduce polynomials mod [n]_q by summing over coefficients

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

intrinsic ReducePolynomial(f::RngUPolElt, n::RngIntElt) -> RngUPolElt
{Reduce a polynomial mod q^n, but keep it in the integer ring}
    q := Parent(f).1;
    return &+[&+[Coefficient(f,m*n+j) : m in [0..Round(Degree(f)/n)+1]]*q^j : j in [0..n-1]];
end intrinsic;

intrinsic ReducePolynomialSequence(F::[RngUPolElt]) -> SeqEnum[RngUPolElt]
{Reduce a sequence of polynomials, where F[n] is reduced mod q^n-1}
    return [ReducePolynomial(F[n], n) : n in [1..#F]];
end intrinsic;