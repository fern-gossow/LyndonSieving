// Construct and test Lyndon families of polynomials

intrinsic QNumber(n :: RngIntElt) -> RngUPolElt
{Return the q-analogue [n]_q=1+q+...+q^(n-1)}
    R<q> := PolynomialRing(Integers());
    if n le 0 then
        return Zero(R);
    end if;
    return &+[q^i : i in [0..n-1]];
end intrinsic;

intrinsic QFactorial := function(n :: RngIntElt) -> RngUPolElt
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
    qmult := &*qnums[1..&+a]/&*[&*qnums[1..x] : x in a | x gt 0];
    // Ensure the divison was completed
    assert Denominator(qmult) eq One(R);
    return Numerator(qmult);
end intrinsic;

intrinsic QBinomial(n :: RngIntElt, k :: RngIntElt) -> RngUPolElt
{Return the q-binomial [(n,k)]_q!=([n]_q...[n-k+1]_q)/([k]_q...[1]_q)}
    return QMultinomial([k,n-k]);
end intrinsic;

