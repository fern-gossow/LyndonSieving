// Functions for constructing Lyndon structures and coloured cyclic compositions

intrinsic RamanujanSum(d::RngIntElt, j::RngIntElt) -> RngIntElt
{Return the Ramanujan sum c_j(d), sum of the jth powers of dth primitve roots of unity}
    requirege d, 1;
    x := ExactQuotient(d,GCD(j,d));
    return MoebiusMu(x)*ExactQuotient(EulerPhi(d),EulerPhi(x));
end intrinsic;

intrinsic TestGaussCongruence(X::UserProgram, n::RngIntElt) -> BoolElt
{Tests where the values of X satisfy the Gauss congruence for n, i.e. Sum(mu(n/d)X(d),d|n)=0 (mod n)}
    require n ge 1: "n must be positive";
    require &and[Type(X(d)) eq RngIntElt : d in Divisors(n)]: "X must return integer values";
    return &+[MoebiusMu(ExactQuotient(n,d))*X(d) : d in Divisors(n)] mod n eq 0;
end intrinsic;

intrinsic TestGaussCongruenceFull(X::UserProgram, n::RngIntElt) -> BoolElt
{Tests where the values of X satisfy the Gauss congruence for n and all values of the Ramanujan sum}
    require n ge 1: "n must be positive";
    require &and[Type(X(d)) eq RngIntElt : d in Divisors(n)]: "X must return integer values";
    return &and[&+[RamanujanSum(ExactQuotient(n,d),j)*X(d) : d in Divisors(n)] mod n eq 0 : j in Divisors(n)];
end intrinsic;