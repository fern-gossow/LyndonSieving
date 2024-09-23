// Functions for constructing Lyndon structures and coloured cyclic compositions

intrinsic RamanujanSum(d::RngIntElt, j::RngIntElt) -> RngIntElt
{Return the Ramanujan sum c_j(d), sum of the jth powers of dth primitve roots of unity}
    requirege d, 1;
    x := ExactQuotient(d,GCD(j,d));
    return MoebiusMu(x)*ExactQuotient(EulerPhi(d),EulerPhi(x));
end intrinsic;

intrinsic IsGaussCongruence(X::[RngIntElt]) -> BoolElt
{Check whether the values of X satisfy the Gauss congruence}
    for n in [1..#X] do
        if not &+[MoebiusMu(ExactQuotient(n,d))*X[d] : d in Divisors(n)] mod n eq 0 then
            return false;
        end if;
    end for;
    // If passes all checks return true
    return true;
end intrinsic;

intrinsic IsGaussCongruence(X::[RngIntElt], j::RngInt) -> BoolElt
{Check whether the values of X satisfy the Ramanujan-sum Gauss congruence at j}
    for n in [1..#X] do
        if not &+[RamanujanSum(ExactQuotient(n,d), j)*X[d] : d in Divisors(n)] mod n eq 0 then
            return false;
        end if;
    end for;
    // If passes all checks return true
    return true;
end intrinsic;

intrinsic LyndonParameters(X::[RngIntElt]) -> SeqEnum[RngIntElt]
{Given the Lyndon structure counts, return the Lyndon parameters}
    require IsGaussCongruence(X): "Values of X must be Gauss-congruent";
    return [ExactQuotient(&+[MoebiusMu(ExactQuotient(n,d))*X[d] : d in Divisors(n)], n) : n in [1..#X]];
end intrinsic;

intrinsic LyndonSizesFromParameters(B::[RngIntElt]) -> SeqEnum[RngIntElt]
{Given a sequence of Lyndon parameters, determine the sizes for the Lyndon structure}
    return [&+[d*B[d] : d in Divisors(n)] : n in [1..#B]];
end intrinsic;

// Given a subset of k-1 bars in n+k-1 objects, return the composition of n
SubsetToComposition := function(set, max)
    if set eq {} then
        return [max];
    end if;
    seq := Sort(Setseq(set));
    assert max ge seq[#seq];
    seq := &cat[[0], seq, [max+1]];
    return [seq[i+1] - seq[i] - 1 : i in [1..#seq-1]];
end function;

intrinsic WeakCompositions(n::RngIntElt, k::RngIntElt) -> SetEnum[SeqEnum[RngIntElt]]
{Return all nonnegative compositions of n with length k}
    require n ge 0: "n must be positive";
    require k ge 0: "k must be nonnegative";
    return {SubsetToComposition(set, n+k-1) : set in Subsets({1..n+k-1},k-1)};
end intrinsic;

intrinsic Compositions(n::RngIntElt, k::RngIntElt) -> SetEnum[SeqEnum[RngIntElt]]
{Return all (positive) compositions of n with length k}
    require n ge 1 and k ge 1 and n ge k: "k must be between 1 and n";
    return {[y + 1 : y in x] : x in WeakCompositions(n-k, k)};
end intrinsic;

intrinsic ColouredCyclicCompositions(C::[RngIntElt], n::RngIntElt, k::RngIntElt) -> RngIntElt
{Return the number of coloured cyclic compositions of n into k parts coloured by C}
    require n le #C: "n must be less than the length of colours provided";
    require 1 le k and k le n: "k must be between 1 and n";
    x := 0;
    for a in Compositions(n,k) do
        x +:= a[1]*&*[C[a[i]] : i in [1..k]];
    end for;
    return x;
end intrinsic;

intrinsic ColouredCyclicCompositions(C::[RngIntElt], n::RngIntElt) -> RngIntElt
{Return the number of coloured cyclic compositions of n coloured by C}
    require n le #C: "n must be less than the length of colours provided";
    return &+[ColouredCyclicCompositions(C,n,k) : k in [1..n]];
end intrinsic;

intrinsic LyndonSizesFromColours(C::[RngIntElt]) -> RngIntElt
{Return the number of coloured cyclic compositions coloured by C}
    // Array for storing the values of X
    X := [0 : c in [1..#C]];
    for n in [1..#C] do
        // Recursive formula for computation
        x := n*C[n];
        for m in [1..n-1] do
            x +:= X[m]*C[n-m];
        end for;
        X[n] := x;
    end for;
    return X;
end intrinsic;

intrinsic LyndonColours(X::[RngIntElt]) -> SeqEnum[RngIntElt]
{Find colouring function so that X which induces the counts of X as a Lyndon structure}
    // Create array of colour counts, initially zero
    C := [0 : n in [1..#X]];
    C[1] := X[1];
    for n in [2..#X] do
        // Find all coloured cyclic compositions made from smaller parts, then subtract from X(n) and divide by n
        new_elts := (X[n] - &+[X[m]*C[n-m] : m in [1..n-1]]);
        require new_elts mod n eq 0: "The values of X must count a coloured cyclic composition";
        C[n] := ExactQuotient(new_elts,n);
    end for;
    return C;
end intrinsic;

intrinsic LyndonPolynomials(X::[RngIntElt]) -> SeqEnum[RngUPolElt]
{Return the reduced polynomials which give a CSP with X}
    R<q> := PolynomialRing(Integers());
    // Array to store the polynomials
    polys := [Zero(R) : x in [1..#X]];
    for n in [1..#X] do
        j_params := [&+[RamanujanSum(ExactQuotient(n,d), j)*X[d] : d in Divisors(n)] : j in [0..n-1]];
        require &and([p mod n eq 0 : p in j_params]): "X must be fully Gauss congruent";
        polys[n] := &+[ExactQuotient(j_params[j],n) * q^(j-1) : j in [1..n]];
    end for;
    return polys;
end intrinsic;