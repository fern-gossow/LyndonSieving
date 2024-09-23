// Functions for constructing Lyndon structures and coloured cyclic compositions

intrinsic RamanujanSum(d::RngIntElt, j::RngIntElt) -> RngIntElt
{Return the Ramanujan sum c_j(d), sum of the jth powers of dth primitve roots of unity}
    requirege d, 1;
    x := ExactQuotient(d,GCD(j,d));
    return MoebiusMu(x)*ExactQuotient(EulerPhi(d),EulerPhi(x));
end intrinsic;

intrinsic IsGaussCongruence(X::UserProgram, n::RngIntElt) -> BoolElt
{Tests where the values of X satisfy the Gauss congruence for n, i.e. Sum(mu(n/d)X(d),d|n)=0 (mod n)}
    require n ge 1: "n must be positive";
    require &and[Type(X(d)) eq RngIntElt : d in Divisors(n)]: "X must return integer values";
    return &+[MoebiusMu(ExactQuotient(n,d))*X(d) : d in Divisors(n)] mod n eq 0;
end intrinsic;

intrinsic IsGaussCongruenceFull(X::UserProgram, n::RngIntElt) -> BoolElt
{Tests where the values of X satisfy the Gauss congruence for n and all values of the Ramanujan sum}
    require n ge 1: "n must be positive";
    require &and[Type(X(d)) eq RngIntElt : d in Divisors(n)]: "X must return integer values";
    return &and[&+[RamanujanSum(ExactQuotient(n,d),j)*X(d) : d in Divisors(n)] mod n eq 0 : j in Divisors(n)];
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

intrinsic ColouredCyclicCompositions(c::UserProgram, n::RngIntElt) -> RngIntElt
{Return the number of coloured cyclic compositions of size n coloured by c}
    require n ge 1: "n must be positive";
    cols := [c(m) : m in [1..n]];
    require &and[Type(col) eq RngIntElt : col in cols]: "c must return integer values";
    x := 0;
    for k in [1..n] do
        for alpha in Compositions(n,k) do
            x +:= alpha[1]*&*[cols[alpha[i]] : i in [1..k]];
        end for;
    end for;
    return x;
end intrinsic;

intrinsic ColouredCyclicCompositions(c::UserProgram, n::RngIntElt, k::RngIntElt) -> RngIntElt
{Return the number of coloured cyclic compositions of size n coloured by c}
    require n ge 1 and k ge 1 and n ge k: "k must be between 1 and n";
    cols := [c(m) : m in [1..n]];
    require &and[Type(col) eq RngIntElt : col in cols]: "c must return integer values";
    x := 0;
    for alpha in Compositions(n,k) do
        x +:= alpha[1]*&*[cols[alpha[i]] : i in [1..k]];
    end for;
    return x;
end intrinsic;

intrinsic IsLyndonStructure(X::UserProgram, max_n::RngIntElt) -> BoolElt, SeqEnum[RngIntElt]
{If X counts a Lyndon structure up to max_n, return true and the colour values. Else, return false}
    require max_n ge 1: "maximum n value must be positive";
    require &and[Type(X(n)) eq RngIntElt : n in [1..max_n]]: "X must return integer values";
    // Create array of irreducible counts, initially zero
    cols := [0 : n in [1..max_n]];
    cols[1] := X(1);
    for n in [2..max_n] do
        // Find all coloured cyclic compositions made from smaller parts, then subtract from X(n) and divide by n
        new_elts := (X(n) - &+[&+[alpha[1]*&*[cols[alpha[i]] : i in [1..#alpha]] : alpha in Compositions(n,k)] : k in [2..n]]);
        if new_elts mod n eq 0 then
            cols[n] := ExactQuotient(new_elts,n);
        else
            return false;
        end if;
    end for;
    return true, cols;
end intrinsic;