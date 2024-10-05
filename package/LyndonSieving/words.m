// Intrinsics for words, rotations, major index and the q-exponential

intrinsic MajorIndex(w::[RngIntElt]) -> RngIntElt
{Returns the major index of a word}
    desc := [i : i in [1..#w-1] | w[i] gt w[i+1]];
    if #desc eq 0 then
        return 0;
    else
        return &+desc;
    end if;
end intrinsic;

intrinsic MajorIndex(w::[RngIntElt]) -> RngIntElt
{Returns the major index of a word}
    desc := [i : i in [1..#w-1] | w[i] gt w[i+1]];
    if #desc eq 0 then
        return 0;
    else
        return &+desc;
    end if;
end intrinsic;

intrinsic MajorIndexPolynomial(W::SetEnum[SeqEnum[RngIntElt]]) -> RngUPolElt
{Returns the major index polynomial of a set of words}
    R<q> := PolynomialRing(Integers());
    if #W eq 0 then
        return Zero(R);
    else
        return &+[q^MajorIndex(w) : w in W];
    end if;
end intrinsic;

intrinsic QExponential(k::RngIntElt, n::RngIntElt) -> RngUPolElt
{Returns the major index polynomial for words of length n over 1..k, which specialises to k^n}
    requirege n, 0;
    requirege k, 0;
    return MajorIndexPolynomial(Subsequences({1..k},n));
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
{Returns all nonnegative compositions of n with length k}
    require n ge 0: "n must be positive";
    require k ge 0: "k must be nonnegative";
    return {SubsetToComposition(set, n+k-1) : set in Subsets({1..n+k-1},k-1)};
end intrinsic;

intrinsic Compositions(n::RngIntElt, k::RngIntElt) -> SetEnum[SeqEnum[RngIntElt]]
{Returns all (positive) compositions of n with length k}
    require n ge 1 and k ge 1 and n ge k: "k must be between 1 and n";
    return {[y + 1 : y in x] : x in WeakCompositions(n-k, k)};
end intrinsic;