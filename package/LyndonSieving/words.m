// Intrinsics for words, rotations, major index and the q-exponential

intrinsic MajorIndex(w::[RngIntElt]) -> RngIntElt
{Return the major index of a word}
    desc := [i : i in [1..#w-1] | w[i] gt w[i+1]];
    if #desc eq 0 then
        return 0;
    else
        return &+desc;
    end if;
end intrinsic;

intrinsic MajorIndex(w::[RngIntElt]) -> RngIntElt
{Return the major index of a word}
    desc := [i : i in [1..#w-1] | w[i] gt w[i+1]];
    if #desc eq 0 then
        return 0;
    else
        return &+desc;
    end if;
end intrinsic;

intrinsic MajorIndexPolynomial(W::SetEnum[SeqEnum[RngIntElt]]) -> RngUPolElt
{Return the major index polynomial of a set of words}
    R<q> := PolynomialRing(Integers());
    if #W eq 0 then
        return Zero(R);
    else
        return &+[q^MajorIndex(w) : w in W];
    end if;
end intrinsic;

intrinsic QExponential(k::RngIntElt, n::RngIntElt) -> RngUPolElt
{Return the major index polynomial for words of length n over 1..k, which specialises to k^n}
    requirege n, 0;
    requirege k, 0;
    return MajorIndexPolynomial(Subsequences({1..k},n));
end intrinsic;