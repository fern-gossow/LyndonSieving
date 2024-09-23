// Functions for constructing Lyndon structures and coloured cyclic compositions

intrinsic RamanujanSum(d::RngIntElt, j::RngIntElt) -> RngIntElt
{Return the Ramanujan sum c_j(d), sum of the jth powers of dth primitve roots of unity}
    requirege d, 1;
    j := j mod d;
    if j eq 0 then
        j := d;
    end if;
    x := ExactQuotient(d,GCD(j,d));
    return MoebiusMu(x)*ExactQuotient(EulerPhi(d),EulerPhi(x));
end intrinsic;

