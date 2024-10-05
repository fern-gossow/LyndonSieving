// Functions for constructing Lyndon structures and coloured cyclic compositions

intrinsic RamanujanSum(d::RngIntElt, j::RngIntElt) -> RngIntElt
{Return the Ramanujan sum c_j(d), sum of the jth powers of dth primitve roots of unity}
    requirege d, 1;
    x := ExactQuotient(d,GCD(j,d));
    return MoebiusMu(x)*ExactQuotient(EulerPhi(d),EulerPhi(x));
end intrinsic;

intrinsic IsGaussCongruence(A::[RngIntElt]) -> BoolElt
{Check whether an integer sequence satisfies the Gauss congruence}
    for n in [1..#A] do
        if not &+[MoebiusMu(ExactQuotient(n,d))*A[d] : d in Divisors(n)] mod n eq 0 then
            return false;
        end if;
    end for;
    // If passes all checks return true
    return true;
end intrinsic;

intrinsic IsGaussCongruence(F::[RngUPolElt]) -> BoolElt
{Check whether an intege polynoimal sequence satisfies the q-Gauss congruence}
    for n in [1..#F] do
        // Substitute powers of the nth root of unity
        A<z> := CyclotomicField(n);
        for d in Divisors(n) do
            if not Evaluate(F[n], z^d) eq Evaluate(F[d], One(A)) then
                return false;
            end if;
        end for;
    end for;
    // If passes all checks return true
    return true;
end intrinsic;

intrinsic IsGaussCongruence(X::[RngIntElt], j::RngInt) -> BoolElt
{Check whether the sequence satisfies the Ramanujan-sum-adjusted Gauss congruence at j}
    for n in [1..#X] do
        if not &+[RamanujanSum(ExactQuotient(n,d), j)*X[d] : d in Divisors(n)] mod n eq 0 then
            return false;
        end if;
    end for;
    // If passes all checks return true
    return true;
end intrinsic;

intrinsic LyndonParameters(A::[RngIntElt]) -> SeqEnum[RngIntElt]
{Given a sequence satisfying Gauss congruence, return the Lyndon parameters B}
    require IsGaussCongruence(A): "Values of A must satisfy Gauss congruence";
    return [ExactQuotient(&+[MoebiusMu(ExactQuotient(n,d))*A[d] : d in Divisors(n)], n) : n in [1..#X]];
end intrinsic;

intrinsic LyndonSizesFromParameters(B::[RngIntElt]) -> SeqEnum[RngIntElt]
{Given a sequence of Lyndon parameters, return the sizes A}
    return [&+[d*B[d] : d in Divisors(n)] : n in [1..#B]];
end intrinsic;

intrinsic LyndonColours(A::[RngIntElt]) -> SeqEnum[RngIntElt]
{Given a sequence satisfying Gauss congruence, return the Lyndon colours C}
    require IsGaussCongruence(A): "Values of A must satisfy Gauss congruence";
    // Create array of colour counts
    C := [0 : n in [1..#A]];
    C[1] := A[1];
    for n in [2..#A] do
        C[n] := ExactQuotient(A[n] - &+[A[m]*C[n-m] : m in [1..n-1]],n);
    end for;
    return C;
end intrinsic;

intrinsic LyndonSizesFromColours(C::[RngIntElt]) -> RngIntElt
{Return the (signed) count of festoons of each length coloured by C}
    // Create array of A values
    A := [0 : c in [1..#C]];
    for n in [1..#C] do
        A[n] = n*C[n] + &+[A[m]*C[n-m] : m in [1..n-1]]
    end for;
    return X;
end intrinsic;

intrinsic LyndonPolynomials(A::[RngIntElt]) -> SeqEnum[RngUPolElt]
{Given a sequence satisfying Gauss congruence, return the associated sequence of reduced polynomials satisfying q-Gauss congruence}
    require IsGaussCongruence(A): "Values of A must satisfy Gauss congruence";
    B := LyndonParameters(A);
    return &+[B[d]*Evaluate(QNumber(d),q^ExactQuotient(n,d)) : d in Divisors(n)];
end intrinsic;

intrinsic LyndonSizesFromPolynomials(F::[RngUPolElt])
{Given a sequence of polynomials satisfying q-Gauss congruence, return their evaluations at q=1}
    require IsGaussCongruence(F): "Values of F must satisfy q-Gauss congruence";
    return [Evaluate(f,1) : f in F];
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