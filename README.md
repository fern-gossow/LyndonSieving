# LyndonSieving
## A MAGMA package for the Lyndon-like cyclic sieving phenomenon

The purpose of this package is to assist in computing examples of the _Lyndon-like cyclic sieving phenomenon_, first introduced in a [paper](https://arxiv.org/abs/1903.01327) by Alexandersson, Linusson and Potka in 2019, and further investigated by myself in 2024 [here](https://arxiv.org/search/math?searchtype=author&query=Gossow,+F).

The MAGMA package and spec file is contained in the folder *LyndonSieving*, and some Python code for scraping the OEIS and finding examples can be found in *oeis-scrape*.

### Mathematical background

We begin by defining Lyndon-like cyclic sieving for families indexed by the positive integers, and explain their connection to Gauss congruence.

Let $X=\{X_n\}$ be a family of finite sets each equipped with a group action of $C_n$. Suppose moreover that for every $d\mid n$, the number of elements in $X_n$ which are fixed by the action of the subgroup $C_d\leq C_n$ for every $d\mid n$ is $|X_{n/d}|$. This is a fairly strict condition, but it naturally arises in the study of necklaces and words under rotation.

For each $n$, suppose we can find a integer polynomial $f_n(q)$ such that $f_n(\omega_d)=|X_{n/d}|$ for every $d\mid n$, where $\omega_d$ is a primitive $d^\text{th}$ root of unity. Then the triple $(X_n,C_n,f_n)$ is a _cyclic sieving triple_, since evaluations at roots of unity count fixed points under the associated subgroup action. Given such a family of polynomials $\{f_n\}$, we say the family of triples $\{(X_n,C_n,f_n)\}$ exhibit the _Lyndon-like cyclic sieving phenomenon_.

If this is the case, we have $f_n(1)=|X_n|$ and $f_n(\omega_d)=f_{n/d}(1)$ for every $d\mid n$. This second condition is equivalent to $q$-Gauss congruence on $f_n(q)$, i.e. that the polynomials satisfy

$$\sum_{d\mid n}\mu(d)f_{n/d}(q^d)\equiv 0\pmod{[n]_q},$$

where $\mu$ is the Möbius function and $[n]_q:=1+q+\cdots+q^{n-1}$. By substituting $q\mapsto 1$, this in turn implies that the values $a_n:=|X_n|$ enumerating the families satisfy the Gauss congruence:

$$\sum_{d\mid n}\mu(d)a_{n/d}\equiv 0\pmod{n}.$$

By Möbius inversion, we can find a family of integers $\{b_n\}$ such that

$$a_n=\sum_{d\mid n}b_dd$$

for every $n\geq 1$. We call these the _Lyndon parameters_ of $X$. Note that $b_n$ counts the number of orbits of (maximal) rank in $X_n$ under the action of $C_n$. If we set

$$f_n(q)=\sum_{d\mid n}b_d[d]_{q^{n/d}},$$

then it turns out that $(X_n,C_n,g_n)$ is a cyclic sieving triple if and only if $g_n(q)\equiv f_n(q)\pmod{[n]_q}$, so we have a canonical way of writing the polynomials which make up our triples.

Our primary example of a Lyndon structure is _festoons_, which generalise combinatorial necklaces. Given a set of integer values $\{c_n\}$, colour the bead of length $n$ in $c_n$ ways. If $X_n$ is the set of ways to place beads around a circle whose lengths sum to $n$, then $X$ is a Lyndon structure with counts

$$a_n=\sum_{k\geq 1}\sum_{n_1+\cdots+n_k=n}n_1c_{n_1}\cdots c_{n_k}.$$

If $c_n$ takes negative values, then $a_n$ represents a signed count of these festoons. We prove in the paper that for _any_ sequence $a_n$ satisfying the Gauss congruence, we can find values $\{c_n\}$ satisfying the formula, and we call these the _colour values_.

### LyndonSieving package

The MAGMA package is split up into three files:
- `structure.m` provides intrinsics for translating between the families $\{a_n\}$, $\{b_n\}$, $\{c_n\}$ and $\{f_n\}$, as well as check whether a given sequence of integers or integer polynomials satisfies the Gauss congruence. These families should always be represented as a sequence of values, where the first element corresponds to $n=1$.
- `polynomials.m` assists in the computation of polynomials which are commonly related to the $q$-Gauss congruence, such as the $q$-binomial coefficients. Given a sequence of polynomials, the function `ReducePolynomialSequence` maps
$$[f_1,f_2,\dots,f_n]\mapsto [f_1\pmod{[1]_q},f_2\pmod{[2]_q},\dots,f_n\pmod{[n]_q}]$$
which is useful for quickly checking whether two polynomial sequences correspond to the same Lyndon structure.
- `words.m` calculates the major index polynomial for a set of words (represented as integer sequences), as well as the $q$-analogue of the exponential function defined in the paper, giving another common class of $q$-Gauss congruences.

### Installation
If you already have a user-package folder for MAGMA, copy the contents of `LyndonSieving/package` into this folder, and append the contents of spec to your spec file. Otherwise, copy the contents of the package folder into a new folder in your MAGMA directory. You will have to tell MAGMA to look for the spec file by creating/setting the system environment variable `MAMGA_USER_SPEC` to the location of this spec file. More information is available on the [MAGMA](http://magma.maths.usyd.edu.au/magma/) page.

### Example

(1) It is a fact of [Zarelua](https://link.springer.com/article/10.1134/S008154380804007X) that for any integer matrix $M$, the sequence $a_n=\mathrm{tr}(A^n)$ satisfies the Gauss congruence. We claim furthermore that the colour values $c_n$ are the negatives of the coefficients in the characteristic polynomials for $M$. We check this assertion using our package.

    > M := Matrix([[1,5,2],[0,-1,1],[3,0,1]]);
    > A := [Trace(M^n) : n in [1..8]];
    > IsGaussCongruence(A);
    true
    > LyndonColours(A);
    [ 1, 1, 14, 0, 0, 0, 0, 0 ]
    > R<q> := PolynomialRing(Integers());
    > CharacteristicPolynomial(M);
    q^3 - q^2 - q - 14

(2) We find the number of binary Lyndon words on $n$ letters, and confirm that the major index polynomial is congruent to the major index polynomial in a particular case.

    > W := [Subsequences({0,1},n) : n in [1..12]];
    > A := [#w : w in W];
    > IsGaussCongruence(A);
    true
    > LyndonParameters(A);
    [ 2, 1, 2, 3, 6, 9, 18, 30 ]
    > F := [MajorIndexPolynomial(w) : w in W];
    > ReducePolynomial(F[4], 4);
    3*q^3 + 4*q^2 + 3*q + 6
    > LyndonPolynomials(A)[4];
    3*q^3 + 4*q^2 + 3*q + 6

(3) Given a set of Lyndon parameters, we determine the enumeration of the Lyndon structure.

    > B := [(-1)^n : n in [1..8]];
    > LyndonSizesFromParameters(B);
    [ -1, 1, -4, 5, -6, 4, -8, 13 ]

(4) Finally, we see that the Catalan numbers do not satisfy the Gauss congruence.

    > A := [Catalan(n) : n in [1..6]];
    > IsGaussCongruence(A)
    false

Trying to determine the Lyndon parameters, colours or polynomials in this case will raise an error.

### Future Work

Generalise to index sets other than the positive integers.
