# Manual

## Finite differences

Consider the analytical function ``f`` tabulated in *forward order*  
(growing index) at ``n`` positions on a *uniform grid*.

**Forward difference notation**

The *forward translation* from position ``n-1`` to position ``n`` on 
the grid is expressed by the relation

```math
f[n] = (1 + Δ) f[n-1] ,
```

where ``Δ`` is the forward difference operator. By formal inversion  
we find

```math
f[n-1]=(1+Δ)^{-1}f[n]=(1-Δ+Δ^2-Δ^3+⋯)f[n],
```

where ``Δ^k`` is the  ``k^{th}``-*order forward difference* defined as
a *weighted sum* over the function values ``f[n:n+k]`` (involving
``k+1`` points appearing in regular=forward order),

```math
\begin{aligned}
Δ^k f[n] &= f[n+k] + c_1^kf[n+k-1] + ⋯  + f[n] \\
         &= f[n] + c_{k-1}^kf[n+1] + ⋯  + f[n+k] = \sum_{j=0}^{k} c_{k-j}^k f[n+j],
\end{aligned}
```
where the ``k+1`` coefficients ``c_{k-j}^k=(-1)^k c_j^k``,  with

```math
c_{j}^{k}=(-1)^{j}\binom{k}{j},
```

are the *summation weights* (short: *weights*) defining the summation, with special values 
``c_{0}^{k}≡1`` and ``c_{k}^{k}≡(-1)^{k}``.

In inner product form the result becomes

```math
Δ^k f[n]=\bar{c}^k \cdot f[n:n+k],
```

where ``\bar{c}^k  ≡ [c_k^k,⋯\ c_0^k]``.

**Backward difference notation**

The *backward translation* from position ``n`` to position ``n-1`` on 
the grid is expressed by the relation

```math
f[n] = (1 - ∇) f[n+1],
```

where ``∇`` is the backward difference operator.  By formal inversion  
we obtain in this case

```math
f[n+1]=(1-∇)^{-1}f[n]=(1+∇+∇^2+∇^3+⋯)f[n],
```

where ``∇^k`` is the  ``k^{th}``-*order backward difference* defined as
a *weighted sum* over the function values ``f[n:-1:n-k]`` (involving 
``k+1`` points appearing in reversed order),

```math
∇^k f[n] = f[n] + c_1^kf[n-1] + ⋯ + f[n-k]
= \sum_{j=0}^{k} c_j^kf[n-j],
```

where the ``k+1`` coefficients ``c_j^k`` are the *summation weights* 
(short: *weights*) defining the summation. 

In inner product form the result becomes

```math
∇^k f[n] = c^k \cdot f[n:-1:n-k] =\bar{c}^k \cdot f[n-k:n],
```
where ``c^k  ≡ [c_0^k,⋯\ c_k^k]``.

Coefficients:

[`fdiff_weight(k, j)`](@ref) ``→ c^k ≡ [c_0^k,⋯\ c_k^k]``

```@docs
# isforward(notation)
# isregular(ordering)
fdiff_weight(k::Int, j::Int)
```

## Finite difference expansions

Finite-difference calculus builds on the *finite-difference expansion*.

**Forward difference notation**

In terms of forward differences the generic form of the finite-difference expansion is given by

```math
\sum_{p=0}^{\infty}α_{p}Δ^{p}f[n]
=\sum_{p=0}^{k}α_{p}Δ^{p}f[n]+⋯.
```

Evaluated to order ``k``, the expansion is defined by ``k+1`` *finite-difference expansion coefficients*, 
supplied by the user in the form of the expansion vector, ``α = [α_{0},⋯\ α_{k}]``. This vector contains the
coefficients in the *regular* ordering of *growing index*. It takes some bookkeeping to rewrite the expansion 
as a *weighted sum* over the ``k+1`` *function values* ``f[n:n+k]`` (note the *regular* ordering of
*growing* grid position). Substituting the definition of the forward difference, ``Δ = f[n+1] - f[n]``, 
the finite-difference expression takes the form

```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n]
=\sum_{p=0}^{k}α_{p}\sum_{j=0}^{p}c_{p-j}^{p}f[n+j]
=\sum_{j=0}^{k}\sum_{p=j}^{k}α_{p}c_{p-j}^{p}f[n+j]
=\sum_{j=0}^{k}F_{j}^{k}f[n+j],
```

where the weighted summation is defined by the *weights*
```math
F_{j}^{k}=\sum_{p=j}^{k}α_{p}c_{p-j}^{p}
=\sum_{p=j}^{k}(-1)^{p+j}\binom{p}{j}α_{p},
```

with ``j=0,⋯\ k``. In inner product form, the expansion becomes
```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n]
=\sum_{j=0}^{k}F_{j}^{k}f[n+j]
=F^{k} \cdot f[n:n+k],
```

where ``F^k  ≡ [F_0^k,⋯\ F_k^k]``.

Coefficients:

[`fdiff_expansion_weights(polynom, fwd, reg)`](@ref)
``→ F^k ≡ [F_0^k,⋯\ F_k^k]``,

where `polynom` is the expansion vector ``α ≡ [α_0,⋯\ α_k]``, which has to be supplied by the user 
to define the expansion under consideration.  

Some common cases are:

  interpolation expansion: [`fdiff_interpolation_expansion_polynom(ξ, k, fwd)`](@ref)

differentiation expansion: [`fdiff_differentiation_expansion_polynom(ξ, k, fwd)`](@ref)

 Adams-Bashford expansion: [`fdiff_adams_bashford_expansion_polynom(k, fwd)`](@ref)

  Adams-Moulton expansion: [`fdiff_adams_moulton_expansion_polynom(k, fwd)`](@ref)

**Backward difference notation**

In terms of backward differences the generic form of the finite-difference expansion is given by

```math
\sum_{p=0}^{\infty}β_{p}∇^{p}f[n]=\sum_{p=0}^{k}β_{p}∇^{p}f[n]+⋯.
```

In this case the ``k^{th}``- order *finite-difference expansion* is defined by the (user-supplied)
expansion vector ``β = [β_{0},⋯\ β_{k}]``, containing the expansion coefficients 
in *regular* ordering (*growing* index). The expansion can written as 
a *weighted sum* over the ``k+1`` *function values* ``f[n:-1:n-k]`` (note *reversed* ordering 
of *decreasing* grid position). Substituting the definition of the backward difference, ``∇ = f[n] - f[n-1]``, 
the finite-difference expression takes the form

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n]
=\sum_{p=0}^{k}β_{p}\sum_{j=0}^{p}c_{j}^{p}f[n-j]
=\sum_{j=0}^{k}\sum_{p=j}^{k}β_{p}c_{j}^{p}f[n-j]
=\sum_{j=0}^{k}B_{j}^{k}f[n-j],
```

where the *weights* are given by

```math
B_{j}^{k}=\sum_{p=j}^{k}β_{p}c_{j}^{p}
=\sum_{p=j}^{k}(-1)^{j}\binom{p}{j}β_{p},
```
with ``j=0,⋯\ k``. In inner product form, the expansion becomes

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n]
=\sum_{j=0}^k B_j^k f[n-j]
=\bar{B}^k \cdot f[n-k:n],
```

where ``\bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]`` is the *weights vector* ``B^{k} ≡ [B_0^k,⋯\ B_k^k]`` tabulated 
in reversed order.

The relation between ``B^k`` and ``F^k`` depends on the relation between the expansion polynoms ``α`` and ``β``,
which is not unique (it depends on the expansion under consideration).

Coefficients:

[`fdiff_expansion_weights(polynom, bwd, rev)`](@ref)
`` → \bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]``,

where `polynom` is the expansion vector ``β ≡ [β_0,⋯\ β_k]``, which has to be supplied by the user 
to define the expansion under consideration. 

Some common cases are:

  interpolation expansion: [`fdiff_interpolation_expansion_polynom(ξ, k, bwd)`](@ref)

differentiation expansion: [`fdiff_differentiation_expansion_polynom(ξ, k, bwd)`](@ref)

 Adams-Bashford expansion: [`fdiff_adams_bashford_expansion_polynom(k, bwd)`](@ref)

  Adams-Moulton expansion: [`fdiff_adams_moulton_expansion_polynom(k, bwd)`](@ref)

```@docs
fdiff_expansion(polynom, f, notation=CamiMath.bwd)
fdiff_expansion_weights(polynom, notation=CamiMath.bwd, ordering=CamiMath.rev)
```

## Lagrange interpolation

The *Lagrange polynomial* of degree ``k`` on a uniform grid is the polynomial running 
through ``k+1`` subsequent points on the [`Grid`](@ref). We derive expressions for 
interpolation in both forward- and backward-difference notation. Beware 
that Lagrange interpolation becomes inaccurate if the tabulated function cannot be 
approximated by a polynomial of degree ``k``.

**Forward difference notation**

Starting from the forward translation expression
```math
f[n]=(1+Δ)f[n-1],
```
we obtain by formal operator inversion
```math
f[n-1] = (1 + Δ)^{-1} f[n] ≡ \sum_{p=0}^{\infty}(-1)^p Δ^p f[n],
```
```math
f[n-2] = (1 + Δ)^{-2} f[n] ≡ \sum_{p=0}^{\infty}(-1)^p pΔ^p f[n],
```
```math
\vdots
```
where ``p`` is called the order of the expansion and ``n`` is the index of the reference
position. For interpolation position ``n-σ`` (where ``σ`` may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n-σ] = (1 + Δ)^{-σ} f[n] ≡ \sum_{p=0}^{\infty} (-1)^p l_p(σ) Δ^p f[n],
```
where ``α_p(σ) = (-1)^p l_p(σ)`` is the ``p^{th}``-order *finite-difference expansion coefficient*
for lagrangian interpolation. Here we define 

```math
l_p(σ) ≡ (σ)_p/p!\,,
``` 
with
```math
(σ)_{p}=\begin{cases}
1 & p=0\\
σ(σ+1)(σ+2)\cdots(σ+p-1) & p>0
\end{cases}
```
being the Pochhammer symbol `CamiMath.pochhammer`. Note that for ``σ = 1`` we find 
``α_p ≡ α_p(1) ≡ (-1)^p``, regaining the expansion coefficients obtained above in the 
example of [`CamiDiff.fdiff_expansion`](@ref). 

For ``-k ≤ σ ≤ 0`` the method can be used for *interpolation* over the grid position interval 
``n ≤ x ≤ n+k``. Outside this interval, in particular for ``σ > 0``, the method amounts to 
*extrapolation*. The method is most accurate for ``-1 ≤ σ ≤ 1`` (corresponding to the grid 
position interval ``n-1 ≤ x ≤ n+1``). Extrapolation to values ``x > n+k`` is not recommended. 

Evaluating the finite-difference expansion up to order ``k`` we obtain  

```math
f[n-σ] =\sum_{p=0}^{k}α_p(σ)Δ^pf[n]
=\sum_{j=0}^{k}F_j^k(σ)f[n+j]
=F^k(σ) \cdot f[n:n+k],
```

where the ``k+1`` *weights*

```math
F_j^k(σ)= \sum_{p=j}^{k} (-1)^k α_p(σ) c_j^p
=\sum_{p=j}^{k} (-1)^j \binom{p}{j}(σ)_p/p!
```
are the *lagrangian interpolation weights* corresponding to the point ``f[n-σ]``.

Symmetry relation:

```math
\bar{F}^k(-k-σ) = F^k(σ)
```

Weight functions:

[`fdiff_expansion_weights(polynom, fwd, reg)`](@ref)
`` → F^k(σ) ≡ [F^k_0(σ),⋯\ F^k_k]``,

where the vector

`polynom = `[`fdiff_interpolation_expansion_polynom(σ, k, fwd)`](@ref)
`` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``  contains the coefficients of the
lagrangian-interpolation expansion.

**Backward difference notation**

Starting from the relation
```math
f[n]=(1-∇)f[n+1].
```
we obtain by formal operator inversion
```math
f[n+1] = (1 - ∇)^{-1} f[n] ≡ \sum_{p=0}^{\infty}∇^p f[n],
```
```math
f[n+2] = (1 - ∇)^{-2} f[n] ≡ \sum_{p=0}^{\infty}p∇^p f[n],
```
```math
\vdots
```

where ``k`` is called the order of the expansion and ``n`` is the reference
index. For interpolation position ``n-σ`` (where σ may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n+σ] = (1 - ∇)^{-σ} f[n] ≡ \sum_{p=0}^{\infty} l_p(σ) ∇^p f[n],
```
where ``β_p(σ) = l_p(σ)`` is the ``p^{th}``-order *finite-difference expansion coefficient* 
for lagrangian interpolation, with ``(σ)_{p}`` being the Pochhammer symbol `CamiMath.pochhammer`.  
Note that for ``σ = 1`` we find ``β_p ≡ β_p(1) ≡ 1``, regaining the expansion coefficients obtained above in the 
example of [`CamiDiff.fdiff_expansion`](@ref). 

For ``-k ≤ σ ≤ 0`` the method can be used for *interpolation* over the grid position interval 
``n-k ≤ x ≤ n``, outside this interval, in particular for ``σ > 0``, the method amounts to 
*extrapolation*. The method is most accurate for ``-1 ≤ σ ≤ 1`` (corresponding to the grid 
position interval ``n-1 ≤ x ≤ n+1``). Extrapolation to values ``x < n-k`` is not recommended. 

Evaluating the finite-difference expansion up to order ``k`` we obtain

```math
f[n+σ] =\sum_{p=0}^{k}β_p(σ)∇^pf[n]
= \sum_{j=0}^{k}B^k_j(σ)f[n-j]
= \bar{B}^k(σ) ⋅ f[n-k:n],
```

where the ``k+1`` *weights*

```math
B^k_j(σ)= \sum_{p=j}^{k} β_p(σ) c_j^p
```

are the corresponding *lagrangian interpolation weights*.  

Symmetry relations:

```math
B^k(σ) = F^k(σ) = \bar{B}^k(-k-σ)
```

```math
\bar{B}^k(σ) = B^k(-k-σ)
```

Weight function:

[`fdiff_expansion_weights(polynom, bwd, rev)`](@ref)
`` → \bar{B}^k(σ) ≡ [B_k^k(σ),⋯\ B_0^k(σ)]``,

where the vector

`polynom = `[`fdiff_interpolation_expansion_polynom(σ, k=3, notation=bwd)`](@ref)
`` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]`` contains the coefficients of the
lagrangian-interpolation expansion.

```@docs
fdiff_interpolation_expansion_polynom(σ::T, notation=bwd; k=3) where T<:Real
fdiff_interpolation_expansion_weights(σ::T, notation=bwd, ordering=rev; k=3) where T<:Real
fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}
```

## Lagrangian differentiation

**Forward difference notation** (`notation = fwd`)

To derive the *lagrangian differentiation* formulas we formally differentiate

```math
f[n-x] = (1+Δ)^{-x} f[n]
```
with respect to ``-x``,

```math
-\frac{df}{dx}[n-x]
=-ln(1+Δ)\ (1+Δ)^{-x}f[n]
=\sum_{q=1}^{k}(-1)^q\tfrac{1}{q}∇^{q}\sum_{p=0}^{k}l_{p}(x)∇^{p}f[n]+⋯.
```

Rewriting the r.h.s. as a single expansion in powers of ``∇``, we obtain

```math
\frac{df}{dx}[n+x]=\sum_{p=1}^{k}β_p(x)∇^{p}f[n]+⋯,
```

where ``β_p(x)`` represents the *finite-difference expansion coefficients*
for *lagrangian differentiation* at position ``n+x``. The coefficients ``β_p(x)`` 
are obtained by polynomial multiplication using the function
[`CamiMath.polynom_product(p1,p2)`](@extref CamiMath.polynom_product), 
where ``p_1`` and ``p_2`` are [`CamiMath.polynom`](@extref CamiMath.polynom) vectors. 
The resulting coefficients are contained in the following [`CamiMath.polynom`](@extref CamiMath.polynom) vector of order ``k``, 

[`fdiff_differentiation_expansion_polynom(k,x)`](@ref) `` → β(x) ≡ [β_0(x),⋯\ β_p(x)]``, with ``β_0(x)≡ 0``.

Substituting the *finite-difference operators*, the *lagrangian derivative* takes the form  

```math
\frac{df}{dx}[n+x]
=\sum_{j=0}^{k}B_j^k(x)f[n-j]
=B^k(x) ⋅ f[n:-1:n-k],
```

where the ``k+1`` *weights*

```math
 B_j^k(x)=\sum_{p=j}^{k}β_p(x)c_{j}^{p}
```

are the ``k^{th}``-order *lagrangian-differentiation weights*

[`fdiff_expansion_weights(β, bwd, reg)`](@ref) `` → B^k(x) ≡ [B^k_0(x),⋯\ B^k_k(x)]``.

After changing dummy index to reverse the summation the expansion becomes

```math
\frac{df}{dx}[n+x]
=\sum_{j=0}^{k}\bar{B}^k_j(x)f[n-k+j]
=\bar{B}^k(x) ⋅ f[n-k:n],
```

where

[`fdiff_expansion_weights(β, bwd, rev)`](@ref) `` → \bar{B}^k(x) ≡ [B^k_k(x),⋯\ B^k_0(x)]``.

**backward difference notation** (`notation = bwd`)

To derive the *lagrangian differentiation* formulas we formally differentiate

```math
f[n+x] = (1 - ∇)^{-x} f[n]
```
with respect to ``x``,

```math
\frac{df}{dx}[n+x]
=-ln(1-∇)\ (1-∇)^{-x}f[n]
=\sum_{q=1}^{k}\tfrac{1}{q}∇^{q}\sum_{p=0}^{k}l_{p}(x)∇^{p}f[n]+⋯.
```

Rewriting the r.h.s. as a single expansion in powers of ``∇``, we obtain

```math
\frac{df}{dx}[n+x]=\sum_{p=1}^{k}β_p(x)∇^{p}f[n]+⋯,
```

where ``β_p(x)`` represents the *finite-difference expansion coefficients*
for *lagrangian differentiation* at position ``n+x``. The coefficients ``β_p(x)`` 
are obtained by polynomial multiplication using the function
[`CamiMath.polynom_product(p1,p2)`](@extref CamiMath.polynom_product), 
where ``p_1`` and ``p_2`` are [`CamiMath.polynom`](@extref CamiMath.polynom) vectors. 
The resulting coefficients are contained in the following [`CamiMath.polynom`](@extref CamiMath.polynom) vector of order ``k``, 

[`fdiff_differentiation_expansion_polynom(k,x)`](@ref) `` → β(x) ≡ [β_0(x),⋯\ β_p(x)]``, with ``β_0(x)≡ 0``.

Substituting the *finite-difference operators*, the *lagrangian derivative* takes the form  

```math
\frac{df}{dx}[n+x]
=\sum_{j=0}^{k}B_j^k(x)f[n-j]
=B^k(x) ⋅ f[n:-1:n-k],
```

where the ``k+1`` *weights*

```math
 B_j^k(x)=\sum_{p=j}^{k}β_p(x)c_{j}^{p}
```

are the ``k^{th}``-order *lagrangian-differentiation weights*

[`fdiff_expansion_weights(β, bwd, reg)`](@ref) `` → B^k(x) ≡ [B^k_0(x),⋯\ B^k_k(x)]``.

After changing dummy index to reverse the summation the expansion becomes

```math
\frac{df}{dx}[n+x]
=\sum_{j=0}^{k}\bar{B}^k_j(x)f[n-k+j]
=\bar{B}^k(x) ⋅ f[n-k:n],
```

where

[`fdiff_expansion_weights(β, bwd, rev)`](@ref) `` → \bar{B}^k(x) ≡ [B^k_k(x),⋯\ B^k_0(x)]``.



```@docs
fdiff_differentiation_expansion_polynom(ξ::T, k=3) where T<:Real
create_lagrange_differentiation_matrix(k::Int)
```

## Integration

```@docs
trapezoidal_epw(k::Int; rationalize=false, devisor=false)
trapezoidal_integration(f, x1, x2, weights)
```

## Adams Method

#### Adams-Bashford expansion

The *Adams-Bashford integration step* is given by the expansion

```math
y[n+1]-y[n] = -\frac{h ∇}{(1-∇)ln(1-∇)}f[n+1]=h (\sum_{p=0}^{\infty}B_p∇^p)f[n+1].
```

A closed expression for the *Adams-Bashford expansion coefficients*, ``B_k``,
is not available. As we already have a finite-difference expansion for the
operator ``(1-∇)^{-1}``,

```math
\frac{1}{1-∇}≡\sum_{p=0}^{\infty}∇^p,
```

we ask for the expansion of

```math
-\frac{∇}{ln(1-∇)}
=(1-\frac{1}{2}∇-\frac{1}{24}∇^2-\frac{1}{12}∇^3+⋯)f[n+1]
= (\sum_{p=0}^{\infty}b_p∇^p)f[n+1].
```

This is known as the *Adams-Moulton expansion*. Its expansion coefficients are
calculated numerically by the function
[`fdiff_adams_moulton_expansion_polynom(k)`](@ref). The *Adams-Bashford expansion* is
obtained as the polynomial product of the two expansions,

```math
(\sum_{p=0}^{\infty}B_p∇^p)f[n+1]
=(\sum_{p=0}^{\infty}∇^p)(\sum_{p=0}^{\infty}b_p∇^p)f[n+1]
=\ ( 1 + \frac{1}{2}∇ + \frac{5}{12}∇^2 + ⋯)f[n+1]
```

where the vector ``β = [B_0,⋯\ B_k]`` contains the *Adams-Bashford expansion coefficients*,
rational numbers generated numerically by the function
[`fdiff_adams_bashford_expansion_polynom(k)`](@ref). Evaluating the finite-difference
expansion up to order ``k`` we obtain (after changing dummy index bring the
summation in forward order)

```math
\sum_{p=0}^{k}B_p∇^pf[n]
=\sum_{p=0}^{k}B_p\sum_{j=0}^{p} c_j^if[n-j]
= \sum_{j=0}^{k}A_j^kf[n-j]
= \sum_{j=0}^{k}A_{k-j}^kf[n-k+j],
```

where the ``A_j^k= \sum_{p=j}^{k} B_pc_j^p`` are the ``(k+1)``-point
*Adams-Bashford integration weights*.

Function:

`β` = [`fdiff_adams_bashford_expansion_polynom(k)`](@ref)
``→ [B_0,⋯\ B_k]``

`adams_bashford_weights`
= [`fdiff_expansion_weights(β, bwd, rev)`](@ref)
 ``→ [A_k^k,⋯\ A_0^k]``

`adams_bashford_weights` = [`create_adams_bashford_weights(k)`](@ref)
``→ [A_k^k,⋯\ A_0^k]``

```@docs
fdiff_adams_bashford_expansion_polynom(k::Int; T=Int, msg=true)
create_adams_bashford_weights(k::Int; rationalize=false, devisor=false, T=Int)
```

### Adams-Moulton expansion

The *Adams-Moulton integration* step is given by the expansion

```math
y[n+1]-y[n]
= -\frac{∇}{ln(1-∇)}f[n+1]
= ( 1 - \frac{1}{2}∇ - \frac{1}{12}∇^2 - \frac{1}{24}∇^3 +⋯)f[n+1].
```

For the evaluation of the integration step we limit the summation to ``k+1``
terms (order ``k``),

```math
y[n+1]-y[n]= (\sum_{p=0}^{k}b_p∇^p)f[n+1]+⋯.
```

where the vector ``β = [b_0,⋯\ b_k]`` contains the *Adams-Moulton expansion coefficients*,
rational numbers generated numerically by the function
[`fdiff_adams_moulton_expansion_polynom(k)`](@ref). Extracting the greatest
common denominator, ``1/D``, the step becomes

```math
y[n+1]-y[n]= \frac{1}{D}(\sum_{p=0}^{k}b_p′∇^p)f[n+1]+⋯,
```

where ``b_0′,⋯\ b_k′`` are integers and
``b_p=b_p′/D``. In practice the expansion is restricted to ``k<18``
(as limited by integer overflow). Note that this limit is much higher than
values used in calculations (typically up to ``k = 10``). Evaluating the
finite-difference expansion up to order ``k`` we obtain (after changing
dummy index bring the summation in forward order)

```math
\sum_{p=0}^{k}b_p∇^pf[n]
=\sum_{p=0}^{k}b_p\sum_{j=0}^{p} c_j^if[n-j]
= \sum_{j=0}^{k}a_j^kf[n-j]
= \sum_{j=0}^{k}a_{k-j}^kf[n-k+j],
```

where the ``a_j^k= \sum_{p=j}^{k} b_pc_j^p`` are the ``(k+1)``-point
*Adams-Moulton integration weights*.

Functions:

`β` = [`fdiff_adams_moulton_expansion_polynom(k)`](@ref) ``→ [b_0,⋯\ b_k]``

`adams_moulton_weights`
= [`fdiff_expansion_weights(β, bwd, rev)`](@ref)
``→ [a_k^k,⋯\ a_0^k]``.

`adams_moulton_weights` = [`create_adams_moulton_weights(k)`](@ref)
``→ [a_k^k,⋯\ a_0^k]``

```@docs
fdiff_adams_moulton_expansion_polynom(k::Int; T=Int, msg=true)
create_adams_moulton_weights(k::Int; rationalize=false, devisor=false, T=Int)
```

# Application

The three elementary operations are:

[`grid_interpolation(f, rval, grid)`](@ref)

[`grid_differentiation(f, grid)`](@ref)

[`grid_integration(f, grid)`](@ref)

```@docs
grid_interpolation(f::Vector{T}, rval::T, grid::Grid{T}; k=5) where T<:Real
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real
grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
```