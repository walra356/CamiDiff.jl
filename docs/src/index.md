```@meta
CurrentModule = CamiDiff
```

# CamiDiff.jl

A [Julia](http://julialang.org) package for finite-difference analysis

---
## Table of contents

```@contents
```

## Introduction

In `CamiDiff` we present general purpose tools for the finite difference analysis 
of *Real* analytic functions of a single variable, which we denote by ``f(x)``. 

The finite-difference analysis starts by discretization of ``f(x)`` onto a grid of ``N`` points,
```math
f(x) ↦ f[n],
```
where ``f[n]`` is a discrete function representing the function ``f`` at position 
```math
x = t(n) ≡ (n−u) * h + x0.
```
Here the *linear function* ``t(n)`` is called the *step function*, with ``u`` the *index base*, 
``h`` the *step size* and ``x0`` the *offset* of the grid. As [Julia](http://julialang.org) is a unit-based 
language (``u == 1``), we have ``f[1] = f(x0)``. 

The discretization map is defined by 
```math
f[n] = r0 g(t[n]),
```
where ``g(t)`` is a (generally nonlinear) function, the *grid function*, and ``r0`` the 
*scaling factor*.

`CamiDiff` was developped for the use of spherical coordinates, i.e., for case of 
zero offset (``x0 = 0``). For this case we write ``f(r)`` rather than ``f(x)``, with
```math
r = t(n) ≡ (n−1) * h.
```

## Grid

The `Grid` object is the backbone for the numerical procedure on a non-uniform
grid. Its principal fields are `grid.r` and `grid.r′`, which are discrete
functions of `N` elements representing the grid function and its derivative.

```@docs
Grid{T}
gridname(ID::Int)
gridfunction(ID::Int, n::Int, h::T; p=5, polynom=[0,1], deriv=0) where T <: Real
castGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, polynom=[0,1], epn=5, k=7, msg=true)
findIndex(rval::T, grid::Grid{T}) where T<:Number
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=3) where T<:Real
grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
```



## Finite-difference methods

### Finite differences

Consider the analytical function ``f`` tabulated in *forward order*  
(growing index) at ``n`` positions on a *uniform grid*.

**Forward difference notation**

In *forward difference* notation, the *finite difference* of two adjacent
values on the grid is defined as

```math
Δ f[n] = f[n+1]-f[n],
```

where ``Δ`` is the forward difference operator. By a formal inversion  
procedure we find

```math
f[n-1]=(1+Δ)^{-1}f[n]=(1-Δ+Δ^2-Δ^3+⋯)f[n],
```

where ``Δ^k`` is the  ``k^{th}``-*order forward difference* defined as
a *weighted sum* over the function values ``f[n:n+k]`` (involving
``k+1`` points),

```math
Δ^k f[n] = c_{k}^kf[n] + c_{k-1}^kf[n+1] + ⋯  + f[n+k]
= \sum_{j=0}^{k} c_{k-j}^k f[n+j].
```
The ``k+1`` coefficients

```math
c_{k-j}^{k}=(-1)^{k-j}\binom{k}{j}
```

are the *summation weights* (short: *weights*) which define the summation.

**Backward difference notation**

In *backward difference* notation, the *finite difference* of two adjacent
values on the grid is defined as

```math
∇ f[n] = f[n]-f[n-1],
```

where ``∇`` is the backward difference operator.  By a formal inversion  
procedure we find

```math
f[n+1]=(1-∇)^{-1}f[n]=(1+∇+∇^2+∇^3+⋯)f[n],
```

where ``∇^k`` is the  ``k^{th}``-*order backward difference* defined as
a *weighted sum* over the function values tabulated in backward order,
``f[n:-1:n-k]`` (involving ``k+1`` points),

```math
∇^k f[n] = f[n] + c_1^kf[n-1] + ⋯ + c_k^kf[n-k]
= \sum_{j=0}^{k} c_j^kf[n-j],
```

where the ``k+1`` coefficients

```math
c_{j}^{k}=(-1)^{j}\binom{k}{j}
```

are the *summation weights* (short: *weights*) which define the summation.
Note the special cases ``c_{0}^{k}≡1``, ``c_{k}^{k}≡(-1)^{k}`` and the symmetry
relation

```math
c_{k-j}^k=(-1)^k c_j^k.
```

Coefficients:  

[`fdiff_weight(k,j)`](@ref) `` → c_j^k=(-1)^j\binom{k}{j}``

```@docs
# isforward(notation)
# isregular(ordering)
fdiff_weight(k::Int, j::Int)
```

### Finite difference expansions

Finite-difference calculus builds on the *finite-difference expansion*.

**Forward difference notation**

In terms of forward differences the expansion takes the form

```math
\sum_{p=0}^{\infty}α_{p}Δ^{p}f[n]
=\sum_{p=0}^{k}α_{p}Δ^{p}f[n]+⋯.
```

A finite-difference expansion truncated at order ``k`` is defined
by ``k+1`` *finite-difference expansion coefficients*, represented by the
vector ``α = [α_{0},⋯\ α_{k}]``. It takes some bookkeeping to rewrite the
expansion as a *weighted sum* over the ``k+1``
*function values in forward tabulated form* ``f[n:n+k]``.
Substituting the finite difference expression for ``Δ^k``, we obtain

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

with ``j=0,⋯\ k``. In inner product form the expansion becomes
```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n]
=\sum_{j=0}^{k}F_{j}^{k}f[n+j]
=F^{k} \cdot f[n:n+k],
```

where ``F^k  ≡ [F_0^k,⋯\ F_k^k]``.

```math
f[n:n+k] = \left[\begin{array}{c}
f[n]\\
\vdots\\
f[n+k]
\end{array}\right].
```

Coefficients:

[`fdiff_expansion_weights(coeffs, fwd, reg)`](@ref)
``→ F^k ≡ [F_0^k,⋯\ F_k^k]``,

where the `coeffs` ``  α ≡ [α_0,⋯\ α_k]`` are user supplied to define the
expansion.

**Backward difference notation**

In terms of backward differences the expansion takes the form

```math
\sum_{p=0}^{\infty}β_{p}∇^{p}f[n]=\sum_{p=0}^{k}β_{p}∇^{p}f[n]+⋯.
```

In this case the ``k^{th}``- order *finite-difference expansion* is defined
by the vector ``β = [β_{0},⋯\ β_{k}]``. The expansion can written as
*weighted sum* over the ``k+1`` *function values in backward tabulated form*
``f[n:-1:n-k]``. Substituting the finite
difference expression for ``∇^k``, we obtain

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
with ``j=0,⋯\ k``. In inner product form the expansion becomes

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n]
=\sum_{j=0}^k B_j^k f[n-j]
=\bar{B}^k \cdot f[n-k:n],
```

where the *weights vector* ``\bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]`` contains
the weights in backward order.

In general there is *no simple symmetry relation* between
``B^k`` and ``F^k``.

Coefficients:

[`fdiff_expansion_weights(coeffs, bwd, rev)`](@ref)
`` → \bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]``,

where the `coeffs`  ``  β ≡ [β_0,⋯\ β_k]`` are user supplied to
define the expansion.

```@docs
fdiff_expansion_weights(coeffs, notation=CamiMath.bwd, ordering=CamiMath.rev)
fdiff_expansion(coeffs, f, notation=CamiMath.bwd)
```

### Lagrange-polynomial interpolation/extrapolation

The Lagrange polynomial of degree k on a uniform grid is the polynomial running through k+1 subsequent points on the grid. We derive expressions for interpolation/extrapolation in both forward- and backward-difference notation. Beware that Lagrange interpolation becomes inaccurate if the tabulated function cannot be approximated by a polynomial of degree k.

**Forward difference notation**

Starting from the relation
```math
f[n]=(1+Δ)f[n+1],
```
we obtain by formal operator inversion
```math
f[n+1] = (1 + Δ)^{-1} f[n] \equiv \sum_{p=0}^{\infty}(-1)^p Δ^p f[n],
```
```math
f[n+2] = (1 + Δ)^{-2} f[n] \equiv \sum_{p=0}^{\infty}(-1)^p pΔ^p f[n],
```
```math
\vdots
```
where ``k`` is called the order of the expansion and ``n`` is the reference
index. For interpolation position ``n-σ`` (where σ may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n-σ] = (1 + Δ)^{-σ} f[n] \equiv \sum_{p=0}^{\infty} α_p(σ) Δ^p f[n],
```
where

```math
α_p(σ) ≡ (-1)^p(σ)_p/p!
```
is the ``p^{th}``-order *finite-difference expansion coefficient*
for lagrangian lagrangian interpolation over the
interval ``-k ≤σ ≤0\ \ (n \le n-σ \le n+k)``,

```math
(σ)_{p}=\begin{cases}
1 & p=0\\
σ(σ+1)(σ+2)\cdots(σ+p-1) & p>0
\end{cases}
```
being the Pochhammer symbol `CamiMath.pochhammer`. For ``σ`` outside 
the interpolation interval the method corresponds to *extrapolation* along the
Lagrange polynomial. Evaluating the finite-difference expansion up to
order ``k`` we obtain  

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
are the *lagrangian interpolation weights* corresponding to the point
``f[n-σ]``.

Symmetry relation:

```math
\bar{F}^k(-k-σ) = F^k(σ)
```

Weight functions:

[`fdiff_expansion_weights(coeffs, fwd, reg)`](@ref)
`` → F^k(σ) ≡ [F^k_0(σ),⋯\ F^k_k]``,

where the vector

`coeffs = `[`fdiff_interpolation_expansion_coeffs(σ, k, fwd)`](@ref)
`` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``  contains the coefficients of the
lagrangian-interpolation expansion.

**Backward difference notation**

Starting from the relation
```math
f[n]=(1-∇)f[n+1].
```
we obtain by formal operator inversion
```math
f[n+1] = (1 - ∇)^{-1} f[n] \equiv \sum_{p=0}^{\infty}∇^p f[n],
```
```math
f[n+2] = (1 - ∇)^{-2} f[n] \equiv \sum_{p=0}^{\infty}p∇^p f[n],
```
```math
\vdots
```

where ``k`` is called the order of the expansion and ``n`` is the reference
index. For interpolation position ``n-σ`` (where σ may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n+σ] = (1 - ∇)^{-σ} f[n] \equiv \sum_{p=0}^{\infty} β_p(σ) ∇^p f[n],
```
where

```math
β_p(σ) ≡ (σ)_p/p! = (-1)^p α_p(σ)
```

is the ``p^{th}``-order *finite-difference expansion coefficient* for
lagrangian interpolation over the interval ``-k ≤σ ≤0\ \ (n-k \le n+σ \le n)``,
with

```math
(σ)_{p}=\begin{cases}
1 & p=0\\
σ(σ+1)(σ+2)\cdots(σ+p-1) & p>0
\end{cases}
```
being the Pochhammer symbol `CamiMath.pochhammer`. For ``σ`` outside 
the interpolation interval the method corresponds to *extrapolation* along the
Lagrange polynomial. Evaluating the finite-difference expansion up to
order ``k`` we obtain

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

[`fdiff_expansion_weights(coeffs, bwd, rev)`](@ref)
`` → \bar{B}^k(σ) ≡ [B_k^k(σ),⋯\ B_0^k(σ)]``,

where the vector

`coeffs = `[`fdiff_interpolation_expansion_coeffs(σ, k=3, notation=bwd)`](@ref)
`` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]`` contains the coefficients of the
lagrangian-interpolation expansion.

```@docs
fdiff_interpolation_expansion_coeffs(ξ::T, k=3, notation=CamiMath.bwd) where T<:Real
fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}
```

### Lagrangian differentiation

To derive the *lagrangian differentiation* formulas we formally differentiate

```math
f[n+x] = (1 - ∇)^{-x} f[n]
```

with respect to ``x``.

```math
\frac{df}{dx}[n+x]
=-ln(1-∇)\ (1-∇)^{-x}f[n]
=\sum_{q=1}^{k}\tfrac{1}{q}∇^{q}\sum_{p=0}^{k}l_{p}(x)∇^{p}f[n]+⋯.
```

Rewriting the r.h.s. as a single summation in powers of ``∇`` for given values
of ``n`` and ``x`` we obtain an expression of the form

```math
\frac{df}{dx}[n+x]=\sum_{p=1}^{k}β_p(x)∇^{p}f[n]+⋯,
```

where ``β_p(x)`` represents the *finite-difference expansion coefficients*
for *lagrangian differentiation* at position ``n+x``. These
coefficients are determined numerically by polynomial multiplication. As the
expansion algorith requires the presentce of a ``β_0(x)`` coefficient we add
a (vanishing) ``p=0`` term, ``β_0(x)\equiv 0``. The corresponding coefficient
vector is given by [`fdiff_differentiation_expansion_coeffs(k,x)`](@ref).
Evaluating the finite-difference expansion up to order ``k`` we obtain

```math
\frac{df}{dx}[n+x]
=\sum_{p=0}^{k}β_p(x)∇^pf[n]
=\sum_{j=0}^{k}B_j^k(x)f[n-j]
=B^k(x) ⋅ f[n:-1:n-k],
```

where the ``k+1`` *weights*

```math
 B_j^k(x)=\sum_{p=j}^{k}β_p(x)c_{j}^{p}
```

are the ``k^{th}``-order lagrangian differentiation weights. After changing
dummy index to reverse the summation the expansion becomes

```math
\frac{df}{dx}[n+x]
=\sum_{j=0}^{k}\bar{B}^k_j(x)f[n-k+j]
=\bar{B}^k(x) ⋅ f[n-k:n].
```

Functions:

[`fdiff_expansion_weights(β, bwd, reg)`](@ref)
`` → B^k(x) ≡ [B^k_0(x),⋯\ B^k_k(x)]``

[`fdiff_expansion_weights(β, bwd, rev)`](@ref)
`` → \bar{B}^k(x) ≡ [B^k_k(x),⋯\ B^k_0(x)]``

where

[`fdiff_differentiation_expansion_coeffs(o, k)`](@ref)
``→ β ≡ [β_0(x),⋯\ β_k(x)]``.

```@docs
fdiff_differentiation_expansion_coeffs(ξ::T, k=3) where T<:Real
fdiff_differentiation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}
create_lagrange_differentiation_matrix(k::Int)
```

### Integration

```@docs
trapezoidal_epw(k::Int; rationalize=false, devisor=false)
trapezoidal_integration(f, x1, x2, weights)
```
