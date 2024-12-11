```@meta
CurrentModule = CamiDiff
```

# CamiDiff.jl

A [Julia](http://julialang.org) package for one-dimensional finite-difference analysis

---

## Install

The package is installed using the Julia package manager

```
julia> using Pkg; Pkg.add("CamiMath")

julia> using CamiMath
```
## Table of contents

```@contents
```

## Introduction

`CamiDiff` is a package for finite-difference analysis of *real analytic functions of a single variable*. 

### Finite differences

Finite-difference calculus is based on the manipulation of finite differences of analytic functions. To introduce the subject we consider the analytic function ``f(x)``. The finite difference of ``f(x+h)`` and ``f(x)`` is called the *forward difference* and is defined as 
```math
Δ f(x)=f(x+h)-f(x).
```
Here, ``Δ`` is called the *forward-difference operator*. Likewise one defines the *backward-difference operator* ``∇``, 
```math
∇ f(x)=f(x)-f(x-h).
```
We first focus on *forward differences*. The derivative of ``f(x)`` is given by 
```math
f^′(x)=\underset{h→0}{\mathrm{lim}}\,\frac{f(x+h)-f(x)}{h}=\underset{Δ x→0}{\mathrm{lim}}\,\frac{Δ f(x)}{Δ x},
```
where ``h ≡ Δx ≥ 0`` is the *difference interval*. Introducing the differential operator, ``f^′(x) ≡ Df(x)``, we have 
```math
D≡\frac{d}{dx}=\underset{Δ x→0}{\mathrm{lim}}\,\frac{Δ}{Δ x}=\underset{h→0}{\mathrm{lim}}\,\frac{Δ}{h}.
```
### Translation operators 

**Forward difference notation**

With regard to *forward differences* we rewrite the forward difference definition in the form of a *forward translation*,
```math
f(x+h)=(1+Δ)f(x),
```
where ``T≡(1+Δ)`` is the *forward translation operator*. This operator shifts the function over the infinitesimal interval ``h`` to larger values of ``x``. The translation operator can be expressed in terms of the differential operator as follows by Taylor expansion of ``f(x)`` about the point ``x``, 
```math
f(x± h)=(1± hD+\tfrac{1}{2}h^2D^2±\tfrac{1}{3!}h^3D^3+⋯)f(x)=e^{± hD}f(x).
```
Comparing the operator expression for the *forward* translation with the Taylor expansion we obtain, by formal inversion of the operator ``T``, an operator identity for the inverted translation operator ``T^{-1}``,    
```math
T≡(1+Δ)=e^{hD}\,\,\,⇒\,\,\,T^{-1}=e^{-hD}=(1+Δ)^{-1}.
```
With this procedure, the explicit dependence on ``h`` can be replaced by an implicit dependence on ``h`` through an expansion in powers of ``Δ`` ,
```math
f(x-h)=(1+Δ)^{-1}f(x)=(1-Δ+Δ^{2}-Δ^3+⋯)f(x).
```
By choosing the proper expansion order, ``f(x-h)`` can be approximated to any desired level of accuracy.

**Backward difference notation**

Likewise, for *backward differences*, we rewrite the backward-difference definition in the form 
```math
f(x-h)=(1-∇)f(x),
```
where ``B≡(1-∇)`` is the *backward translation operator*. Comparing this *backward* translation with the Taylor expansion we obtain, by formal inversion of the operator ``B``, an operator identity for the *forward* translation operator ``T``, 
```math
B≡(1-∇)=e^{-hD}=T^{-1}\,\,\,⇒\,\,\,T=e^{hD}=(1+∇)^{-1}.
```
Note how the *backward* translation operator was identified with the inverse *forward* translation operator, ``B=T^{-1}``. When using backward differences, the explicit dependence on ``h`` can be replaced by an implicit dependence on ``h`` through an expansion in powers of ``∇``, 
```math
f(x+h)=(1-∇)^{-1}f(x)=(1+∇+∇^{2}+∇^3+⋯)f(x).
```
By choosing the proper expansion order, ``f(x+h)`` can be approximated to any desired level of accuracy. 

### Discretization

Finite-difference analysis starts by discretization of the function ``f(x)`` onto a [`Grid`](@ref) of ``N`` points, 
which is based on the map ``n ↦ x`` and defined by the discrete function
```math
x[n] = s_0 * g(t[n]).
```
Here ``g(t)`` is called the [`gridfunction`](@ref) and ``s_0`` the *scaling factor*. The [`gridfunction`](@ref) 
is a (generally nonlinear) analytic function *running through the origin*, ``g(0) = 0``. Its argument is 
the *ticks function*
```math
t[n] ≡ (n−u) * h,
```
which is a *linear* function, with ``u`` the *index base* and ``h`` the *step size*. Writing
```math
f[n] = f(x[n]),
```
we recognize in ``f[n]`` a discrete function representing the function ``f(x)`` at position ``x[n]``, with ``n = 1, ⋯ N``. 

Note that ``h`` determines the *coarseness* of the [`Grid`](@ref). The results of a finite-difference calculation 
on a coarse grid will be less accurate than those on a fine grid, but the algorithm remains the same because 
finite-difference expansions only depend on ``h`` *implicitely*. Since [Julia](http://julialang.org) uses 
unit-based indexing (``u = 1``), the index convention implies ``f[1] = f(0)``.  

NB. The current implementation of `CamiDiff` was developped for grid functions defined on the domain ``[0, ∞)``. 
For this case we use the variable ``r`` rather than ``x``, writing ``f(r)`` rather than ``f(x)``, with ``r ≥ 0``.

## Grid

The [`Grid`](@ref) object is the backbone for the numerical procedure on a (generally) non-uniform
grid. Its principal fields are `grid.r`, `grid.r′` and `grid.r′′` which are discrete
functions of `N` elements representing the grid function and its first two derivatives.

```@docs
Grid{T}
castGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, polynom=[0,1], epn=5, k=7, msg=true)
gridfunction(ID::Int, n::Int, h::T; p=5, polynom=[0,1], deriv=0) where T <: Real
gridname(ID::Int)
findIndex(rval::T, grid::Grid{T}) where T<:Number
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real
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

[`fdiff_expansion_weights(polynom, fwd, reg)`](@ref)
``→ F^k ≡ [F_0^k,⋯\ F_k^k]``,

where `polynom` is the [`CamiMath.polynom`](https://walra356.github.io/CamiMath.jl/stable/#CamiMath.polynom) vector ``  α ≡ [α_0,⋯\ α_k]``. This `polynom` has to be supplied by the user to define the expansion under consideration. Some common cases are:

  interpolation expansion: [`fdiff_interpolation_expansion_polynom(ξ, k, fwd)`](@ref)

differentiation expansion: [`fdiff_differentiation_expansion_polynom(ξ, k, fwd)`](@ref)

 Adams-Bashford expansion: [`fdiff_adams_bashford_expansion_polynom(k, fwd)`](@ref)

  Adams-Moulton expansion: [`fdiff_adams_moulton_expansion_polynom(k, fwd)`](@ref)

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

[`fdiff_expansion_weights(polynom, bwd, rev)`](@ref)
`` → \bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]``,

where `polynom` is the [`CamiMath.polynom`](https://walra356.github.io/CamiMath.jl/stable/#CamiMath.polynom) vector  ``  β ≡ [β_0,⋯\ β_k]``. This `polynom` has to be supplied by the user to define the expansion under consideration. Some common cases are:

  interpolation expansion: [`fdiff_interpolation_expansion_polynom(ξ, k, bwd)`](@ref)

differentiation expansion: [`fdiff_differentiation_expansion_polynom(ξ, k, bwd)`](@ref)

 Adams-Bashford expansion: [`fdiff_adams_bashford_expansion_polynom(k, bwd)`](@ref)

  Adams-Moulton expansion: [`fdiff_adams_moulton_expansion_polynom(k, bwd)`](@ref)

```@docs
fdiff_expansion_weights(polynom, notation=CamiMath.bwd, ordering=CamiMath.rev)
fdiff_expansion(polynom, f, notation=CamiMath.bwd)
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
f[n+1] = (1 + Δ)^{-1} f[n] ≡ \sum_{p=0}^{\infty}(-1)^p Δ^p f[n],
```
```math
f[n+2] = (1 + Δ)^{-2} f[n] ≡ \sum_{p=0}^{\infty}(-1)^p pΔ^p f[n],
```
```math
\vdots
```
where ``k`` is called the order of the expansion and ``n`` is the reference
index. For interpolation position ``n-σ`` (where σ may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n-σ] = (1 + Δ)^{-σ} f[n] ≡ \sum_{p=0}^{\infty} α_p(σ) Δ^p f[n],
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
f[n+σ] = (1 - ∇)^{-σ} f[n] ≡ \sum_{p=0}^{\infty} β_p(σ) ∇^p f[n],
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

[`fdiff_expansion_weights(polynom, bwd, rev)`](@ref)
`` → \bar{B}^k(σ) ≡ [B_k^k(σ),⋯\ B_0^k(σ)]``,

where the vector

`polynom = `[`fdiff_interpolation_expansion_polynom(σ, k=3, notation=bwd)`](@ref)
`` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]`` contains the coefficients of the
lagrangian-interpolation expansion.

```@docs
fdiff_interpolation_expansion_polynom(ξ::T, k=3, notation=CamiMath.bwd) where T<:Real
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
a (vanishing) ``p=0`` term, ``β_0(x)≡ 0``. The corresponding coefficient
vector is given by [`fdiff_differentiation_expansion_polynom(k,x)`](@ref).
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

[`fdiff_differentiation_expansion_polynom(o, k)`](@ref)
``→ β ≡ [β_0(x),⋯\ β_k(x)]``.

```@docs
fdiff_differentiation_expansion_polynom(ξ::T, k=3) where T<:Real
create_lagrange_differentiation_matrix(k::Int)
```

### Integration

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

## Index

```@index
```