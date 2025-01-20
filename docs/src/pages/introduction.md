# Introduction

## Finite differences

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
## Translation operators 

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

## Discretization

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

## Finite difference expansions

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

