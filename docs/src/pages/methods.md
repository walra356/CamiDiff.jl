## Finite-difference methods

### Finite differences

Consider the analytical function ``f`` tabulated in *forward order*  
(growing index) at ``n`` positions on a *uniform grid*.

**Forward difference notation**

The *forward translation* from position ``n-1`` to position ``n`` on 
the grid is expressed by the relation

```math
f[n] = (1 + Δ) f[n-1] ,
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
Δ^k f[n] = f[n] + c_{k-1}^kf[n+1] + ⋯  + f[n+k]
= \sum_{j=0}^{k} c_{k-j}^k f[n+j],
```
where the ``k+1`` coefficients

```math
c_{k-j}^{k}=(-1)^{k-j}\binom{k}{j}
```

are the *summation weights* (short: *weights*) which define the summation, 
with special values ``c_{0}^{k}≡1``, ``c_{k}^{k}≡1`` 

**Backward difference notation**

The *backward translation* from position ``n-1`` to position ``n`` on 
the grid is expressed by the relation

```math
f[n-1] = (1 - ∇) f[n],
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
∇^k f[n] = f[n] + c_1^kf[n-1] + ⋯ + f[n-k]
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

