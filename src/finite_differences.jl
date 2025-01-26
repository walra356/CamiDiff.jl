# SPDX-License-Identifier: MIT

# Copyright (c) 2024 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ==============================================================================
#                               finite_differences.jl
# ==============================================================================

# ------------------------------------------------------------------------------
#                               fdiff_weight(k, j)
# ------------------------------------------------------------------------------

@doc raw"""
    fdiff_weight(k::Int, j::Int)

Finite difference weight coefficient,

```math
c_{j}^{k}=(-1)^{j}\binom{k}{j}.
```
#### Example:
```
julia> c(k,j) = fdiff_weight(k,j);

julia> a = [[c(k,j) for j=0:k] for k=0:3]
4-element Vector{Vector{Int64}}:
 [1]
 [1, -1]
 [1, -2, 1]
 [1, -3, 3, -1]

julia> b = [[c(k,k-j) for j=0:k] for k=0:3]
4-element Vector{Vector{Int64}}:
 [1]
 [-1, 1]
 [1, -2, 1]
 [-1, 3, -3, 1]

julia> b == reverse.(a)
true
```
"""
function fdiff_weight(k::Int, j::Int)

    return Base.iseven(j) ? Base.binomial(k,j) : -Base.binomial(k,j)

end

# ------------------------------------------------------------------------------
#                  fdiff_expansion_weights(polynom, bwd, rev)
# ------------------------------------------------------------------------------

function reg_fwd_expansion_weights(polynom)

    k = Base.length(polynom)-1

    o = [Base.sum([polynom[p+1] * fdiff_weight(p, p-j) for p=j:k]) for j=0:k]

    return o

end
# ..............................................................................
function rev_fwd_expansion_weights(polynom)

    k = Base.length(polynom)-1

    o = [Base.sum([polynom[p+1] * fdiff_weight(p, p-j) for p=j:k]) for j=k:-1:0]

    return o

end
# ..............................................................................
function reg_bwd_expansion_weights(polynom)

    k = Base.length(polynom)-1

    o = [sum([polynom[p+1] * fdiff_weight(p, j) for p=j:k]) for j=0:k]

    return o

end
# ..............................................................................
function rev_bwd_expansion_weights(polynom)

    k = Base.length(polynom)-1

    o = [sum([polynom[p+1] * fdiff_weight(p, j) for p=j:k]) for j=k:-1:0]

    return o

end
# ..............................................................................
@doc raw"""
    fdiff_expansion_weights(polynom [, notation=bwd [, ordering=rev]])

Weights vector corresponding to the expansion coefficient vector `polynom`
of a (user-defined) finite-difference expansion.

**Forward-difference notation** (`notation = fwd`)

The weights vector ``F^k ≡ [F_k^k,⋯\ F_0^k]`` corresponds in this case to the expansion coefficient 
vector `` α ≡ [α_0,⋯\ α_k]`` of the ``k^{th}``-order *forward-difference* expansion (`polynom = α`).

```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n]
=\sum_{j=0}^{k}F_{j}^{k}f[n+j]
=F^{k} \cdot f[n:n+k],
```

where ``f[n:n+k]`` are elements of the
analytic function ``f`` tabulated in regular=*forward* order.

[`fdiff_expansion_weights(α, fwd, reg)`](@ref)
``→ F^k ≡ [F_0^k,⋯\ F_k^k]``,

where `` α ≡ [α_0,⋯\ α_k]`` has to be supplied by the user in combination with `fwd` to
indicate that the weights must be evaluated in forward-difference notation.

**Backward difference notation** (`notation = bwd`)

The weights vector ``\bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]`` corresponds in this case to the expansion c
oefficient vector ``β ≡ [β_0,⋯\ β_k]`` of the ``k^{th}``-order *backward-difference* expansion (`polynom = β`).

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n]
=\sum_{j=0}^{k}B_{j}^kf[n-j]
=\bar{B}^k \cdot f[n-k:n].
```

where ``f[n-k:n]`` are elements of the
analytic function ``f`` tabulated in *forward* order.

[`fdiff_expansion_weights(β, bwd, rev)`](@ref)
`` → \bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]``,

where `` β ≡ [β_0,⋯\ β_k]`` has to be supplied by the user in combination with `bwd` to
indicate that the weights must be evaluated in backward-difference notation.
#### Example:

Consider the expansions,

```math
\begin{aligned}
f[n-1]=(1+Δ)^{-1}f[n]=(1-Δ+Δ^2-Δ^3+⋯)f[n]&=F^{k} \cdot f[n:n+k],\\
f[n+1]=(1-∇)^{-1}f[n]=(1+∇+∇^2+∇^3+⋯)f[n]&=\bar{B}^k \cdot f[n-k:n].
\end{aligned}
```
```
julia> f = [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100];

julia> α = [1,-1,1,-1,1];

julia> Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
Fk = [5, -10, 10, -5, 1]

julia> β = [1,1,1,1,1];

julia> revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
revBk = [1, -5, 10, -10, 5]

julia> sum(Fk .* f[7:11]) # inner product
25

julia> sum(revBk .* f[1:5]) # inner product
25

julia> f[6]
25
```
"""
function fdiff_expansion_weights(polynom, notation=bwd, ordering=rev)

    if CamiMath.isforward(notation)

        o = CamiMath.isregular(ordering) ? reg_fwd_expansion_weights(polynom) :
                                           rev_fwd_expansion_weights(polynom)
    else

        o = CamiMath.isregular(ordering) ? reg_bwd_expansion_weights(polynom) :
                                           rev_bwd_expansion_weights(polynom)
    end

    return o

end

# ------------------------------------------------------------------------------
#                  fdiff_expansion(polynom, f, notation=bwd)
# ------------------------------------------------------------------------------

@doc raw"""
    fdiff_expansion(polynom, f [, notation=bwd])

Finite difference expansion evaluated to ``k^{th}`` order for the analytical 
function ``f``, tabulated in *regular order* (growing index) at ``k+1`` positions 
on a [`Grid`](@ref). The expansion coefficients are specified by the vector 
`polynom`. By default the expansion is calculated in backward-difference notation (`bwd`). 

**Forward difference notation** (`notation = fwd`)
```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n] = F^{k} \cdot f[n:n+k],
```
where ``f[n:n+k]`` are elements of the analytical function ``f`` (tabulated in 
*forward* order) and `polynom` is the (user-supplied) expansion coefficient vector 
``α ≡ [α_0,⋯\ α_k]``.

**Backward difference notation** (`notation = bwd`)
```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n] = \bar{B}^k \cdot f[n-k:n].
```
where ``f[n-k:n]`` are elements of the
analytical function ``f`` (tabulated in *forward* order) and `polynom` is the (user-supplied)
expansion coefficient vector ``β ≡ [β_0,⋯\ β_k]``.

NB. The vector `polynom` determines the order of the expansion, 
``k+1 = \rm{length}(α) = \rm{length}(β)``. The weights vectors ``F^k`` and 
``\bar{B}^k`` are *internally generated* by the function
[`fdiff_expansion_weights(polynom, notation, ordering)`](@ref).

#### Examples:
Consider the function ``f(x)=x^2`` and the expansions,
```math
\begin{aligned}
f[n-1]=(1+Δ)^{-1}&=(1-Δ+Δ^2-Δ^3+⋯)f[n]=F^{k} \cdot f[n:n+k],\\
f[n]=(1+Δ)^{-1}&=(1-Δ+Δ^2-Δ^3+⋯)f[n+1]=F^{k} \cdot f[n+1:n+k+1],
\end{aligned}
```
```math
\begin{aligned}
f[n+1]=(1-∇)^{-1}&=(1+∇+∇^2+∇^3+⋯)f[n]=\bar{B}^k \cdot f[n-k:n],\\
f[n]=(1-∇)^{-1}&=(1+∇+∇^2+∇^3+⋯)f[n-1]=\bar{B}^k \cdot f[n-k-1:n-1]
\end{aligned}
```
Note that, in these examples, to fourth order in the expansion ``(k=4)``, 
the forward- and backward-difference coefficient vectors `polynom` are 
given by ``α=[1,-1,1,-1,1]`` and ``β=[1,1,1,1,1]``, respectively. 
```
julia> f = [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100];

julia> α = [1,-1,1,-1,1];

julia> β = [1,1,1,1,1];

julia> fdiff_expansion(α,f[7:11],fwd)
25

julia> fdiff_expansion(β,f[1:5],fwd)
25

julia> f[6]
25
```
In these cases the results are exact because the function is quadratic and
the expansion is third order (based on the polynomial of ``k^{th}`` degree
running through the ``k+1`` points of the tabulated function). Compare with 
the example of [`fdiff_interpolation(f, v, k=3)`](@ref).
"""
function fdiff_expansion(polynom, f, notation=bwd)

    ordering = CamiMath.isforward(notation) ? reg : rev
    w = fdiff_expansion_weights(polynom, notation, ordering)

    #return LinearAlgebra.dot(w, f)
    return sum(w .* f)


end

# ------------------------------------------------------------------------------
#           fdiff_interpolation_expansion_polynom(k, x, notation=bwd)
# ------------------------------------------------------------------------------

function fwd_interpolation_expansion_polynom(σ::T; k=3) where T<:Real

    o = Base.ones(T,k+1)
    σ == 0 ? (for p=2:k+1; o[p] = 0 end) :
             (for p=1:k; o[p+1] = o[p]*(-σ-p+1)/p end)

    return o

end
#...............................................................................
function bwd_interpolation_expansion_polynom(σ::T; k=3) where T<:Real

    o = Base.ones(T,k+1)
    σ == 0 ? (for p=2:k+1; o[p] = 0 end) :
             (for p=1:k; o[p+1] = o[p]*(σ+p-1)/p end)

    return o

end
function fwd_interpolation_expansion_polynom1(ξ::T, k=3) where T<:Real

    o = Base.ones(T,k+1)
    ξ == 0 ? (for p=2:k+1; o[p] = 0 end) :
             (for p=1:k; o[p+1] = o[p]*(ξ-p+1)/p end)

    return o

end
#...............................................................................
function bwd_interpolation_expansion_polynom1(ξ::T, k=3) where T<:Real

    o = Base.ones(T,k+1)
    ξ == 0 ? (for p=2:k+1; o[p] = 0 end) :
             (for p=1:k; o[p+1] = o[p]*(p-ξ-1)/p end)

    return o

end
#...............................................................................
@doc raw"""
    fdiff_interpolation_expansion_polynom(σ::T [, notation=bwd [; k=3]]) where T<:Real

Finite-difference expansion coefficient vector defining the ``k^{th}``-order
(default *third* order) Lagrange-polynomial interpolation of a tabulated
analytic function ``f[n]`` at offset ``σ`` with respect to index
position ``n``, which is positive for increasing index and negative for
decreasing index.

**Forward difference notation** (`notation = fwd`)

In this case we consider the tabulated interval ``f[n:n+k]``. The interpolated
value ``f[n-σ]`` is given by the forward-difference expansion

```math
f[n-σ] = \sum_{p=0}^k α_p(σ) Δ^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_polynom(σ, fwd; k=3)`](@ref)
`` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``. 

Application: This polynom can serve to predict `f[n-1]` by *extrapolation* (using ``σ=1``) 
if `f[n:n+k]` are known. More generally, it can serve to *interpolate* to (real) positions 
``n ≤ v ≤ n+k`` (using ``-k ≤ σ ≤ 0``) and predict `f[n-σ]` by *extrapolation* to (real) 
positions ``v<n`` (using ``σ > 0``) or ``v>n+k`` (using ``σ < -k``).  NB. The forward offset 
is given by ``σ ≡ n-v``.

**Backward difference notation** (`notation = bwd`)

In this case we consider the tabulated interval ``f[n-k:n]``. The interpolated
value ``f[n+σ]`` is given by the backward-difference expansion

```math
f[n+σ] = \sum_{p=0}^k β_p(σ) ∇^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_polynom(σ, bwd; k=3)`](@ref)
`` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``. 

Application: This polynom can serve to predict `f[n+1]` by *extrapolation* (using ``σ=1``) 
if `f[n-k:n]` are known. More generally, it can serve to *interpolate* to (real) positions 
``n-k ≤ v ≤ n`` (using ``-k ≤ σ ≤ 0``) and predict `f[n+σ]` by *extrapolation* to (real) 
positions ``v<n`` (using ``σ > 0``) or ``v>n+k`` (using ``σ < -k``). NB. The backward offset 
is given by ``σ ≡ -(n-v)``.

#### Example 1 - extrapolation on a uniform grid:
```
julia> σ = 1;

julia> f = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100];

julia> α = fdiff_interpolation_expansion_polynom(σ, fwd; k=5); println("α = $α")
α = [1, -1, 1, -1, 1, -1]

julia> Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
Fk = [6, -15, 20, -15, 6, -1]

julia> Fk ⋅ f[5:10] == f[4] == 16
true

julia> β = fdiff_interpolation_expansion_polynom(σ, bwd; k=5); println("β = $β")
β = [1, 1, 1, 1, 1, 1]

julia> revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
revBk = [-1, 6, -15, 20, -15, 6]

julia> revBk ⋅ f[1:6] == f[7] == 49
true

```
#### Example 2 - interpolation on a uniform grid:
```
julia> f = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

julia> v = 7.5; n = 5; k = 5;

julia> σ = n-v # forward difference offset
-2.5

julia> α = fdiff_interpolation_expansion_polynom(σ, fwd; k); println("α = $α")
α = [1.0, 2.5, 1.875, 0.3125, -0.0390625, 0.01171875]

julia> Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
Fk = [0.01171875, -0.09765625, 0.5859375, 0.5859375, -0.09765625, 0.01171875]

julia> Fk ⋅ f[n:n+k] ≈ 7.5
true

julia> v = 7.5; n = 9; k = 5;

julia> σ = -(n-v) # backward difference offset
-1.5

julia> β = fdiff_interpolation_expansion_polynom(σ, bwd; k); println("β = $β")
β = [1.0, -1.5, 0.375, 0.0625, 0.0234375, 0.01171875]

julia> revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
revBk = [-0.01171875, 0.08203125, -0.2734375, 0.8203125, 0.41015625, -0.02734375]

julia> revBk ⋅ f[n-k:n] ≈ v
true
```
"""
function fdiff_interpolation_expansion_polynom(σ::T, notation=bwd; k=3) where T<:Real

    o = CamiMath.isforward(notation) ? fwd_interpolation_expansion_polynom(σ; k) :
                                       bwd_interpolation_expansion_polynom(σ; k)
    return o

end

# ------------------------------------------------------------------------------
#                  fdiff_interpolation_expansion_weights(polynom)
# ------------------------------------------------------------------------------

function fwd_interpolation_expansion_weights(σ::T, ordering=reg; k=3) where T<:Real

    α = fdiff_interpolation_expansion_polynom(σ, fwd; k)
    o = fdiff_expansion_weights(α, fwd, ordering)

    return o

end
#...............................................................................
function bwd_interpolation_expansion_weights(σ::T, ordering=rev; k=3) where T<:Real

    β = fdiff_interpolation_expansion_polynom(σ, bwd; k)
    o = fdiff_expansion_weights(β, bwd, ordering)

    return o

end
#...............................................................................
@doc raw"""
    fdiff_interpolation_expansion_weights(σ::T, [, notation=bwd [, ordering=rev [, k=3]]]) where T<:Real
    fdiff_interpolation_expansion_weights(polynom [, notation=bwd [, ordering=rev]])

Finite-difference expansion weights vector defining the ``k^{th}``-order
(default *third* order) Lagrange-polynomial interpolation of a tabulated
analytic function ``f[n]`` at offset ``σ`` with respect to index
position ``n``, which is positive for increasing index and negative for
decreasing index.

**Forward difference notation** (`notation = fwd`)

In this case we consider the tabulated interval ``f[n:n+k]``. The interpolated
value ``f[n-σ]`` is given by the forward-difference expansion

```math
f[n-σ] = \sum_{p=0}^k α_p(σ) Δ^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_polynom(σ, fwd; k=3)`](@ref)
`` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``. In this notation the range
``-k ≤ σ ≤ 1`` corresponds to interpolation and the ranges ``σ < -k`` and
``σ > 1k`` to extrapolation.

**Backward difference notation** (`notation = bwd`)

In this case we consider the tabulated interval ``f[n-k:n]``. The interpolated
value ``f[n+σ]`` is given by the backward-difference expansion

```math
f[n+σ] = \sum_{p=0}^k β_p(σ) ∇^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_polynom(σ, bwd; k=3)`](@ref)
`` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``. In this notation the range
``-k ≤ σ ≤ 1`` corresponds to interpolation and the ranges ``σ < -k`` and
``σ > 1k`` to extrapolation.

#### Examples:
```
julia> Fk1 = fdiff_interpolation_expansion_weights(1, fwd, reg; k=4); println("Fk1 = $(Fk1)")
Fk1 = [5, -10, 10, -5, 1]

julia> revBk1 = fdiff_interpolation_expansion_weights(1, bwd, rev; k=4); println("revBk1 = $(revBk1)")
revBk1 = [1, -5, 10, -10, 5]

julia> α = fdiff_interpolation_expansion_polynom(1, fwd; k=4); println("α = $α")
α = [1, -1, 1, -1, 1]

julia> Fk1 = fdiff_interpolation_expansion_weights(α, fwd, reg); println("Fk1 = $(Fk1)")
Fk1 = [5, -10, 10, -5, 1]

julia> β = fdiff_interpolation_expansion_polynom(1, bwd; k=4); println("β = $β")
β = [1, 1, 1, 1, 1]

julia> revBk1 = fdiff_interpolation_expansion_weights(β, bwd, rev); println("revBk1 = $(revBk1)")
revBk1 = [1, -5, 10, -10, 5]
```
"""
function fdiff_interpolation_expansion_weights(σ::T, notation=bwd, ordering=rev; k=3) where T<:Real

    o = CamiMath.isforward(notation) ? fwd_interpolation_expansion_weights(σ, ordering; k) :
                                       bwd_interpolation_expansion_weights(σ, ordering; k)
    return o

end
function fdiff_interpolation_expansion_weights(polynom, notation=bwd, ordering=rev)

    if CamiMath.isforward(notation)

        o = CamiMath.isregular(ordering) ? reg_fwd_expansion_weights(polynom) :
                                           rev_fwd_expansion_weights(polynom)
    else

        o = CamiMath.isregular(ordering) ? reg_bwd_expansion_weights(polynom) :
                                           rev_bwd_expansion_weights(polynom)
    end

    return o

end

# ------------------------------------------------------------------------------
#                          fdiff_interpolation(f, x; k=3)
# ------------------------------------------------------------------------------

@doc raw"""
    fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}

Finite difference lagrangian interpolation (by default *third* order) at real
position v (in index units) with respect to the elements of the uniformly
tabulated analytic function `f[1:N]`. The interpolation points are situated on a
Lagrange polynomial of degree ``k`` (by default *third* degree) running through
``k+1`` subsequenct points of the tabulated function. Outside the tabulated
range, the method represents extrapolation on the lagrangian polynomial defined
by the first/last ``k+1`` tabulated points.

Beware that the interpolation becomes inaccurate if the tabulated function
cannot be approximated by a polynomial of degree ``k``.
#### Examples:
```
f = [1,2,3,4,5,6,7]
[fdiff_interpolation(f, v; k=3) for v=1:0.5:7]
  [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]

f = [1,4,9,16,25,36,49]
[fdiff_interpolation(f, v; k=3) for v=1:0.5:7]
 [1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0, 30.25, 36.0, 42.25, 49.0]

 f = [x^3 for x=-4:2]
 f1(v) = fdiff_interpolation(f, v; k=1)
 f2(v) = fdiff_interpolation(f, v; k=2)
 f3(v) = fdiff_interpolation(f, v; k=3)
 [[f1(v),f2(v),f3(v)] for v=1:0.5:9]
   [[-64.0, -64.0, -64.0], [-45.5, -43.25, -42.875], [-27.0, -27.0, -27.0],
   [-17.5, -16.0, -15.625], [-8.0, -8.0, -8.0], [-4.5, -3.75, -3.375],
   [-1.0, -1.0, -1.0], [-0.5, -0.5, -0.125], [0.0, 0.0, 0.0],
   [0.5, -0.25, 0.125], [1.0, 1.0, 1.0], [4.5, 3.75, 3.375], [8.0, 8.0, 8.0],
   [11.5, 13.75, 15.625], [15.0, 21.0, 27.0], [18.5, 29.75, 42.875],
   [22.0, 40.0, 64.0]]
```
The result for f3(v) is exact because the function is cubic and
the expansion is third order - see Figure below. The tabulated function is
given by the black points. The interpolation and extrapolation points are red.

![Image](../assets/lagrangian_interpolation.png)
"""
function fdiff_interpolation(f::Vector{T}, v::V; k=4) where {T<:Real, V<:Real}

    l = length(f)
    k = min(k,l-1)
    k > 0 || error("Error: k ≥ 1 required for lagrangian interpolation")
    n = v < 1 ? 1 : v < l-k ? floor(Int,v) : l-k
    α = fdiff_interpolation_expansion_polynom(n-v, k, fwd)
    o = fdiff_expansion(α, f[n:n+k], fwd)

    return o

end

# ------------------------------------------------------------------------------
#                fdiff_differentiation_expansion_polynom(ξ, k)
# ------------------------------------------------------------------------------

function fwd_differentiation_expansion_polynom(ξ::T, k=3) where T<:Real
# ==============================================================================
#   Forward difference expansion coefficients for differentiation of an
#   analytic function f(x) tabulated under the convention f[n,n+k] and
#   evaluated at the interpolation position n+ξ.
# ==============================================================================
    Float = (Float64, Float32, Float16, BigFloat) #######################################################################################

    a = T ∈ Float ? Base.ones(T,k+1) :
                    T <: Rational{}  ? Base.ones(T,k+1) :
                                       Base.ones(Rational{T},k+1)
    for i=1:k
        a[i+1] = iseven(i) ? -i : i
    end

    a = 1 ./ a; a[1] = 0

    ξ == 0 && return a
    α = fdiff_interpolation_expansion_polynom(-ξ, k, fwd)

    a,α = Base.promote(a,α)

return CamiMath.polynom_product_expansion(a, α, k)

end
#...............................................................................
function bwd_differentiation_expansion_polynom(ξ::T, k=5) where T<:Real
# ==============================================================================
#   Backward difference expansion coefficients for differentiation of an
#   analytic function f(x) tabulated under the convention f[n-k,n] and
#   evaluated at the interpolation position n+ξ.
# ==============================================================================
    Float = (Float64, Float32, Float16, BigFloat) #######################################################################################

    b = T ∈ Float ? Base.ones(T,k+1) :
                    T <: Rational{}  ? Base.ones(T,k+1) :
                                       Base.ones(Rational{T},k+1)
    for i=1:k
        b[i+1] = i
    end

    b = 1 ./ b; b[1] = 0

    ξ == 0 && return b
    β = fdiff_interpolation_expansion_polynom(ξ, k, bwd)

    b,β = Base.promote(b,β)

    return CamiMath.polynom_product_expansion(b, β, k)

end
#...............................................................................
@doc raw"""
    fdiff_differentiation_expansion_polynom(ξ::T [, k=3 [, notation=bwd]]) where T<:Real

Finite-difference expansion coefficient vector defining ``k^{th}``-order
*lagrangian differentiation* of the tabulated analytic function ``f[n]``
at offset ``ξ`` (with respect to index position ``n``), which is positive
for increasing index and negative for decreasing index.

**Forward difference notation** (`notation = fwd`)

```math
\frac{df}{dξ}[n+ξ]=\sum_{p=0}^kα_p(ξ)Δ^{p}f[n]
```

Offset convention: ``ξ = -σ`` with respect to index ``n`` in tabulated
interval ``f[n:n+k]``

**Backward difference notation** (`notation = bwd`)

```math
\frac{df}{dξ}[n+ξ]=\sum_{p=0}^kβ_p(ξ)∇^{p}f[n]
```
where ``β(ξ) ≡ [β_0(ξ),\ ⋯,\ β_p(ξ)]``

Offset convention: ``ξ = σ`` with respect to index ``n`` in tabulated
interval ``f[n-k:n]``
#### Example:
```
k = 2; ξ = 0
o = fdiff_differentiation_expansion_polynom(ξ, k); println(o)
 [0.0, 1.0, -1.5]
```
"""
function fdiff_differentiation_expansion_polynom(ξ::T, k=3, notation=bwd) where T<:Real

    o = CamiMath.isforward(notation) ? fwd_differentiation_expansion_polynom(ξ, k) :
                                       bwd_differentiation_expansion_polynom(ξ, k)
    return o

end

# ------------------------------------------------------------------------------
#                   create_lagrange_differentiation_matrix(k)
# ------------------------------------------------------------------------------

@doc raw"""
    create_lagrange_differentiation_matrix(k::Int)

Lagrange differentiation matrix, ``m[i,j]=s_{k-j}^k(i)``, for ``k^{th}``-order
lagrangian differentiation,
```math
\frac{dy}{dx}[i]= \sum_{j=0}^{k}m[i,j]y[j],
```
#### Example:
```
k = 3
create_lagrange_differentiation_matrix(k)
 4×4 Matrix{Rational{Int64}}:
  -11//6   3//1  -3//2   1//3
   -1//3  -1//2   1//1  -1//6
    1//6  -1//1   1//2   1//3
   -1//3   3//2  -3//1  11//6
```
"""
function create_lagrange_differentiation_matrix(k::Int)

    m = Matrix{Rational{Int}}(undef,k+1,k+1)

    for i=0:k
        polynom = CamiDiff.fdiff_differentiation_expansion_polynom(k-i,k)
        m[1+i,1:k+1] = fdiff_expansion_weights(polynom)
    end

    return m

end

# ------------------------------------------------------------------------------
#             trapezoidal_epw(k; rationalize=false, devisor=false)
# ------------------------------------------------------------------------------

@doc raw"""
    trapezoidal_epw(k::Int [; rationalize=false [, devisor=false]])

Endpoint weights vector ``a=[a_1,⋯\ a_k]`` of trapeziodal rule
optimized for functions of polynomial form,
```math
    ∫_0^n f(x) dx = a_1 (f_0+f_n) + ⋯ + a_k (f_{k-1}+f_{n-k+1})
                                                         + (f_k+⋯+f_{n-k}),
```
where ``k`` is *odd*. The rule is exact for polynonials of degree ``d=0,\ 1,
⋯,\ k-1``. For ``k=1`` the rule reduces to the ordinary trapezoidal rule.
By default the output is in Float64, optionally the output is rational, with or
without specification of the gcd devisor.
#### Example:
```
[trapezoidal_epw(k; rationalize=true, devisor=true) for k=1:2:9]
5-element Vector{Tuple{Int64, Int64, Vector{Int64}}}:
  (1, 2, [1])
  (3, 24, [9, 28, 23])
  (5, 1440, [475, 1902, 1104, 1586, 1413])
  (7, 120960, [36799, 176648, 54851, 177984, 89437, 130936, 119585])
  (9, 7257600, [2082753, 11532470, 261166, 16263486, -1020160, 12489922,
                                                     5095890, 7783754, 7200319])
```
"""
function trapezoidal_epw(k::Int; rationalize=false, devisor=false)
# ==============================================================================
# trapezoidal_epw(k; rationalize=false, devisor=false)
# ==============================================================================
    strWarn = "Warning: k = $(k) → $(k-1) (trapezoidal rule requires odd k)"
    Base.isodd(k) ? true : (k=k+1; println(strWarn))

    l = k - 1
    σ = Base.Matrix{Int}(undef,k,k)

    for i=0:k-1
        σ[i+1,:] = CamiMath.polynom_power([i,-1],l) # corresponds to a0,...ak
        σ[i+1,1] = σ[i+1,1] + i^l                  # correction of coeff a0
    end

    F = CamiMath.faulhaber_polynom(k)
    s = [CamiMath.polynom_power([-k,1],p) for p=0:k] .* F
    s[1][1] = -CamiMath.faulhaber_summation(k-1,l)//1

    c = [Base.sum([s[p+1][i+1] for p=i:k]) for i=0:k][1:end-1]

    σ = Base.inv(Base.transpose(σ))

    o = -σ * c  # solving matrix equation results in trapezoidal_epw

    if rationalize
        a = fdiff_adams_moulton_expansion_polynom(k)
        D = Base.denominator(Base.gcd(a))          # == Adams-Moulton devisor
        o = devisor ? (k, D, Base.round.(Int, o * D)) :
                      Base.round.(Int, o * D) // D  # convert to Rational{Int}
    end

    return o

end

# ------------------------------------------------------------------------------
#                 trapezoidal_integration(f, domain; k=5)
# ------------------------------------------------------------------------------

@doc raw"""
    trapezoidal_integration(f, x1, x2, weights)

Integral of the tabulated function ``f=[f_0,⋯\ f_n]`` over the `domain`
``x1 ≤ x ≤ x2`` using the optimized trapezoidal rule with endpoint correction
by the weights vector `weights`,
```math
    ∫_0^n f(x) dx = a_1 (f_0+f_n) + ⋯ + a_k (f_{k-1}+f_{n-k+1})
                                                         + (f_k+⋯+f_{n-k}).
```
The rule is exact for polynonials of degree ``d=0,\ 1,⋯\ k-1``.
For ``k=1`` the rule reduces to the ordinary trapezoidal rule (weights = [1/2]).
#### Examples::
```
p = 3
c = [1 for i=0:p]
pol = ImmutablePolynomial(c,:z)
Ipol = integrate(pol)
n = 10

x1=0.0
x2=5.0
x = collect(range(x1, x2, n))
f = pol.(x .-2.5)

w3 = trapezoidal_epw(3)
trapezoidal_integration(f, x1, x2, w3)
 15.416666666666673

Ipol(2.5)-Ipol(-2.5)
 15.41666666666666
```
"""
function trapezoidal_integration(f, x1, x2, weights)

    n = Base.length(f)
    k = Base.length(weights)
    s = (x2-x1)/(n-1)
    a = Base.ones(n); a[1:k] = weights; a[end-k+1:end] = Base.reverse(weights)
    #o = LinearAlgebra.dot(f, a) * s
    o = sum(f .* a) * s

    return o

end
