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

The weights vector ``F^k ≡ [F_k^k,⋯\ F_0^k]`` corresponds to the expansion coefficient 
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

The weights vector ``\bar{B}^{k} ≡ [B_k^k,⋯\ B_0^k]`` corresponds to the expansion coefficient 
vector ``β ≡ [β_0,⋯\ β_k]`` of the ``k^{th}``-order *backward-difference* expansion (`polynom = β`).

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
positions ``v<n`` (using ``σ > 0``) or ``v>n+k`` (using ``σ < -k``).  

NB. The forward offset is defined as ``σ ≡ n-v``.

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
positions ``xv<n`` (using ``σ > 0``) or ``v>n+k`` (using ``σ < -k``). 

NB. The backward offset is defined as ``σ ≡ -(n-v)``.

#### Example:
```
julia> σ = 1; # offset correponding to extrapolation

julia> α = fdiff_interpolation_expansion_polynom(σ, fwd; k=5); println("α = $α")
α = [1, -1, 1, -1, 1, -1]

julia> Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
Fk = [6, -15, 20, -15, 6, -1]

julia> β = fdiff_interpolation_expansion_polynom(σ, bwd; k=5); println("β = $β")
β = [1, 1, 1, 1, 1, 1]

julia> revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
revBk = [-1, 6, -15, 20, -15, 6]
```
"""
function fdiff_interpolation_expansion_polynom(σ::T, notation=bwd; k=3) where T<:Real

    o = CamiMath.isforward(notation) ? fwd_interpolation_expansion_polynom(σ; k) :
                                       bwd_interpolation_expansion_polynom(σ; k)
    return o

end


# ------------------------------------------------------------------------------
#                fdiff_differentiation_expansion_polynom(σ, k)
# ------------------------------------------------------------------------------

function fwd_differentiation_expansion_polynom(σ::T; k=5) where T<:Real
# ==============================================================================
#   Forward difference expansion coefficients for differentiation of an
#   analytic function f(x) tabulated under the convention f[n,n+k] and
#   evaluated at the interpolation position n+σ.
# ==============================================================================
    Float = (Float64, Float32, Float16, BigFloat) #######################################################################################

    a = T ∈ Float ? Base.ones(T,k+1) :
                    T <: Rational{}  ? Base.ones(T,k+1) :
                                       Base.ones(Rational{T},k+1)
    for i=1:k
        a[i+1] = iseven(i) ? -i : i
    end

    a = 1 ./ a; a[1] = 0

    σ == 0 && return a
    α = fdiff_interpolation_expansion_polynom(σ, fwd; k)

    a,α = Base.promote(a,α)

return CamiMath.polynom_product_expansion(a, α, k)

end
#...............................................................................
function bwd_differentiation_expansion_polynom(α::T; k=5) where T<:Real
# ==============================================================================
#   Backward difference expansion coefficients for differentiation of an
#   analytic function f(x) tabulated under the convention f[n-k,n] and
#   evaluated at the interpolation position n+σ.
# ==============================================================================
    Float = (Float64, Float32, Float16, BigFloat) #######################################################################################

    b = T ∈ Float ? Base.ones(T,k+1) :
                    T <: Rational{}  ? Base.ones(T,k+1) :
                                       Base.ones(Rational{T},k+1)
    for i=1:k
        b[i+1] = i
    end

    b = 1 ./ b; b[1] = 0

    α == 0 && return b
    β = fdiff_interpolation_expansion_polynom(α, bwd; k)

    b,β = Base.promote(b,β)

    return CamiMath.polynom_product_expansion(b, β, k)

end
#...............................................................................
@doc raw"""
    fdiff_differentiation_expansion_polynom(σ::T [, notation=bwd [; k=5]]) where T<:Real

Finite-difference expansion coefficient vector defining ``k^{th}``-order
*lagrangian differentiation* of the tabulated analytic function ``f[n]``
at position  ``v=n-σ``.

**Forward difference notation** (`notation = fwd`)

```math
\frac{df}{dσ}[n+σ]=\sum_{p=0}^kα_p(σ)Δ^{p}f[n]
```

Offset convention: ``σ = n-v`` with respect to index ``n`` in tabulated
interval ``f[n:n+k]``

**Backward difference notation** (`notation = bwd`)

```math
\frac{df}{dσ}[n+σ]=\sum_{p=0}^kβ_p(σ)∇^{p}f[n]
```
where ``β(σ) ≡ [β_0(σ),\ ⋯,\ β_p(σ)]``

Offset convention: ``σ = -(n-v)`` with respect to index ``n`` in tabulated
interval ``f[n-k:n]``
#### Example:
```
julia> σ = 0; # offset correponding to differentiation

julia> o = fdiff_differentiation_expansion_polynom(0, fwd; k=5); println(o)
Rational{Int64}[0, 1, -1//2, 1//3, -1//4, 1//5]
```
"""
function fdiff_differentiation_expansion_polynom(σ::T, notation=bwd; k=5) where T<:Real

    o = CamiMath.isforward(notation) ? fwd_differentiation_expansion_polynom(σ; k) :
                                       bwd_differentiation_expansion_polynom(σ; k)
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
        polynom = CamiDiff.fdiff_differentiation_expansion_polynom(i-k, bwd; k)
        m[1+i,1:k+1] = fdiff_expansion_weights(polynom, bwd, rev)
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
