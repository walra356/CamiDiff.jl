# Documentation

## Discretization

Mathematically, the discretization is based on the map ``n ↦ x,`` which defines the discrete function

```math
x[n] = s_0 * g(t[n]).
```

Here ``g(t)`` is called the [`gridfunction`](@ref) and ``s_0`` the *scaling factor*. The [`gridfunction`](@ref) 
is defined as a (generally nonlinear) function *running through the origin*: ``g(0) = 0``. Its argument is 
the *ticks function*

```math
t[n] ≡ (n−u) * h,
```

which is a discrete *linear* function, where ``u`` is called the *index base* and ``h`` the *step size*. Writing

```math
f[n] = f(x[n]),
```

we recognize in ``f[n]`` a discrete function representing the function ``f(x)`` at position ``x[n]``. This represents 
the tabulated function to be provided by the user. 

NB. The discrete function ``f[n]`` is defined on the *grid of natural numbers*, a uniform grid with unit step size. 
This uniformity greatly simplifies the numerical analysis. The stepsize of the *ticks function*, ``h``, determines 
the *coarseness* of the grid. The results of a finite-difference calculation on a coarse grid will be less accurate 
than those on a fine grid, but the algorithm is identical, because the relevant finite-difference expansions only 
depend *implicitly* on ``h``. Since [Julia](http://julialang.org) uses unit-based indexing (``u = 1``), the index 
convention implies ``f[1] = f(0)``.  

## Finite differences

Let ``f(x)`` be a real, regular function of the variable ``x``. 
The *forward difference* between ``f(x+h)`` and ``f(x)`` is defined as 

```math
Δ f(x)=f(x+h)-f(x).
```

Here, ``Δ`` is called the *forward-difference operator*. Likewise, one defines backward differences 
with the *backward-difference operator* ``∇``, 

```math
∇ f(x)=f(x)-f(x-h).
```
**Forward difference notation**

We first focus on *forward differences*. The derivative of ``f(x)`` is given by 

```math
f^′(x)=\underset{h→0}{\mathrm{lim}}\,\frac{f(x+h)-f(x)}{h}=\underset{Δ x→0}{\mathrm{lim}}\,\frac{Δ f(x)}{Δ x},
```

where ``h ≡ Δx ≥ 0`` is the *difference interval*. Introducing the differential operator, ``f^′(x) ≡ Df(x)``, we have 

```math
D≡\frac{d}{dx}=\underset{Δ x→0}{\mathrm{lim}}\,\frac{Δ}{Δ x}=\underset{h→0}{\mathrm{lim}}\,\frac{Δ}{h}.
```
**Backward difference notation**

In terms of *backward differences*, the derivative of ``f(x)`` is given by 

```math
f^′(x)=\underset{h→0}{\mathrm{lim}}\,\frac{f(x)-f(x-h)}{h}=\underset{Δ x→0}{\mathrm{lim}}\,\frac{∇ f(x)}{Δ x},
```

where ``h ≡ Δx ≥ 0`` is the *difference interval*. Introducing the differential operator, ``f^′(x) ≡ Df(x)``, we have 

```math
D≡\frac{d}{dx}=\underset{Δ x→0}{\mathrm{lim}}\,\frac{∇}{Δ x}=\underset{h→0}{\mathrm{lim}}\,\frac{∇}{h}.
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
Since ``Δ`` depends implicitly on ``h``, the explicit dependence on ``h`` can be replaced by an implicit dependence on ``h`` through an expansion in powers of ``Δ`` ,
```math
f(x-h)=(1+Δ)^{-1}f(x)=(1-Δ+Δ^{2}-Δ^3+⋯)f(x).
```
By choosing the proper expansion order, ``f(x-h)`` can be approximated to any desired level of accuracy.

**Backward difference notation**

For *backward differences*, we rewrite the backward-difference definition in the form 
```math
f(x-h)=(1-∇)f(x),
```
where ``B≡(1-∇)`` is the *backward translation operator*. Comparing this *backward* translation with the Taylor expansion we obtain, by formal inversion of the operator ``B``, an operator identity for the *forward* translation operator ``T``, 
```math
B≡(1-∇)=e^{-hD}=T^{-1}\,\,\,⇒\,\,\,T=e^{hD}=(1+∇)^{-1}.
```
Note how the *backward* translation operator was identified with the inverse *forward* translation operator, ``B=T^{-1}``. As with forward differences, the explicit dependence on ``h`` can be replaced by an implicit dependence on ``h`` through an expansion in powers of ``∇``, 
```math
f(x+h)=(1-∇)^{-1}f(x)=(1+∇+∇^{2}+∇^3+⋯)f(x).
```
By choosing the proper expansion order, ``f(x+h)`` can be approximated to any desired level of accuracy. 

## Finite differences

Let ``f[n]`` be the 'function of interest', tabulated in *forward order*  
(growing index) on the *grid of natural numbers*.

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

---

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

Finite-difference calculus builds on *finite-difference expansions* on a uniform grid.

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
*growing* grid position). After substituting the forward difference, ``Δ = f[n+1] - f[n]``, 
the finite-difference expansion takes the form

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

---

Summary:

In `CamiDiff`, any finite-difference expansion is defined by a (user-supplied) `polynom` 
(the vector of expansion coefficients). In *fwd-difference* notation we write

`polynom` `` → α ≡ [α_0,⋯\ α_k]``.

Examples are given below. Once we have the *coefficients* (in the form of the `polynom`), 
we can calculate the *weights* (the *fwd-difference* weights vector) in reg-fwd-notation,

`weights =` [`fdiff_expansion_weights(polynom, fwd, reg)`](@ref) `` → F^k ≡ [F^k_0,⋯\ F^k_k]``

and the result of the expansion is obtained by evaluating the following inner product (in *fwd-difference* notation)  

```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n] = F^k \cdot f[n:n+k].
```

---

Examples:

The `polynom`s of four common expansions in *fwd-difference* notation are:

interpolation: [`fdiff_interpolation_expansion_polynom(σ, k, fwd)`](@ref) `` → α ≡ [α_0(σ),⋯\ α_k(σ)]``

differentiation: [`fdiff_differentiation_expansion_polynom(σ, k, fwd)`](@ref) `` → α ≡ [α_0(σ),⋯\ α_k(σ)]``

Adams-Bashford: [`fdiff_adams_bashford_expansion_polynom(k, fwd)`](@ref) `` → α ≡ [α_0,⋯\ α_k]``

Adams-Moulton: [`fdiff_adams_moulton_expansion_polynom(k, fwd)`](@ref) `` → α ≡ [α_0,⋯\ α_k]``

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

---

Summary:

In `CamiDiff`, any finite-difference expansion is defined by a (user-supplied) `polynom` 
(the vector of expansion coefficients). In *bwd-difference* notation we write

`polynom` `` → β ≡ [β_0,⋯\ β_k]``.

Examples are given below. Once we have the *coefficients* (in the form of `polynom`) we can calculate the *weights*
(the *bwd-difference* weights vector) in rev-bwd-notation,

`weights =` [`fdiff_expansion_weights(polynom, bwd, rev)`](@ref) `` → \bar{B}^k(σ) ≡ [B_k^k(σ),⋯\ B_0^k(σ)]``,

and the result of the expansion is obtained by evaluating the following inner product (in *bwd-difference* notation),

```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n] = \bar{B}^k(σ) \cdot f[n-k:n].
```

Examples:

The `polynom`s of four common expansions in *bwd-difference* notation are:

interpolation: [`fdiff_interpolation_expansion_polynom(σ, k, bwd)`](@ref) `` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``

differentiation: [`fdiff_differentiation_expansion_polynom(σ, k, bwd)`](@ref) `` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``

Adams-Bashford: [`fdiff_adams_bashford_expansion_polynom(k, bwd)`](@ref) `` → β(σ) ≡ [β_0,⋯\ β_k]``

Adams-Moulton: [`fdiff_adams_moulton_expansion_polynom(k, bwd)`](@ref)`` → β(σ) ≡ [β_0,⋯\ β_k]``

```@docs
# fdiff_expansion(polynom, f, notation=bwd)
fdiff_expansion_weights(polynom, notation=CamiMath.bwd, ordering=CamiMath.rev)
```

## Lagrange interpolation

The *Lagrange polynomial* of degree ``k`` on a uniform grid is the polynomial running 
through ``k+1`` subsequent points on the grid. We derive expressions for 
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
f[n-1] = (1 + Δ)^{-1} f[n] ≡ \sum_{p=0}^{\infty}(-1)^p Δ^p f[n]+⋯,
```
```math
f[n-2] = (1 + Δ)^{-2} f[n] ≡ \sum_{p=0}^{\infty}(-1)^p pΔ^p f[n],
```
```math
\vdots
```
where ``p`` is called the order of the expansion and ``n`` is the index of the reference
position. For interpolation position ``v=n-σ`` (where ``σ`` may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n-σ] = (1 + Δ)^{-σ} f[n] ≡ \sum_{p=0}^{\infty} (-1)^p\ l_p(σ) Δ^p f[n],
```
where ``α_p(σ) = (-1)^p\ l_p(σ)`` is the ``p^{th}``-order *finite-difference expansion coefficient*
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
``α_p ≡ α_p(1) ≡ (-1)^p``, which are the coefficients for the 'next-point' expansion. 

For ``-k ≤ σ ≤ 0`` the method can be used for *interpolation* over the grid position interval 
``n ≤ v ≤ n+k``. Outside this interval, in particular for ``σ > 0``, the method amounts to 
*extrapolation*. The method is most accurate for ``-1 ≤ σ ≤ 1`` (corresponding to the grid 
position interval ``n-1 ≤ v ≤ n+1``). Extrapolation to values ``v > n+k`` is not recommended. 

NB. The forward offset is defined as ``σ ≡ n-v``.

---

Summary:

In `CamiDiff`, the `polynom` of the *fwd-interpolation* expansion is calculated with

`polynom =` [`fdiff_interpolation_expansion_polynom(σ, k, fwd)`](@ref) `` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``,
where ``α_p(σ) = (-1)^p\ l_p(σ)``.

Once we have the *coefficients* (in the form of `polynom`) we can calculate the *weights* 
(the *fwd-interpolation* expansion weights vector) in reg-fwd-notation,

`weights =` [`fdiff_expansion_weights(polynom, fwd, reg)`](@ref) ``→ F^k(σ) ≡ [F_0^k(σ),⋯\ F_k^k(σ)]``,

and the *interpolated value* at grid position `n-σ` evaluates (in fwd-difference notation) to

```math
f[n-σ] = F^{k}(σ) \cdot f[n:n+k].
```

---

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
index. For interpolation position ``v=σ-n`` (where σ may be *real* valued in
index units) these expansions can be generalized to the form of
*lagrangian interpolation*,

```math
f[n+σ] = (1 - ∇)^{-σ} f[n] ≡ \sum_{p=0}^{\infty} l_p(σ) ∇^p f[n]⋯,
```
where ``β_p(σ) = l_p(σ)`` is the ``p^{th}``-order *finite-difference expansion coefficient* 
for lagrangian interpolation, with ``(σ)_{p}`` being the Pochhammer symbol `CamiMath.pochhammer`.  
Note that for ``σ = 1`` we find ``β_p ≡ β_p(1) ≡ 1``, which are the coefficients for the 
'next-point' expansion. 

For ``-k ≤ σ ≤ 0`` the method can be used for *interpolation* over the grid position interval 
``n-k ≤ v ≤ n``, outside this interval, in particular for ``σ > 0``, the method amounts to 
*extrapolation*. The method is most accurate for ``-1 ≤ σ ≤ 1`` (corresponding to the grid 
position interval ``n-1 ≤ v ≤ n+1``). Extrapolation to values ``v < n-k`` is not recommended. 

NB. The backward offset is defined as ``σ ≡ -(n-v)``.

---

Summary:

In `CamiDiff`, the `polynom` of the *bwd-interpolation* expansion is calculated with

`polynom =` [`fdiff_interpolation_expansion_polynom(σ, k, bwd)`](@ref) `` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``,
where ``β_p(σ) = l_p(σ)``.

Once we have the *coefficients* (in the form of `polynom`) we can calculate the *weights*
(the *bwd-interpolation* expansion weights vector) in rev-bwd-notation,

`weights =` [`fdiff_expansion_weights(polynom, bwd, rev)`](@ref) `` → \bar{B}^k(σ) ≡ [B_k^k(σ),⋯\ B_0^k(σ)]``,

and the *interpolated value* at grid position `n+σ` evaluates  (in *bwd-difference* notation) to
```math
f[n+σ] = \bar{B}^k(σ) \cdot f[n-k:n].
```
##### Example 1:

Demonstration of forward-difference *extrapolation* to 'next point' (grid position 'v=n-1')
```
julia> n=5; v=4; k=5;

julia> σ = n-v # fwd offset to 'next point'
1

julia> α = fdiff_interpolation_expansion_polynom(σ, k, fwd); println("α = $α")
α = [1, -1, 1, -1, 1, -1]

julia> Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
Fk = [6, -15, 20, -15, 6, -1]

julia> f = [v^2 for v=1:10]; println("f = $f")
f = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]

julia> Fk ⋅ f[n:n+k] == f[n-1] == v^2
true
```
##### Example 2:

Demonstration of backward-difference *interpolation* to grid position 'v=6.25`
```
julia> n=9; v=6.25; k=5;

julia> σ = -(n-v) # bwd-offset of interpolation position
-1.25

julia> β = fdiff_interpolation_expansion_polynom(σ, k, bwd); println("β = $β")
β = [1.0, -2.75, 2.40625, -0.6015625, -0.03759765625, -0.0093994140625]

julia> revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
revBk = [0.0093994140625, -0.0845947265625, 0.845947265625, 0.281982421875, -0.0604248046875, 0.0076904296875]

julia> f = [v^2 for v=1:10]; println("f = $f")
f = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]

julia> revBk ⋅ f[n-k:n] ≈ v^2
true
```

```@docs
fdiff_interpolation_expansion_polynom(σ::T, k=3, notation=bwd) where T<:Real
# fdiff_interpolation_expansion_weights(σ::T, notation=bwd, ordering=rev; k=3) where T<:Real
# fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}
```

## Lagrangian differentiation

**Forward difference notation** (`notation = fwd`)

To derive the *lagrangian differentiation* formulas we formally differentiate

```math
f[n-σ] = (1+Δ)^{-σ} f[n]
```
with respect to ``-σ``,

```math
-\frac{df}{dσ}[n-σ]
=ln(1+Δ)\ (1+Δ)^{-σ}f[n]
=\sum_{q=1}^{k}(-1)^q\tfrac{1}{q}Δ^{q}\sum_{p=0}^{k}(-1)^p\ l_{p}(σ)Δ^{p}f[n]+⋯,
```

where ``l_p(σ) ≡ (σ)_p/p!``, with ``(σ)_{p}`` being the Pochhammer symbol `CamiMath.pochhammer`.  
Rewriting the r.h.s. as a single expansion in powers of ``Δ``, we obtain to order ``k``

```math
-\frac{df}{dσ}[n-σ]=\sum_{p=1}^{k}α_p(σ)Δ^{p}f[n]=F^k(σ) ⋅ f[n:n+k],
```

where ``α_p(σ)`` represents the *finite-difference expansion coefficients*
for *lagrangian differentiation* at position ``n-σ``. 

In the general case ``(σ ≠ 0)``, the coefficients ``α_p(σ)`` are obtained numerically 
by polynomial multiplication using the function
[`CamiMath.polynom_product(p_a, p_b)`](@extref CamiMath.polynom_product), where

```math
\begin{aligned}
p_a &= [0, -1, 1/2, ⋯ (-1)^k\ 1/k] \\
p_b &= [l_0(σ), -l_1(σ), ⋯ (-1)^k l_k(σ)].
\end{aligned}
```
Special case:

In the special case ``σ = 0``, we have ``p_b = [ 1, 0, ⋯, 0]`` and the *fwd-difference* 
coefficient vector for *lagrangian differentiation* at position ``n`` reduces to

```math
α^k = α^k(0) = p_a = [0, -1, 1/2, ⋯ (-1)^k\ 1/k].
```

---

Summary:
 
In `CamiDiff`, the `polynom` of the *forward-differentiation* expansion is calculated by

`polynom =` [`fdiff_differentiation_expansion_polynom(σ, k, fwd)`](@ref) `` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``, 
with ``α_0(σ)≡ 0``.

Once we have the *coefficients* (in the form of `polynom`) we can calculate the *weights* 
(the *fwd-differentiation* weights vector) in reg-fwd-notation,

`weights =` [`fdiff_expansion_weights(polynom, fwd, reg)`](@ref) `` → F^k(σ) ≡ [F^k_0(σ),⋯\ F^k_k(σ)]``,

and the *derivative* at grid position `n-σ` in fwd-difference notation evaluates to

```math
-\frac{df}{dσ}[n-σ]
=\sum_{j=0}^{k}F_j^k(σ)f[n+j]
= F^k(σ) ⋅ f[n:n+k],
```

##### Example:
First derivative of the tabulated function ``f[n:n+k]`` at the position ``v = 6.5`` (zero-offset)
```
n=5; v=6.5; k=5;

julia> σ = n-v # forward offset at grid position ``v``.
-1.5

julia> α = fdiff_differentiation_expansion_polynom(σ, k, fwd); println("α = $α")
α = [0.0, 1.0, 1.0, -0.041666666666666685, 0.0, 0.004687500000000011]

julia> Fk = fdiff_expansion_weights(α, fwd, reg); println("Fk = $(Fk)")
Fk = [0.036979166666666674, -1.1015625, 1.078125, 0.005208333333333426, -0.023437500000000056, 0.004687500000000011]

julia> f = [v^2 for v=1:10]; println("f = $f")
f = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]

julia> Fk ⋅ f[n:n+k] ≈ 2v
true

julia> f = [v^3 for v=1:10]; println("f = $f")
f = [1, 8, 27, 64, 125, 216, 343, 512, 729, 1000]

julia> Fk ⋅ f[n:n+k] ≈ 3v^2
true
```

---

**backward difference notation** (`notation = bwd`)

To derive the *lagrangian differentiation* formulas we formally differentiate

```math
f[n+σ] = (1 - ∇)^{-σ} f[n]
```
with respect to ``σ``,

```math
\frac{df}{dσ}[n+σ]
=-ln(1-∇)\ (1-∇)^{-σ}f[n]
=\sum_{q=1}^{k}\tfrac{1}{q}∇^{q}\sum_{p=0}^{k}l_{p}(σ)∇^{p}f[n]+⋯,
```

where ``l_p(σ) ≡ (σ)_p/p!``, with ``(σ)_{p}`` being the Pochhammer symbol `CamiMath.pochhammer`.
Rewriting the r.h.s. as a single expansion in powers of ``∇``, we obtain to order ``k``

```math
\frac{df}{dσ}[n+σ]=\sum_{p=1}^{k}β_p(σ)∇^{p}f[n]=\bar{B}^k(σ) ⋅ f[n-k:n],
```

where ``β_p(σ)`` represents the *finite-difference expansion coefficients*
for *lagrangian differentiation* at position ``n+σ``. 

In the general case ``(σ ≠ 0)``, the coefficients ``β_p(σ)`` 
are obtained numerically by polynomial multiplication using the function
[`CamiMath.polynom_product(p_a, p_b)`](@extref CamiMath.polynom_product), where 

```math
\begin{aligned}
p_a &= [0, 1, 1/2, ⋯ 1/k]\\
p_b &= [l_0(σ), l_1(σ), ⋯ l_k(σ)].
\end{aligned}
```
Special case:

In the special case ``σ = 0``, we have ``p_b = [ 1, 0, ⋯, 0]`` and the *bwd-difference* 
coefficient vector for *lagrangian differentiation* at position ``n`` reduces to

```math
β^k = β^k(0) = p_a = [0, 1, 1/2, ⋯  1/k].
```

---

Summary:

In `CamiDiff`, the `polynom` of the *backward-differentiation* expansion is calculated by

`polynom =` [`fdiff_differentiation_expansion_polynom(σ, k, fwd)`](@ref) `` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``, 
with ``β_0(σ)≡ 0``.

Once we have the *coefficients* (in the form of `polynom`) we can calculate the *weights* 
(the *bwd-differentiation* weights vector) in rev-bwd-notation,

``\bar{B}^k(σ) =`` [`fdiff_expansion_weights(polynom, bwd, rev)`](@ref) `` → \bar{B}^k(σ) ≡ [B_k^k(σ),⋯\ B_0^k(σ)]``,

and the *derivative* at grid position `n+σ` in *bwd-difference* notation evaluates to 

```math
\frac{df}{dσ}[n+σ]=\bar{B}^k(σ) ⋅ f[n-k:n],
```

##### Example:
First derivative of the tabulated function ``f[n-k:n]`` at the position ``v = n`` (zero-offset)
```
julia> n=9; v=9; k=5;

julia> σ = -(n-v) # backward offset at grid position ``v``.
0

julia> β = fdiff_differentiation_expansion_polynom(σ, k, bwd); println("β = $β")
β = Rational{Int64}[0, 1, 1//2, 1//3, 1//4, 1//5]

julia> revBk = fdiff_expansion_weights(β, bwd, rev); println("revBk = $(revBk)")
revBk = Rational{Int64}[-1//5, 5//4, -10//3, 5, -5, 137//60]

julia> f = [v^2 for v=1:10]; println("f = $f")
f = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]

julia> revBk ⋅ f[n-k:n] ≈ 2v
true
```

```@docs
fdiff_differentiation_expansion_polynom(σ::T, k=5, notation=bwd) where T<:Real
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
