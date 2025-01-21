```@meta
CurrentModule = CamiDiff
```

# Home

`CamiDiff.jl` is a [Julia](http://julialang.org) package for finite-difference analysis of *real analytic functions of a single variable*. 

---

## Install

The package is installed using the Julia package manager

```
julia> using Pkg; Pkg.add("CamiDiff")

julia> using CamiDiff
```

# Introduction

`CamiDiff` has been developped to study *continuously differentiable* functions, provided by the user in *tabulated form*. 
The package is based on the methods of finite-difference analysis. 

Throughout the documentation ``f(x)`` will be the function under investigation. The tabulated form can be regarded as 
the result of discretization of ``f(x)``, onto a [`Grid`](@ref) of ``N`` points, addressable by the *gridindex* ``n = 1, ⋯ N``.
The [`Grid`](@ref) can be linear or non-linear as specified by a [`gridfunction`](@ref) - see [Discretization](@ref).

The current implementation of `CamiDiff` was developped for *real functions of a single variable*, using 
[`gridfunction`](@ref)s defined on the domain ``[0, ∞)``. A set of 4 predefined types is included: exponential, 
quasi-exponential, linear and polynomial. To underline the restriction to the non-negative domain, we shall often 
use the variable ``r`` rather than ``x``, writing ``f(r)`` rather than ``f(x)``, with the implicit condition ``r ≥ 0``.

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

we recognize in ``f[n]`` a discrete function representing the function ``f(x)`` at position ``x[n]``. 

Note that ``h`` determines the *coarseness* of the [`Grid`](@ref). The results of a finite-difference calculation 
on a coarse grid will be less accurate than those on a fine grid, but the algorithm is identical, because the relevant
finite-difference expansions only depend on ``h`` *implicitely*. Since [Julia](http://julialang.org) uses 
unit-based indexing (``u = 1``), the index convention implies ``f[1] = f(0)``.  

## Grid

The [`Grid`](@ref) object is the backbone for the numerical procedure on a (generally) non-uniform
grid. Its principal fields are `grid.r`, `grid.r′` and `grid.r′′` which are discrete
functions of `N` elements representing the grid function and its first two derivatives.

```@docs
Grid{T}
castGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001,  p=5, polynom=[0,1], epn=5, k=7, msg=true)
gridfunction(ID::Int, n::Int, h::T; p=5, polynom=[0,1], deriv=0) where T <: Real
gridname(ID::Int)
findIndex(rval::T, grid::Grid{T}) where T<:Real
findΔn(n::Int, rval::T, grid::Grid{T}; ϵ = 1e-8, k = 7) where T<:Real
grid_interpolation(f::Vector{T}, rval::T, grid::Grid{T}; k=5) where T<:Real
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real
grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
```

## Finite differences

Having discretized the *analytic* function ``f(x)``, one defines finite differences. 
The *forward difference* of ``f(x+h)`` and ``f(x)`` is defined by 

```math
Δ f(x)=f(x+h)-f(x).
```

Here, ``Δ`` is called the *forward-difference operator*. Likewise one defines backward differences 
with the *backward-difference operator* ``∇``, 

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