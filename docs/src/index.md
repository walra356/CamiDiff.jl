```@meta
CurrentModule = CamiDiff
```

# Home

`CamiDiff.jl` is a [Julia](http://julialang.org) package for one-dimensional finite-difference analysis. 

---

## Install

The package is installed using the Julia package manager

```
julia> using Pkg; Pkg.add("CamiDiff")

julia> using CamiDiff
```

## Introduction

`CamiDiff` has been developed to analyze mathematical functions, provided by the user 
*in tabulated form*. The package is based on finite-difference analysis methods in one dimension. 

Throughout the documentation ``f(x)`` will be the 'function of interest' under investigation. In tabulated 
form we write ``f[n]``. This form can be regarded as the result of the discretization of ``f(x)`` 
on a [`Grid`](@ref) of ``N`` points, addressable by the *grid index* ``n = 1, ⋯ N``. The [`Grid`](@ref) 
can be linear (uniform) or non-linear as specified by a [`gridfunction`](@ref) - see [Discretization](@ref).

The current implementation of `CamiDiff` was developed for *real functions of a single variable*. 
A set of four predefined [`gridfunction`](@ref)`s` is included in the package: 'exponential', 
'truncated-exponential', 'linear' and 'polynomial'. These are restricted to the domain ``[0, ∞)``. 
When applicable, we shall underline this restriction to the non-negative domain by using the variable ``r`` 
rather than ``x``, writing ``f(r)`` rather than ``f(x)``, with the implicit condition ``r ≥ 0``.

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

## Grid

The [`Grid`](@ref) object is the backbone for numerical procedures on the real domain ``[0, ∞)``. Its principal fields 
are `grid.r`, `grid.r′` and `grid.r′′`. These are discrete functions of `N` elements representing the grid function 
and its first two derivatives. The function ``f[n]`` is tabulated on this [`Grid`](@ref) and the function 
``r[n]``` = grid.r` represents the transformation by the [`gridfunction`](@ref). 

Once the [`Grid`](@ref) is specified, three [`Grid`](@ref) operations are at our disposal - see [Documentation](@ref)

[`grid_interpolation(f, rval, grid)`](@ref)

[`grid_differentiation(f, grid)`](#ref)

[`grid_integration(f, grid)`](@ref)

```@docs
Grid{T}
castGrid(ID::Int, N::Int, T::Type; h=1, rmax=10, p=5, polynom=[0,1], epn=5, k=5, msg=false)
gridfunction(ID::Int, n::Int, T::Type; h=1, p=5, polynom=[0,1], deriv=0)
gridtypename(ID::Int)
gridtypeID(name::String)
gridPos(rval::T, grid::Grid{T}) where T<:Real
fracPos(n::Int, rval::T, grid::Grid{T}; ϵ = 1e-8, k = 7) where T<:Real
```