# Grid

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
