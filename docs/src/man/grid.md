# Grid

The [`Grid`](@ref) object is the backbone for numerical procedures on the real domain ``[0, ∞)``. Its principal fields 
are `grid.r`, `grid.r′` and `grid.r′′`. These are discrete functions of `N` elements representing the grid function 
and its first two derivatives. The function ``f[n]`` is tabulated on this [`Grid`](@ref) and the function 
``r[n]``` = grid.r` represents the transformation by the [`gridfunction`](@ref). 

Once the [`Grid`](@ref) is specified, three basic [`Grid`](@ref) operations are at our disposal: for *interpolation*, *integration* and *differentiation of functions tabulated on this grid - see [Applications](@ref)

```@docs
Grid{T}
castGrid(ID::Int, N::Int, T::Type; h=1, rmax=10, p=5, polynom=[0,1], epn=5, k=5, msg=false)
```

## Grid functions

```@docs
gridfunction(ID::Int, n::Int, T::Type; h=1, p=5, polynom=[0,1], deriv=0)
```

## Grid navigation

```@docs
gridtypename(ID::Int)
gridtypeID(name::String)
gridPos(rval::T, grid::Grid{T}) where T<:Real
fracPos(n::Int, rval::T, grid::Grid{T}; ϵ = 1e-8, k = 7) where T<:Real
```

# Applications

Three basic operations are available for functions defined on a [`Grid`](@ref).

* *interpolation* of the tabulated function `f` at position `r` between two subsequent points on the [`Grid`](@ref) - see [`grid_interpolation(f, grid, r)`](@ref)

* *integration* of the tabulated function `f` over the full (or part) of the [`Grid`](@ref) - see [`grid_integration(f, grid)`](@ref)

* *differentiation* of the tabulated function `f` on all points of the [`Grid`](@ref), on part of the [`Grid`](@ref), at a given point on the [`Grid`](@ref) or at a position between two subsequent points on the [`Grid`](@ref) - see [`grid_differentiation(f, grid)`](@ref)

```@docs
grid_interpolation(f::Vector{T}, grid::Grid{T}, rv::T, notation=fwd; k=5) where T<:Real
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real
grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real 
```