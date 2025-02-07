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