# Applications

Three basic operations are available for functions defined on a [`Grid`](@ref).

* interpolation of the tabulated 'function of interest' at a position between two subsequent points on the grid - see [`grid_interpolation(f, grid, r)`](@ref)

* integration of the tabulated 'function of interest' over the full (or part) of the grid - see [`grid_integration(f, grid)`](@ref)

* differentiation of the tabulated 'function of interest' on all points of the grid, on part of the grid, at a given point on the grid or at a position between two subsequent points on the grid - see [`grid_differentiation(f, grid)`](@ref)

```@docs
grid_interpolation(f::Vector{T}, grid::Grid{T}, rv::T, notation=fwd; k=5) where T<:Real
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real
grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real 
```