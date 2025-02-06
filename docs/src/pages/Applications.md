# Applications

The three elementary operations are:

[`grid_interpolation(f, rval, grid)`](@ref)

[`grid_differentiation(f, grid)`](@ref)

[`grid_integration(f, grid)`](@ref)

```@docs
grid_interpolation(f::Vector{T}, grid::Grid{T}, rv::T, notation=fwd; k=5) where T<:Real
grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real
grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
```