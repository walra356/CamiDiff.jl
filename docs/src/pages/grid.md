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
