# SPDX-License-Identifier: MIT

# author: Jook Walraven - 18-11-2024

# ==============================================================================
#                               grid.jl
# ============================================================================== 

# ------------------------------------------------------------------------------
#           Grid (ID, name, Type, N, r, r′, r′′, h, r0, epn, epw, k)
# ------------------------------------------------------------------------------

"""
Grid(ID, name, T, N, r, r′, h, r0, epn, epw, k)

Type with fields:
* `.ID`:    grid identifer name (`::Int`)
* `.name`:  grid identifer name (`::String`)
* `.T`:     gridType (`::Type`)
* `.N`:     number of grid points (`::Int`)
* `.r  `:   tabulated grid function (`::Vector{T}`)
* `.r′ `:   tabulated first derivative of grid function (`::Vector{T}`)
* `.r′′`:   tabulated second derivative of grid function (`::Vector{T}`)
* `.h` :    grid step multiplyer (`::T`)
* `.r0`:    grid scale factor (`::T`)
* `.epn`:   number of endpoints used for trapezoidal endpoint correction (must be odd) (`::Int`)
* `.epw`:   trapezoidal endpoint weights for n=1:epn (`::Vector{Vector{T}}`)
* `.k`:     finite-difference order (`::Int`)
The object `Grid` is best created with the function [`castGrid`](@ref).
"""
struct Grid{T}
ID::Int
name::String
T::Type
N::Int
r ::Vector{T}
r′::Vector{T}
r′′::Vector{T}
h::T
r0::T
epn::Int
epw::Vector{Vector{T}}
k::Int
end

# ------------------------------------------------------------------------------
#                        gridfunction(n, h; deriv=0)
# ------------------------------------------------------------------------------

function _walterjohnson(n::Int, h::T; deriv=0) where T <: Real
# ==============================================================================
#  gridfunction(n, h) = (exp((n-1) * h)-1.0) # gridfunction from Walter Johnson
# ==============================================================================
    deriv ≥ 0 || return 0.0
    
    f = deriv > 0 ? h^(deriv)*exp(n*h) : exp(n*h)-1
    
    return f
    
end
# ..............................................................................    
function _jw_gridfunction(n::Int, h::T; p=5, deriv=0) where T <: Real
# ==============================================================================
# jw_gridfunction(n, h [; p=5[, deriv=0]]) based on truncated exponential 
# ==============================================================================
    deriv ≥ 0 || return T(0)
    deriv ≤ p || return T(0)
    
    nul = T(0)
    
    f = deriv > 0 ? h^(deriv)*CamiMath.texp(n*h, nul, p-deriv) : 
                    CamiMath.texp(n*h, nul, p) - 1  # note: texp() not exp()
    
return f
    
end
# ..............................................................................     
function _linear_gridfunction(n::Int, h::T; deriv=0) where T <: Real
# ==============================================================================
#  linear_gridfunction(n, h; deriv) = n * h
# ==============================================================================
    deriv ≥ 0 || return T(0)
    deriv ≤ 1 || return T(0)

    f = deriv > 0 ? h : h * n
    
    return f
    
end
# ..............................................................................     
function _polynomial_gridfunction(n::Int, h::T; polynom=[0,1], deriv=0) where T <: Real
# ==============================================================================
#  polynomial_gridfunction(n, h; deriv) 
# ==============================================================================
    p = length(polynom) - 1
    deriv ≥ 0 || return T(0)
    deriv ≤ p || return T(0)

    f = deriv > 0 ? h^(deriv) * CamiMath.polynomial(polynom, n*h; deriv) :
                                CamiMath.polynomial(polynom, n*h)
    return f
    
end
# ..............................................................................
@doc raw"""
    gridfunction(ID::Int, n::Int, h::T; p=5, polynom=[0,1], deriv=0) where T <: Real

* `ID = 1`: exponential grid function,
```math
    g(t) = e^t - 1.0
```
* `ID = 2`: quasi-exponential grid function of degree `p` (linear grid for `p = 1`),
```math
    g(t) = t + \frac{1}{2}t^2 + ⋯ + \frac{1}{p!}t^p
```
* `ID = 3`: linear grid function,
```math
    g(t) = t
```
* `ID = 4`: polynomial grid function of degree `p = length(c)-1` defined by its [`CamiMath.polynom`](@reg) vector ``c = [c_0, c_1,c_2,⋯\ c_p]``,
```math
    g(t) = c_0 + c_1 t + c_2 t^2 + ⋯ + c_p t^p,
```
with ``c_0 ≡ 0`` because grid functions are defined to run through the origin, ``g(0) = 0``. 

The actual grid is given by 
```math
    x[n] = r_0 * g(t[n]),
```
where ``t[n] = (n-1) * h`` is the *ticks function* for the unit-based-array convention. 
NB. Note that ``x[1] = 0`` for all grid functions.
#### Examples:
```
h = 0.1
r = [gridfunction(1, n-1, h) for n=1:5]                            # exponential
 [0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032, 0.49182469764127035]

r = [gridfunction(2, n-1, h; p = 4) for n=1:5]  # quasi exponential (degree p=4)
 [0.0, 0.10517083333333321, 0.22140000000000004, 0.3498375, 0.49173333333333336]

r = [gridfunction(3, n-1, h) for n=1:5]              # linear
  [0.0, 0.1, 0.2, 0.3, 0.4]

r′= [gridfunction(3, n-1, h; deriv=1) for n=1:5]     # linear (first derivative)
   [0.1, 0.1, 0.1, 0.1, 0.1]

  r = [gridfunction(4, n-1, h; polynom = [0,1,1/2,1/6,1/24]) for n=1:5]  # polynomial of degree 4)
   [0.0, 0.10517083333333334, 0.2214, 0.3498375000000001, 0.49173333333333336]
```
"""
function gridfunction(ID::Int, n::Int, h::T; p=5, polynom=[0,1], deriv=0) where T <: Real

    return  ID == 1 ? _walterjohnson(n, h; deriv) :
            ID == 2 ? _jw_gridfunction(n, h; deriv, p) :
            ID == 3 ? _linear_gridfunction(n, h; deriv) :
            ID == 4 ? _polynomial_gridfunction(n, h; polynom, deriv) : throw(DomainError(ID, "unknown gridfunction"))

end

# ------------------------------------------------------------------------------
#            castGrid(ID, T, N; h=1, r0=0.01,  p=5, polynom=[0,1], epn=5, k=7)
# ------------------------------------------------------------------------------
    
function _gridspecs(ID::Int, N::Int, T::Type; h=1, r0=0.001, rmax=0, p=5, polynom=[0,1], epn=5, k=5, msg=true)
    
    rmax = ID == 1 ? r0 * _walterjohnson(N, h) :
           ID == 2 ? r0 * _jw_gridfunction(N, h; p) :
           ID == 3 ? r0 * _linear_gridfunction(N, h)  :
           ID == 4 ? r0 * _polynomial_gridfunction(N, h; polynom) : throw(DomainError(ID, "unknown gridfunction"))

    ID = ID ≠ 2 ? ID : p == 1 ? 3 : 2
    name = gridname(ID::Int)
    str_h = repr(h, context=:compact => true)
    str_r0 = repr(r0, context=:compact => true)
    str_rmax = repr(rmax, context=:compact => true)
    strA = "Grid created: $(name), $(T), rmax = "  * str_rmax * ", Ntot = $N, "
    
    return ID == 1 ? strA * "h = " * str_h * ", r0 = " * str_r0 :
           ID == 2 ? strA * "p = $p, h = " * str_h * ", r0 = " * str_r0 :
           ID == 3 ? strA * "p = 1, h = " * str_h * ", r0 = " * str_r0 :
           ID == 4 ? strA * "polynom = $(polynom), h = " * str_h * ", r0 = " * str_r0 : throw(DomainError(ID, "unknown gridfunction"))
    
end
# ..............................................................................
@doc raw"""
    castGrid(ID::Int, N::Int, T::Type; h=1, r0=1,  p=5, polynom=[0,1], epn=5, k=7, msg=false)
    castGrid(name::String, N::Int, T::Type; h=1, r0=1,  p=5, polynom=[0,1], epn=5, k=7, msg=false)

Method to create the Grid object

`ID = 1`: exponential grid,
`ID = 2`: quasi-exponential grid,
`ID = 3`: linear grid
`ID = 4`: polynomial grid
#### Examples:
```
julia> grid = castGrid(1, 4, Float64; h = 0.1, r0 = 1.0, msg=true);
Grid created: exponential, Float64, Rmax = 0.491825 a.u., Ntot = 4, h = 0.1, r0 = 1.0

julia> grid = castGrid(2, 4, Float64; p = 4, h = 0.1, r0 = 1.0, msg=true);
Grid created: quasi-exponential, Float64, Rmax = 0.491733 a.u., Ntot = 4, p = 4, h = 0.1, r0 = 1.0

julia> grid = castGrid(3, 4, Float64; h = 0.1, r0 = 1.0, msg=true);
Grid created: linear (uniform), Float64, Rmax = 0.4 a.u., Ntot = 4, p = 1, h = 0.1, r0 = 1.0

julia> grid.r′
4-element Vector{Float64}:
 0.1
 0.1
 0.1
 0.1

julia> grid = castGrid(4, 4, Float64; polynom=[0, 1, 1/2, 1/6, 1/24], h = 0.1, r0 = 1.0, msg=true);
Grid created: polynomial, Float64, Rmax = 0.491733 a.u., Ntot = 4, polynom = [0.0, 1.0, 0.5, 0.16666666666666666, 0.041666666666666664], h = 0.1, r0 = 1.0
 
```
"""
function castGrid(ID::Int, N::Int, T::Type; h=1, r0=0.001, p=5, polynom=[0,1], epn=5, k=5, msg=false)
# ==============================================================================
#  castGrid: creates the grid object
# ==============================================================================
    h = convert(T, h)
    r0 = convert(T, r0)
    polynom = convert.(T, polynom)
    epw = [convert.(T, trapezoidal_epw(n; rationalize=true)) for n=1:2:epn]
    name = gridname(ID)

    r  = r0 .* [gridfunction(ID, n-1, h; p, polynom) for n=1:N]
    r′ = r0 .* [gridfunction(ID, n-1, h; p, polynom, deriv=1) for n=1:N]     # r′= dr/dn
    r′′= r0 .* [gridfunction(ID, n-1, h; p, polynom, deriv=2) for n=1:N]     # r′′= d2r/dn2

    r[1] = T == BigFloat ? T(eps(Float64)) : T(eps(Float64))
    rmax = r[N]

    msg && println(_gridspecs(ID, N, T; h, r0, rmax, p, polynom, epn, k, msg))

    return Grid(ID, name, T, N, r, r′, r′′, h, r0, epn, epw, k)

end

# ........................ gridname(ID) ........................................
@doc raw"""
    gridname(ID::Int)
    
Name corresponding to the grid ID.
#### Example:
```
julia> gridname(2)
"quasi-exponential"
```
"""
function gridname(ID::Int)
# ==============================================================================
#  Name used for `Grid` of given `grid.ID`
# ==============================================================================
    
    return ID == 1 ? "exponential" :
           ID == 2 ? "quasi-exponential" :
           ID == 3 ? "linear (uniform)" :
           ID == 4 ? "polynomial" : throw(DomainError(ID, "unknown gridfunction"))
   
end

# =============== findIndex(rval, grid) ========================================

@doc raw"""
    findIndex(rval::T, grid::Grid{T}) where T<:Number

The grid index corresponding to the position `rval` on the `grid`.
#### Example:
```
julia> h = 0.1; r0 = 1.0;
julia> grid = castGrid(1, 4, Float64; h, r0);

julia> r = grid.r; println("r[3] = $(r[3])")
r[3] = 0.22140275816016985

julia> findIndex(0.222, grid)
3
```
"""
function findIndex(rval::T, grid::Grid{T}) where T<:Number
# kanweg
# ==============================================================================
#  grid index of rval, e.g., rval -> classical turning point
# ==============================================================================
    N = grid.N
    r = grid.r

    r[1] ≤ rval ≤ r[end] || throw(DomainError(rval, "rval outside grid range"))

    n = N
    while rval < r[n]     # below classical threshhold
        n > 1 ? n -= 1 : break
    end

    return n

end

# ------------------------------------------------------------------------------
#                       grid_differentiation(f, grid; k=3)
# ------------------------------------------------------------------------------

@doc raw"""
    grid_differentiation(f::Vector{T}, grid::Grid{T}; k=3) where T<:Real
    grid_differentiation(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int; k=3) where T<:Real
    grid_differentiation(f::Vector{T}, grid::Grid{T}, itr::UnitRange; k=3) where T<:Real

``k^{th}``-order lagrangian *differentiation* of the analytic function ``f``,
tabulated in forward order on a [`Grid`](@ref) of ``n`` points, ``f[1:n]``.
#### Example:
```
julia> ID = 3; # linear grid
julia> f = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0];
julia> grid = castGrid(ID, length(f), Float64; r0=1.0, h=1.0, k=3, msg=true);
Grid created: linear, Float64, Rmax = 6.0 a.u., Ntot = 6, p = 1, h = 1.0, r0 = 1.0

julia> f′= grid_differentiation(f, grid; k=3); println("f′= $(f′)")
f′= [0.0, 1.9999999999999991, 4.0, 6.000000000000001, 8.0, 10.0]
```
"""
function grid_differentiation(f::Vector{T}, grid::Grid{T}; k=3) where T<:Real

    r′= grid.r′

    f′= [fdiff_differentiation(f, T(i); k) for i ∈ eachindex(f)]

    return f′ ./ r′

end
function grid_differentiation(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int; k=3) where T<:Real

    f = f[n1:n2]
    r′= grid.r′[n1:n2]

    l = length(f)
    f′ = [fdiff_differentiation(f, T(v); k) for v=1:l]

    return f′ ./ r′

end
function grid_differentiation(f::Vector{T}, grid::Grid{T}, itr::UnitRange; k=3) where T<:Real

    return grid_differentiation(f, grid, itr.start, itr.stop; k)

end

# =============== grid_integration(f, grid, n1, n2) ===================

@doc raw"""
    grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
    grid_integration(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int) where T<:Real
    grid_integration(f::Vector{T}, grid::Grid{T}, itr::UnitRange) where T<:Real

Integral of the function ``f=[f_0,⋯\ f_n]`` tabulated on a [`Grid`](@ref)
using the trapezoidal rule optimized with endpoint correction by the
weightsvector `grid.epw`,
```math
    ∫_{0}^{r_n} f(r) dr = ∫_{0}^{n} f(x) r^{\prime}(x) dx,
```
where the latter integral corresponds to the optimized trapezoidal rule for a
uniform grid (see [`trapezoidal_integration`](@ref)). The rule is exact for
polynonials of degree ``d=0,\ 1,⋯\ k-1``, where ``k=`` `grid.epn`.
For ``k=1`` the rule reduces to the ordinary trapezoidal rule (weights = [1/2]).
#### Examples:
```
julia> ff(r) = sqrt(2.0/π) * exp(-r^2/2.0);

julia> grid1 = castGrid(1, 1000, Float64; h = 0.005, r0 = 0.1, msg=true);
Grid created: exponential, Float64, rmax = 14.7413, Ntot = 1000, h = 0.005, r0 = 0.1

julia> grid2 = castGrid(2, 1000, Float64; h = 0.005, r0 = 0.1, p=5, msg=true);
Grid created: quasi-exponential, Float64, rmax = 9.04167, Ntot = 1000, p = 5, h = 0.005, r0 = 0.1

julia> grid3 = castGrid(3, 1000, Float64; h = 0.1, r0 = 0.1, msg=true);
Grid created: linear (uniform), Float64, rmax = 10.0, Ntot = 1000, p = 1, h = 0.1, r0 = 0.1

julia> grid4 = castGrid(4, 1000, Float64; h = 0.1, r0 = 0.001, polynom=[0,0,1], msg=true);
Grid created: polynomial, Float64, rmax = 10.0, Ntot = 1000, polynom = [0.0, 0.0, 1.0], h = 0.1, r0 = 0.001

julia> r1 = grid1.r;
julia> r2 = grid2.r;
julia> r3 = grid3.r;
julia> r4 = grid4.r;
julia> f1 = [ff(r1[n]) for n=1:grid1.N];
julia> f2 = [ff(r2[n]) for n=1:grid2.N];
julia> f3 = [ff(r3[n]) for n=1:grid3.N];
julia> f4 = [ff(r4[n]) for n=1:grid4.N];
julia> o1 = grid_integration(f1, grid1);
julia> o2 = grid_integration(f2, grid2);
julia> o3 = grid_integration(f3, grid3, 1:900);
julia> o4 = grid_integration(f4, grid4, 1:900);

julia> println("integral on " * grid1.name * " grid = ", o1)
integral on exponential grid = 1.0

julia> println("integral on " * grid2.name * " grid = ", o2)
integral on quasi-exponential grid: 1.0

julia> println("integral on " * grid3.name * " grid = ", o3)
integral on linear (uniform) grid = 1.000000000000003

julia> println("integral on " * grid3.name * " grid = ", o4)
integral on polynomial grid = 1.0000000000000013
```
"""
function grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
# ==============================================================================
#  trapezoidal integral over the grid indices [n1:n2] with 1 ≤ n1,n2 ≤ N
# ==============================================================================

    r′= grid.r′
    N = grid.N

    epn = grid.epn   # endpoint number
    epw = grid.epw   # endpoint weights array
    
    if N ≥ 2epn
        epi = epn÷2+1        # index endpoint weights
    else
        epn = Base.isodd(N÷2) ? (N÷2) : N÷2-1           # endpoint number
        epi = Base.isodd(N÷2) ? (N÷2)÷2+1 : (N÷2-1)÷2+1 # index endpoint weights
    end

    w = Base.ones(T,N)
    w[1:epn] = epw[epi]
    w[N-epn+1:N] = Base.reverse(epw[epi])

    return LinearAlgebra.dot(f .* r′, w)

end
function grid_integration(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int) where T<:Real
# ==============================================================================
#  trapezoidal integral over the grid indices [n1:n2] with 1 ≤ n1,n2 ≤ N
# ==============================================================================
    f = f[n1:n2]
    r′= grid.r′[n1:n2]
    n = n2 - n1 + 1

    n > 1 || return 0

    epn = grid.epn   # endpoint number
    epw = grid.epw   # endpoint weights array

    if n ≥ 2epn
        epi = epn÷2+1        # index endpoint weights
    else
        epn = Base.isodd(n÷2) ? (n÷2) : n÷2-1           # endpoint number
        epi = Base.isodd(n÷2) ? (n÷2)÷2+1 : (n÷2-1)÷2+1 # index endpoint weights
    end

    w = Base.ones(T,n)
    w[1:epn] = epw[epi]
    w[end-epn+1:end] = Base.reverse(epw[epi])

    return LinearAlgebra.dot(f .* r′, w)

end
function grid_integration(f::Vector{T}, grid::Grid{T}, itr::UnitRange) where T<:Real

    return grid_integration(f, grid, itr.start, itr.stop)

end
