# SPDX-License-Identifier: MIT

# Copyright (c) 2024 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ==============================================================================
#                               grid.jl
# ============================================================================== 

# ------------------------------------------------------------------------------
#    Grid (ID, name, Type, N, r, r′, r′′, r′′, h, r0, epn, epw, k, p, polynom)
# ------------------------------------------------------------------------------

"""
    Grid(ID, name, T, N, r, r′, r′′, h, r0, epn, epw, k, p, polynom)

Type with fields:
* `.ID`:    grid identifer name (`::Int`)
* `.name`:  grid identifer name (`::String`)
* `.T`:     gridtypename (`::Type`)
* `.N`:     number of grid points (`::Int`)
* `.r  `:   tabulated grid function (`::Vector{T}`)
* `.r′ `:   tabulated first derivative of grid function (`::Vector{T}`)
* `.r′′`:   tabulated second derivative of grid function (`::Vector{T}`)
* `.h` :    grid step multiplyer (`::T`)
* `.r0`:    grid scale factor (`::T`)
* `.epn`:   number of endpoints used for trapezoidal endpoint correction (must be odd) (`::Int`)
* `.epw`:   trapezoidal endpoint weights for n=1:epn (`::Vector{Vector{T}}`)
* `.k`:     finite-difference order (`::Int`)
* `.p`:     only for quasi-exponential grid; truncation power (`::Int`)
* `.polynom`: only for polynomial grid: polynom (`::Vector{T}`)
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
p::Int
polynom::Vector{T}

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

`CamiDiff` offers three internal grid functions:

* `ID = 1`: exponential grid function,
```math
    g(t) = e^t - 1
```
* `ID = 2`: quasi-exponential grid function of degree `p` (linear grid for `p = 1`),
```math
    g(t) = t + \frac{1}{2}t^2 + ⋯ + \frac{1}{p!}t^p
```
* `ID = 3`: linear grid function,
```math
    g(t) = t
```
* `ID = 4`: polynomial grid function of degree `p = length(c)-1` defined by its [`CamiMath.polynom`](@extref) vector ``c = [c_0, c_1,c_2,⋯\ c_p]``,
```math
    g(t) = c_0 + c_1 t + c_2 t^2 + ⋯ + c_p t^p,
```
with ``c_0 ≡ 0`` because, *by definition*, all [`gridfunction`](@ref)`s` run through the origin, ``g(0) = 0``. 

The actual [`Grid`](@ref) is given by 
```math
    r[n] = r_0 * g(t[n]),
```
where ``t[n] = (n-1) * h`` is the *ticks function* for the unit-based indexing of [Julia](http://julialang.org).

NB. All [`gridfunction`](@ref)`s` satisfy the properties ``t[1] = 0`` and ``r[1] = 0``.
#### Examples:
```
julia> h = 0.1; r0=1.0; N=4;

julia> r = r0 .* [gridfunction(1, n-1, h) for n=1:N]
[0.0, 0.10517091807564771, 0.22140275816016985, 0.3498588075760032]

julia> r′ = r0 .* [gridfunction(1, n-1, h; deriv=1) for n=1:N]
[0.1, 0.11051709180756478, 0.122140275816017, 0.13498588075760032]

julia> r′′= r0 .* [gridfunction(1, n-1, h; deriv=2) for n=1:N]
[0.010000000000000002, 0.011051709180756479, 0.012214027581601701, 0.013498588075760034]

julia> r = r0 .* [gridfunction(4, n-1, h; polynom=[0,0,1]) for n=1:N]
[0.0, 0.010000000000000002, 0.04000000000000001, 0.09000000000000002]

julia> r′ = r0 .* [gridfunction(4, n-1, h; polynom=[0,0,1], deriv=1) for n=1:N]
[0.0, 0.020000000000000004, 0.04000000000000001, 0.06000000000000001]

julia> r′′= r0 .* [gridfunction(4, n-1, h; polynom=[0,0,1], deriv=2) for n=1:N]
[0.020000000000000004, 0.020000000000000004, 0.020000000000000004, 0.020000000000000004]
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
    
function _gridspecs(ID::Int, N::Int, T::Type; h=1, r0=1, rmax=0, p=5, polynom=[0,1], epn=5, k=5)
    
    rmax = ID == 1 ? r0 * _walterjohnson(N, h) :
           ID == 2 ? r0 * _jw_gridfunction(N, h; p) :
           ID == 3 ? r0 * _linear_gridfunction(N, h)  :
           ID == 4 ? r0 * _polynomial_gridfunction(N, h; polynom) : 
                           throw(DomainError(ID, "unknown gridfunction"))

    ID = ID ≠ 2 ? ID : p == 1 ? 3 : 2
    str_d = ID == 4 ? "of degree $(length(polynom)-1)" : nothing
    name = ID == 4 ? gridtypename(ID::Int) * " of degree $(length(polynom)-1)" : gridtypename(ID::Int)
    str_h = repr(h, context=:compact => true)
    str_r0 = repr(r0, context=:compact => true)
    str_rmax = repr(rmax, context=:compact => true)

    strA = "Grid: $(name), $(T), rmax = "  * str_rmax * ", Ntot = $N, "
    
    return ID == 1 ? strA * "h = " * str_h * ", r0 = " * str_r0 :
           ID == 2 ? strA * "p = $p, h = " * str_h * ", r0 = " * str_r0 :
           ID == 3 ? strA * "p = 1, h = " * str_h * ", r0 = " * str_r0 :
           ID == 4 ? strA * "polynom = $(polynom), h = " * str_h * ", r0 = " * str_r0 : 
                             throw(DomainError(ID, "unknown gridfunction"))
    
end
# ..............................................................................
@doc raw"""
    castGrid(ID::Int, N::Int, T::Type; h=1, r0=1,  p=5, polynom=[0,1], epn=5, k=7, msg=false)
    castGrid(name::String, N::Int, T::Type; h=1, r0=1, p=5, polynom=[0,1], epn=5, k=7, msg=false)

Method to create the [`Grid`](@ref) object

`ID = 1`: exponential,
`ID = 2`: quasi-exponential,
`ID = 3`: linear (uniform)
`ID = 4`: polynomial
#### Examples:
```
julia> grid = castGrid(1, 1000, Float64; h = 0.005, r0 = 0.1, msg=true);
Grid: exponential, Float64, rmax = 14.7413, Ntot = 1000, h = 0.005, r0 = 0.1

julia> grid = castGrid("exponential", 1000, Float64; h = 0.005, r0 = 0.1, msg=true);
Grid: exponential, Float64, rmax = 14.7413, Ntot = 1000, h = 0.005, r0 = 0.1

julia> grid = castGrid(2, 1000, Float64; h = 0.005, r0 = 0.1, p=5, msg=true);
Grid: quasi-exponential, Float64, rmax = 9.04167, Ntot = 1000, p = 5, h = 0.005, r0 = 0.1

julia> grid = castGrid(3, 1000, Float64; h = 0.1, r0 = 0.1, msg=true);
Grid: linear (uniform), Float64, rmax = 10.0, Ntot = 1000, p = 1, h = 0.1, r0 = 0.1

julia> grid = castGrid(4, 1000, Float64; h = 0.1, r0 = 0.001, polynom=[0,0,1], msg=true);
Grid: polynomial of degree 2, Float64, rmax = 10.0, Ntot = 1000, polynom = [0.0, 0.0, 1.0], h = 0.1, r0 = 0.001

julia> grid.r[1:4]
4-element Vector{Float64}:
 [2.220446049250313e-16, 1.0000000000000003e-5, 4.000000000000001e-5, 9.000000000000003e-5]

julia> grid.r′[1:4]
4-element Vector{Float64}:
 [0.0, 2.0000000000000005e-5, 4.000000000000001e-5, 6.0000000000000015e-5]
 
julia> grid.r′′[1:4]
4-element Vector{Float64}:
 [2.0000000000000005e-5, 2.0000000000000005e-5, 2.0000000000000005e-5, 2.0000000000000005e-5] 
```
"""
function castGrid(ID::Int, N::Int, T::Type; h=1, r0=1, p=5, polynom=[0,1], epn=5, k=5, msg=false)
# ==============================================================================
#  castGrid: creates the grid object
# ==============================================================================
    h = convert(T, rationalize(h))
    r0 = convert(T, rationalize(r0))
    polynom = convert(Vector{T}, rationalize.(polynom))
    epw = [convert.(T, trapezoidal_epw(n; rationalize=true)) for n=1:2:epn]
    name = gridtypename(ID)

    r  = r0 .* T[gridfunction(ID, n-1, h; p, polynom) for n=1:N]
    r′ = r0 .* T[gridfunction(ID, n-1, h; p, polynom, deriv=1) for n=1:N]     #  r′= dr/dn
    r′′= r0 .* T[gridfunction(ID, n-1, h; p, polynom, deriv=2) for n=1:N]     # r′′= d²r/dn²

    # r[1] = T == BigFloat ? T(eps(Float64)) : T(eps(Float64)) # quasi zero - kanweg
    rmax = r[N]

    msg && println(_gridspecs(ID, N, T; h, r0, rmax, p, polynom, epn, k))

    return Grid(ID, name, T, N, r, r′, r′′, h, r0, epn, epw, k, p, polynom)

end
function castGrid(name::String, N::Int, T::Type; h=1, r0=1, p=5, polynom=[0,1], epn=5, k=5, msg=false)
    
    ID = name == "exponential" ? 1 :
         name == "quasi-exponential" ? 2 :
         name == "linear" ? 3 :
         name == "polynomial" ? 4 : throw(DomainError(ID, "unknown gridfunction"))

    return castGrid(ID, N, T; h, r0, p, polynom, epn, k, msg)
   
end

# ........................ gridtypename(ID) ........................................
@doc raw"""
    gridtypename(ID::Int)
    
Name corresponding to the [`Grid`](@ref) ID.
#### Example:
```
julia> gridtypename(2)
"quasi-exponential"
```
"""
function gridtypename(ID::Int)
# ==============================================================================
#  Name used for `Grid` of given `grid.ID`
# ==============================================================================
    
    return ID == 1 ? "exponential" :
           ID == 2 ? "quasi-exponential" :
           ID == 3 ? "linear (uniform)" :
           ID == 4 ? "polynomial" : throw(DomainError(ID, "unknown gridfunction"))
   
end

# ........................ gridtypeID(name) ........................................
@doc raw"""
    gridtypeID(name::String)
    
ID corresponding to the [`gridtypename`](@ref).
#### Example:
```
julia> gridtypeID("quasi-exponential")
2
```
"""
function gridtypeID(name::String)
# ==============================================================================
#  ID corresponding to the gridtypename 'name'
# ==============================================================================
    
    return name == "exponential" ? 1 :
           name == "quasi-exponential" ? 2 :
           name == "linear (uniform)" ? 3 :
           name == "polynomial" ? 4 : throw(DomainError(name, "unknown gridfunction"))
   
end

# ------------------------------------------------------------------------------
#                       gridPos(rval, grid)
# ------------------------------------------------------------------------------

@doc raw"""
    gridPos(rval::T, grid::Grid{T}) where T<:Number

The approximate grid position defined as the largest integer `n` satisfying the 
condition `grid.r[n] < rval` on the [`Grid`](@ref).
#### Example:
```
julia> h = 0.1; r0 = 1.0;

julia> grid = castGrid(1, 4, Float64; h, r0);

julia> r = grid.r; println("r[3] = $(r[3])")
r[3] = 0.22140275816016985

julia> gridPos(0.222, grid)
3
```
"""
function gridPos(rval::T, grid::Grid{T}) where T<:Real
    
    N = grid.N
    r = grid.r

    r[1] ≤ rval ≤ r[end] || throw(DomainError(rval, "rval = $(rval) outside range $(r[1]) ≤ rval ≤ $(r[end])"))
    
    n = 1
    m = N
    Δn = (m-n)÷2
    
    while Δn ≥ 1
        if r[n+Δn] ≤ rval
            n += Δn 
        elseif r[m-Δn] ≥ rval
            m -= Δn
        else 
            error("Error: search for intersection point undecided")
        end
        Δn = (m-n)÷2
    end

    return n

end

# ------------------------------------------------------------------------------
#                       fracPos(n, rval, grid; ϵ = 1e-8, k = 7)
# ------------------------------------------------------------------------------

@doc raw"""
    fracPos(n::Int, rval::T, grid::Grid{T}; ϵ = 1e-8, k = 7) where T<:Real

Fractional grid offset with respect to [`Grid`](@ref) position `n`.
"""
function fracPos(n::Int, rval::T, grid::Grid{T}; ϵ = T(1e-8), k = 7) where T<:Real
    
    nul = T(0)
    one = T(1)
    two = T(2)

    ID = grid.ID
    r0 = grid.r0
     h = grid.h
     p = grid.p
    polynom = grid.polynom 

    polynom = [r0*gridfunction(ID, n-1, h; p, polynom, deriv=d)/factorial(d) for d=0:k]

    imin = nul   
    imax = one   
    r1 = CamiMath.polynomial(polynom, imin)
    r2 = CamiMath.polynomial(polynom, imax) 
    r1 ≤ rval ≤ r2 || throw(DomainError(rval, "rval = $(rval) outside range $(r1) ≤ rval ≤ $(r2)"))
    
    Δn = (imax+imin)/two
    while imax-imin > ϵ
        if CamiMath.polynomial(polynom, Δn) ≤ rval
            imin = Δn
        else
            imax = Δn
        end
        Δn = (imax+imin)/two
    end

    return Δn
    
end

# ------------------------------------------------------------------------------
#                       grid_interpolation(f, rval, grid; k=5)
# ------------------------------------------------------------------------------

@doc raw"""
    grid_interpolation(f::Vector{T}, rval::T, grid::Grid{T}; k=5) where T<:Real

Interpolated value for f(rval), with `rval` not on the [`Grid`](@ref) (`rval ∉ grid.r`).
"""
function grid_interpolation(f::Vector{T}, rval::T, grid::Grid{T}; k=5) where T<:Real

    N = grid.N
    k = min(k,N÷2)
    k > 0 || throw(DomainError("k = $k violates k ≥ 1 as required for lagrangian interpolation"))

    n = gridPos(rval, grid)
    σ = fracPos(n, rval, grid; ϵ = 1e-12, k = 7)
    α = fdiff_interpolation_expansion_polynom(σ, k, fwd)
    Fk = convert(Vector{T}, fdiff_expansion_weights(α, fwd, reg))
    o = LinearAlgebra.dot(Fk, f[n:n+k])
    
    return o

end

# ------------------------------------------------------------------------------
#                       grid_differentiation(f, grid; k=5)
# ------------------------------------------------------------------------------

function _regularize_origin(f′::Vector{T}, r′::Vector{T}, k::Int) where T<:Real

    o = f′ ./ r′

    if isinf(o[1])
        α = fdiff_interpolation_expansion_polynom(1, k-1, fwd) 
        Fk = convert(Vector{T}, fdiff_expansion_weights(α, fwd, reg))
        o[1] = LinearAlgebra.dot(Fk, o[2:k+1]) 
    end
    
    return o

end
@doc raw"""
    grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:real

``k^{th}``-order lagrangian *derivative* of the function ``f(r)``,
* `f[1:N]` : the function `f(r)` tabulated in forward order on a [`Grid`](@ref) of ``N`` points.
#### Example:
```
julia> grid = castGrid(3, 1001, Float64; h=2π/1000.0, r0=1.0, msg=false);

julia> f = [sin(grid.r[i]) for i=1:grid.N]

julia> f′ = [cos(grid.r[i]) for i=1:grid.N]

julia> grid_differentiation(f, grid) ≈ f′
true
```
    grid_differentiation(f::Vector{T}, grid::Grid{T}, r::T, notation=fwd; k=5) where T<:Real

``k^{th}``-order lagrangian *derivative* of the function ``f(r)`` at position `r`. 
* `f[1:N]` : the function `f(r)` tabulated in forward order on a [`Grid`](@ref) of ``N`` points.
* `fwd` using fwd-difference notation  
* `bwd` using bwd-difference notation

#### Example:
```
julia> grid = castGrid(3, 1001, Float64; h=2π/1000.0, r0=1.0, msg=false);

julia> f = [sin(grid.r[i]) for i=1:grid.N];

julia> grid_differentiation(f, grid, float(π), fwd) ≈ cos(π)
true
```
    grid_differentiation(f::Vector{T}, grid::Grid{T}, n::Int, notation=fwd; k=5) where T<:Real

``k^{th}``-order lagrangian *derivative* of the function ``f(r)`` at grid position `n`. 
* `f[1:N]` : the function `f(r)` tabulated in forward order on a [`Grid`](@ref) of ``N`` points.
* `fwd` using fwd-difference notation  
* `bwd` using bwd-difference notation

#### Example:
```
julia> grid = castGrid(3, 1001, Float64; h=2π/1000.0, r0=1.0, msg=false);

julia> f = [sin(grid.r[i]) for i=1:grid.N];

julia> f′ = [cos(grid.r[i]) for i=1:grid.N]

julia> grid_differentiation(f, grid, 500, fwd) ≈ f′[500]
true
```
    grid_differentiation(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int; k=5) where T<:Real
    grid_differentiation(f::Vector{T}, grid::Grid{T}, itr::UnitRange; k=5) where T<:Real

``k^{th}``-order lagrangian *derivative* of the function ``f(r)`` over the grid range `n1 ≤ n ≤ n2`.
* `f[1:N]` : the function `f(r)` tabulated in forward order on a [`Grid`](@ref) of ``N`` points.
* `n1=itr.start`, `n2=itr.stop`.  


"""
function grid_differentiation(f::Vector{T}, grid::Grid{T}; k=5) where T<:Real

    N = grid.N
    r′= grid.r′
    k = min(k,N÷2)
    k > 0 || throw(DomainError("k = $k violates k ≥ 1 as required for lagrangian interpolation"))
    
    α = fdiff_differentiation_expansion_polynom(0, k, fwd)
    β = fdiff_differentiation_expansion_polynom(0, k, bwd)
    Fk = convert(Vector{T}, fdiff_expansion_weights(α, fwd, reg))
    Bkrev = convert(Vector{T}, fdiff_expansion_weights(β, bwd, rev))

    f′= [LinearAlgebra.dot(Fk, f[n:n+k]) for n=1:N-k]
    g′= [LinearAlgebra.dot(Bkrev, f[n-k:n]) for n=N-k+1:N]
    f′= append!(f′, g′)
    
    return _regularize_origin(f′, r′, k)

end
function grid_differentiation(f::Vector{T}, grid::Grid{T}, n::Int, notation=fwd; k=5) where T<:Real

    N = grid.N
    r′= grid.r′
    k = min(k,N÷2)
    k > 0 || throw(DomainError("k = $k violates k ≥ 1 as required for lagrangian interpolation"))
    
    if CamiMath.isforward(notation)
        α = fdiff_differentiation_expansion_polynom(0, k, fwd)
        Fk = convert(Vector{T}, fdiff_expansion_weights(α, fwd, reg))
        f′= LinearAlgebra.dot(Fk, f[n:n+k])/r′[n]
    else
        β = fdiff_differentiation_expansion_polynom(0, k, bwd)
        Bkrev = convert(Vector{T}, fdiff_expansion_weights(β, bwd, rev))
        f′= LinearAlgebra.dot(Bkrev, f[n-k:n])/r′[n]
    end
 
    return f′

end
function grid_differentiation(f::Vector{T}, grid::Grid{T}, rval::T, notation=fwd; k=5) where T<:Real

    n = gridPos(rval, grid)
    Δn = fracPos(n, rval, grid; ϵ = T(1e-8), k)
    v = n + Δn

    N = grid.N
    r = grid.r
    r′= grid.r′
    k = min(k,N÷2)
    k > 0 || throw(DomainError("k = $k violates k ≥ 1 as required for lagrangian interpolation"))
    
    if CamiMath.isforward(notation)
        σ = n-v
        α = fdiff_differentiation_expansion_polynom(σ, k, fwd)
        Fk = convert(Vector{T}, fdiff_expansion_weights(α, fwd, reg))
        r′= LinearAlgebra.dot(Fk, r[n:n+k])
        f′= LinearAlgebra.dot(Fk, f[n:n+k])/r′
    else
        σ = v-n
        β = fdiff_differentiation_expansion_polynom(σ, k, bwd)
        Bkrev = convert(Vector{T}, fdiff_expansion_weights(β, bwd, rev))
        r′= LinearAlgebra.dot(Bkrev, r[n-k:n])
        f′= LinearAlgebra.dot(Bkrev, f[n-k:n])/r′
    end
 
    return f′

end
function grid_differentiation(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int; k=5) where T<:Real

    N = grid.N
    r′= [grid.r′[n] for n=n1:n2]

    k = min(k,N÷2)
    1 ≤ n1 ≤ n2 ≤ N || throw(DomainError(n1,n2))
    k > 0 || throw(DomainError("k = $k violates k ≥ 1 as required for lagrangian interpolation"))
    
    α = fdiff_differentiation_expansion_polynom(0, k, fwd)
    β = fdiff_differentiation_expansion_polynom(0, k, bwd)
    Fk = convert(Vector{T}, fdiff_expansion_weights(α, fwd, reg))
    Bkrev = convert(Vector{T}, fdiff_expansion_weights(β, bwd, rev))

    N-k ≥ n2 || println("Warning: N-k < n2, (n2+k-N = ", n2+k-N, ")")

    n3 = min(n2,N-k)
    
 #   f′= [LinearAlgebra.dot(Fk, f[n:n+k]) for n=n1:n3]
 #   g′= [LinearAlgebra.dot(Bkrev, f[n-k:n]) for n=n2-k+1:n2]
 #   n0=1
 #   if n0+k ≤ n1
 #       println("a")
 #       f′ = [sum(Bkrev .* f[n-k+1:n]) for n=n1:n2]
 #   else
 #       println("a")
 #       f′ = [sum(Fk .* f[n:n+k]) for n=n1:n1+k]
 #       g′ = [sum(Bkrev .* f[n-k+1:n]) for n=n1+k+1:n2]
 #       append!(f′, g′)
 #   end    

    f′= [sum(Fk .* f[n:n+k]) for n=n1:n3]
    g′= [sum(Bkrev .* f[n-k:n]) for n=n2-k+1:n2]
#    println("[f′[end] = ", f′[])
    f′= append!(f′, g′)[1:n2-n1+1]
    
    return _regularize_origin(f′, r′, k)

end
function grid_differentiation(f::Vector{T}, grid::Grid{T}, itr::UnitRange; k=5) where T<:Real

    return grid_differentiation(f, grid, itr.start, itr.stop; k)

end

# =============== grid_integration(f, grid, n1, n2) ===================

@doc raw"""
    grid_integration(f::Vector{T}, grid::Grid{T}) where T<:Real
    grid_integration(f::Vector{T}, grid::Grid{T}, n1::Int, n2::Int) where T<:Real
    grid_integration(f::Vector{T}, grid::Grid{T}, itr::UnitRange) where T<:Real

Integral of the analytic function ``f(r)`` tabulated on a generally nonlinear 
[`Grid`](@ref) and evaluated with the generalized trapezoidal rule optimized 
with endpoint correction by the weightsvector `grid.epw`,
```math
    ∫_{0}^{r_n} f(r) dr = ∫_{0}^{n} f(x) r^{\prime}(x) dx.
```
Here, the latter integral corresponds to the optimized trapezoidal rule for a
uniform grid (see [`trapezoidal_integration`](@ref)). The rule is exact for
polynomial functions of degree ``d=0,\ 1,⋯\ k-1``, where ``k=`` `grid.epn`.
For ``k=1`` the rule reduces to the ordinary trapezoidal rule (weights = [1/2]).
#### Examples:
```
julia> ftest(r) = sqrt(2.0/π) * exp(-r^2/2.0);

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
julia> f1 = [ftest(r1[n]) for n=1:grid1.N];
julia> f2 = [ftest(r2[n]) for n=1:grid2.N];
julia> f3 = [ftest(r3[n]) for n=1:grid3.N];
julia> f4 = [ftest(r4[n]) for n=1:grid4.N];
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