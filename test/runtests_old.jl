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

using CamiDiff
using Test # to use runtests with @testset

import CamiMath
fwd = CamiMath.fwd
bwd = CamiMath.bwd
reg = CamiMath.reg
rev = CamiMath.rev

@testset "CamiDiff.jl" begin 

    grid1 = castGrid(1, 4, Float64; h = 0.1, r0 = 2.0);
    grid2 = castGrid(2, 4, Float64; h = 0.1, r0 = 2.0, p=1);
    grid3 = castGrid(3, 4, Float64; h = 0.1, r0 = 2.0);
    grid4 = castGrid("polynomial", 4, Float64; h = 0.1, r0 = 2.0, polynom=[0, 1], msg=true);
    @test [grid2.r, grid2.r′, grid2.r′′] ≈ [grid3.r, grid3.r′, grid3.r′′] ≈ [grid4.r, grid4.r′, grid4.r′′]
    @test [grid1.r, grid1.r′, grid1.r′′] == [[0.0, 0.21034183615129542, 0.4428055163203397, 0.6997176151520064], [0.2, 0.22103418361512955, 0.244280551632034, 0.26997176151520064], [0.020000000000000004, 0.022103418361512958, 0.024428055163203403, 0.02699717615152007]]
    @test grid1.name == "exponential"
    @test gridPos(0.25, grid1) == 2
    
    @test gridtypename(1) == "exponential"
    @test gridtypeID("exponential") == 1
    @test_throws DomainError gridPos(100.0, grid1) == 220
    @test_throws DomainError castGrid(5, 1000, Float64)
    @test_throws DomainError gridtypename(5) 
    @test_throws DomainError gridtypeID("logarithm") 
#   -----------------------------------------------------------------------------------------
    gaussian(r) = sqrt(2.0/π) * exp(-r^2/2.0);
    grid1 = castGrid(1, 1000, Float64; h = 0.005, r0 = 0.1, msg=false);
    grid2 = castGrid(2, 1000, Float64; h = 0.005, r0 = 0.1, p=5, msg=false);
    grid3 = castGrid(3, 1000, Float64; h = 0.1, r0 = 0.1, msg=false);
    grid4 = castGrid(4, 1000, Float64; h = 0.1, r0 = 0.001, polynom=[0,0,1], msg=false);
    r1 = grid1.r;
    r2 = grid2.r;
    r3 = grid3.r;
    r4 = grid4.r;
    f1 = [gaussian(r1[n]) for n=1:grid1.N];
    f2 = [gaussian(r2[n]) for n=1:grid2.N];
    f3 = [gaussian(r3[n]) for n=1:grid3.N];
    f4 = [gaussian(r4[n]) for n=1:grid4.N];
#   -----------------------------------------------------------------------------------------
    o1 = grid_integration(f1, grid1);
    o2 = grid_integration(f2, grid2);
    o3 = grid_integration(f3, grid3);
    o4 = grid_integration(f4, grid4);
    @test o1 ≈ 1.0
    @test o2 ≈ 1.0
    @test o3 ≈ 1.0
    @test o4 ≈ 1.0
    o1 = grid_integration(f1, grid1, 1:900);
    o2 = grid_integration(f2, grid2, 1:900);
    o3 = grid_integration(f3, grid3, 1:900);
    o4 = grid_integration(f4, grid4, 1:900);
    @test o1 ≈ 1.0
    @test o2 ≈ 1.0
    @test o3 ≈ 1.0
    @test o4 ≈ 1.0
#   -----------------------------------------------------------------------------------------
    f′1 = -r1 .* f1;
    f′2 = -r2 .* f2;
    f′3 = -r3 .* f3;
    f′4 = -r4 .* f4;
    o1 = grid_differentiation(f1, grid1);
    @test f′1 ≈ o1
    a1 = f′1 ./ o1
    println("a1: ", a1[992:1000])
    println("f′1: ", f′1[992:1000])
    println("o1: ", o1[992:1000])
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);
    @test f′2 ≈ o2
    @test f′3 ≈ o3
    @test f′4 ≈ o4 
 #   a2 = f′2 ./ o2
 #   a3 = f′3 ./ o3
 #   a4 = f′4 ./ o4
#    println("a2: ", a2[992:1000])
#    println("a3: ", a3[992:1000])
#    println("a4: ", a4[992:1000])
    o1 = grid_differentiation(f1, grid1, 1:900);
    o2 = grid_differentiation(f2, grid2, 1:900);
    o3 = grid_differentiation(f3, grid3, 1:900);
    o4 = grid_differentiation(f4, grid4, 1:900);
    @test f′1[1:900] ≈ o1
    @test f′2[1:900] ≈ o2
    @test f′3[1:900] ≈ o3
    @test f′4[1:900] ≈ o4    
#   -----------------------------------------------------------------------------------------
    exponential(r) = exp(-r);
    grid1 = castGrid(1, 1000, Float64; h = 0.01, r0 = 0.001, msg=true);
    grid2 = castGrid(2, 1000, Float64; h = 0.01, r0 = 0.02, p=6, msg=true);
    grid3 = castGrid(3, 1000, Float64; h = 0.01, r0 = 2.0, msg=true);
    grid4 = castGrid(4, 1000, Float64; h = 0.01, r0 = 0.2, polynom=[0,0,1], msg=true);
    r1 = grid1.r;
    r2 = grid2.r;
    r3 = grid3.r;
    r4 = grid4.r;
    f1 = [exponential(r1[n]) for n=1:grid1.N];
    f2 = [exponential(r2[n]) for n=1:grid2.N];
    f3 = [exponential(r3[n]) for n=1:grid3.N];
    f4 = [exponential(r4[n]) for n=1:grid4.N];
#   -----------------------------------------------------------------------------------------
    rval=2.1
    o1 = grid_interpolation(f1, rval, grid1; k=5)
    o2 = grid_interpolation(f2, rval, grid2; k=5)
    o3 = grid_interpolation(f3, rval, grid3; k=5)
    o4 = grid_interpolation(f4, rval, grid4; k=5)
    @test o1 ≈ exponential(rval)
    @test o2 ≈ exponential(rval)
    @test o3 ≈ exponential(rval)
    @test o4 ≈ exponential(rval)
    o1 = grid_integration(f1, grid1);
    o2 = grid_integration(f2, grid2);
    o3 = grid_integration(f3, grid3);
    o4 = grid_integration(f4, grid4);
    @test o1 ≈ 1.0
    @test o2 ≈ 1.0
    @test o3 ≈ 1.0
    @test o4 ≈ 1.0
    o1 = grid_integration(f1, grid1, 1:990);
    o2 = grid_integration(f2, grid2, 1:950);
    o3 = grid_integration(f3, grid3, 1:950);
    o4 = grid_integration(f4, grid4, 1:990);
    @test o1 ≈ 1.0
    @test o2 ≈ 1.0
    @test o3 ≈ 1.0
    @test o4 ≈ 1.0
#   -----------------------------------------------------------------------------------------
    f′1 = -f1;
    f′2 = -f2;
    f′3 = -f3;
    f′4 = -f4;
    o1 = grid_differentiation(f1, grid1);
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);
    @test f′1 ≈ o1
    @test f′2 ≈ o2
    @test f′3 ≈ o3
    @test f′4 ≈ o4
#    a1 = f′1 ./ o1
#    a2 = f′2 ./ o2
#    a3 = f′3 ./ o3
#    a4 = f′4 ./ o4
#    println("a1: ", a1[992:1000])
#    println("a2: ", a2[992:1000])
#    println("a3: ", a3[992:1000])
#    println("a4: ", a4[992:1000])
    o1 = grid_differentiation(f1, grid1, 10:900); 
    o2 = grid_differentiation(f2, grid2, 10:900);
    o3 = grid_differentiation(f3, grid3, 10:900);
    o4 = grid_differentiation(f4, grid4, 10:900);
    @test f′1[10:900] ≈ o1 
    @test f′2[10:900] ≈ o2
    @test f′3[10:900] ≈ o3
    @test f′4[10:900] ≈ o4
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, 5, Float64; h = 0.01, r0 = 0.001, msg=true);
    grid2 = castGrid(2, 5, Float64; h = 0.01, r0 = 0.02, p=5, msg=true);
    grid3 = castGrid(3, 5, Float64; h = 0.01, r0 = 2.0, msg=true);
    grid4 = castGrid(4, 5, Float64; h = 0.01, r0 = 0.2, polynom=[0,0,1], msg=true);
    r1 = grid1.r;
    r2 = grid2.r;
    r3 = grid3.r;
    r4 = grid4.r;
    f1 = [exponential(r1[n]) for n=1:grid1.N];
    f2 = [exponential(r2[n]) for n=1:grid2.N];
    f3 = [exponential(r3[n]) for n=1:grid3.N];
    f4 = [exponential(r4[n]) for n=1:grid4.N];
#   -----------------------------------------------------------------------------------------
    o1 = grid_integration(f1, grid1);
    o2 = grid_integration(f2, grid2);
    o3 = grid_integration(f3, grid3);
    o4 = grid_integration(f4, grid4);
    @test o1 ≈ 4.081028048564596e-5
    @test o2 ≈ 0.000815888856047799
    @test o3 ≈ 0.0768862163847329
    @test o4 ≈ 0.0003199456063034396
    o1 = grid_integration(f1, grid1, 1:5);
    o2 = grid_integration(f2, grid2, 1:5);
    o3 = grid_integration(f3, grid3, 1:5);
    o4 = grid_integration(f4, grid4, 1:5);
    @test o1 ≈ 4.081028048564596e-5
    @test o2 ≈ 0.000815888856047799
    @test o3 ≈ 0.0768862163847329
    @test o4 ≈ 0.0003199456063034396
    #   -----------------------------------------------------------------------------------------
    f′1 = -f1;
    f′2 = -f2;
    f′3 = -f3;
    f′4 = -f4;
    o1 = grid_differentiation(f1, grid1);
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);
    @test o1 ≈ [-0.9999665172233384, -0.9999564684914404, -0.9999463188148598, -0.999936563671822, -0.9999262091385633]
    @test o2 ≈ [-0.9999684369099704, -0.9997674803104111, -0.9995645451149201, -0.9993600409372991, -0.9991530903064654]
    @test o3 ≈ [-0.9998686481325542, -0.9800699223805415, -0.960663237665263, -0.9416370637447824, -0.9229914006191109]
    @test o4 ≈ [-1.0001249955642058, -1.0000499920670425, -0.9999749885698793, -0.9998500107336976, -0.9997125393138017]
    o1 = grid_differentiation(f1, grid1, 1:5);
    o2 = grid_differentiation(f2, grid2, 1:5);
    o3 = grid_differentiation(f3, grid3, 1:5);
    o4 = grid_differentiation(f4, grid4, 1:5);
    @test o1 ≈ [-0.9999665172233384, -0.9999564684914404, -0.9999463188148598, -0.999936563671822, -0.9999262091385633]
    @test o2 ≈ [-0.9999684369099704, -0.9997674803104111, -0.9995645451149201, -0.9993600409372991, -0.9991530903064654]
    @test o3 ≈ [-0.9998686481325542, -0.9800699223805415, -0.960663237665263, -0.9416370637447824, -0.9229914006191109]
    @test o4 ≈ [-1.0001249955642058, -1.0000499920670425, -0.9999749885698793, -0.9998500107336976, -0.9997125393138017]
#   -----------------------------------------------------------------------------------------
    @test fdiff_interpolation_expansion_polynom(-1, 5) == [1, 1, 1, 1, 1, 1]
    polynom = fdiff_interpolation_expansion_polynom(-1, 5);
    @test fdiff_interpolation_expansion_weights(polynom) ==  [-1, 6, -15, 20, -15, 6]
    @test fdiff_interpolation_expansion_polynom(1, 5, fwd) == [1, -1, 1, -1, 1, -1]
    polynom = fdiff_interpolation_expansion_polynom(1, 5, fwd);
    fdiff_interpolation_expansion_weights(polynom, fwd, reg) == [6, -15, 20, -15, 6, -1]
    fdiff_interpolation_expansion_weights(polynom, fwd, rev) == [-1, 6, -15, 20, -15, 6]
    fdiff_interpolation_expansion_weights(polynom, bwd, reg) == [0, 3, -6, 7, -4, 1]
    fdiff_interpolation_expansion_weights(polynom, bwd, rev) == [1, -4, 7, -6, 3, 0]
    @test fdiff_interpolation_expansion_weights(1, 5, fwd, reg) == [0, 1, 0, 0, 0, 0]
    @test fdiff_interpolation_expansion_weights(1, 5, fwd, rev) == [0, 0, 0, 0, 1, 0]
    @test fdiff_interpolation_expansion_weights(-1, 5, bwd, reg) == [6, -15, 20, -15, 6, -1]
    @test fdiff_interpolation_expansion_weights(-1, 5, bwd, rev) == [-1, 6, -15, 20, -15, 6]
    @test fdiff_interpolation_expansion_weights(1, 5, fwd) == [0, 0, 0, 0, 1, 0]
    @test [fdiff_interpolation([ν^3 for ν = -5:2], ν; k=3) for ν = 1:0.5:8] == [-125.0, -91.125, -64.0, -42.875, -27.0, -15.625, -8.0, -3.375, -1.0, -0.125, 0.0, 0.125, 1.0, 3.375, 8.0]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], fwd, reg) == [-3, 15, -33, 37, -21, 5]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], fwd, rev) == [5, -21, 37, -33, 15, -3]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], bwd, rev) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], bwd, reg) == [15, -55, 85, -69, 29, -5]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5]) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion([1, -1, 1, -1], [1, 4, 9, 16], fwd) == 0
    @test fdiff_expansion([1, 1, 1, 1], [1, 4, 9, 16], bwd) == 25
    @test fdiff_expansion([1, 1, 1, 1], [1, 4, 9, 16]) == 25
    @test fdiff_differentiation_expansion_polynom(0, bwd, 3) == [0 // 1, 1 // 1, 1 // 2, 1 // 3]
    @test fdiff_differentiation_expansion_polynom(1, bwd, 3) == [0 // 1, 1 // 1, -1 // 2, -1 // 6]
    @test create_lagrange_differentiation_matrix(3) == [-11//6 3//1 -3//2 1//3; -1//3 -1//2 1//1 -1//6; 1//6 -1//1 1//2 1//3; -1//3 3//2 -3//1 11//6]
    @test trapezoidal_epw(5; rationalize=true) == [95 // 288, 317 // 240, 23 // 30, 793 // 720, 157 // 160]
    @test trapezoidal_integration([1.0, 4.0, 15.0, 40.0, 85.0, 156.0], 0.0, 5.0, [3 // 8, 7 // 6, 23 // 24]) ≈ 215.4166666
    @test create_adams_moulton_weights(3; rationalize=true) == [1 // 24, -5 // 24, 19 // 24, 3 // 8]
    @test create_adams_moulton_weights(5; rationalize=false) == [0.01875, -0.12013888888888889, 0.3347222222222222, -0.5541666666666667, 0.9909722222222223, 0.3298611111111111]
    @test fdiff_adams_moulton_expansion_coeff(20; msg=false) == -12365722323469980029//4817145976189747200000
    @test fdiff_adams_bashford_expansion_coeff(20; msg=false) == 8136836498467582599787//33720021833328230400000
    @test create_adams_bashford_weights(5; rationalize=true, devisor=false) == Rational{Int64}[-95//288, 959//480, -3649//720, 4991//720, -2641//480, 4277//1440]
    @test create_adams_bashford_weights(5; rationalize=false) == [-0.3298611111111111, 1.9979166666666666, -5.0680555555555555, 6.9319444444444445, -5.502083333333333, 2.970138888888889]
    
end


@doc raw"""
    fdiff_expansion(polynom, f, notation=bwd)

Finite difference expansion evaluated to ``k^{th}`` order for the analytical 
function ``f``, tabulated in *regular order* (growing index) at ``k+1`` positions 
on a [`Grid`](@ref). The expansion coefficients are specified by the vector 
`polynom`. By default the expansion is calculated in backward-difference notation (`bwd`). 

**Forward difference notation** (`notation = fwd`)
```math
\sum_{p=0}^{k}α_{p}Δ^{p}f[n] = F^{k} \cdot f[n:n+k],
```
where ``f[n:n+k]`` are elements of the analytical function ``f`` (tabulated in 
*forward* order) and `polynom` is the (user-supplied) expansion coefficient vector 
``α ≡ [α_0,⋯\ α_k]``.

**Backward difference notation** (`notation = bwd`)
```math
\sum_{p=0}^{k}β_{p}∇^{p}f[n] = \bar{B}^k \cdot f[n-k:n].
```
where ``f[n-k:n]`` are elements of the
analytical function ``f`` (tabulated in *forward* order) and `polynom` is the (user-supplied)
expansion coefficient vector ``β ≡ [β_0,⋯\ β_k]``.

NB. The vector `polynom` determines the order of the expansion, 
``k+1 = \rm{length}(α) = \rm{length}(β)``. The weights vectors ``F^k`` and 
``\bar{B}^k`` are *internally generated* by the function
[`fdiff_expansion_weights(polynom, notation, ordering)`](@ref).

#### Examples:
Consider the function ``f(x)=x^2`` and the expansions,
```math
\begin{aligned}
f[n-1]=(1+Δ)^{-1}&=(1-Δ+Δ^2-Δ^3+⋯)f[n]=F^{k} \cdot f[n:n+k],\\
f[n]=(1+Δ)^{-1}&=(1-Δ+Δ^2-Δ^3+⋯)f[n+1]=F^{k} \cdot f[n+1:n+k+1],
\end{aligned}
```
```math
\begin{aligned}
f[n+1]=(1-∇)^{-1}&=(1+∇+∇^2+∇^3+⋯)f[n]=\bar{B}^k \cdot f[n-k:n],\\
f[n]=(1-∇)^{-1}&=(1+∇+∇^2+∇^3+⋯)f[n-1]=\bar{B}^k \cdot f[n-k-1:n-1]
\end{aligned}
```
Note that, in these examples, to fourth order in the expansion ``(k=4)``, 
the forward- and backward-difference coefficient vectors `polynom` are 
given by ``α=[1,-1,1,-1,1]`` and ``β=[1,1,1,1,1]``, respectively. 
```
julia> f = [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100];

julia> α = [1,-1,1,-1,1];

julia> β = [1,1,1,1,1];

julia> fdiff_expansion(α,f[7:11],fwd)
25

julia> fdiff_expansion(β,f[1:5],fwd)
25

julia> f[6]
25
```
In these cases the results are exact because the function is quadratic and
the expansion is third order (based on the polynomial of ``k^{th}`` degree
running through the ``k+1`` points of the tabulated function). Compare with 
the example of [`fdiff_interpolation(f, v, k=3)`](@ref).
"""
function fdiff_expansion(polynom, f, notation=bwd)

    ordering = CamiMath.isforward(notation) ? reg : rev
    w = fdiff_expansion_weights(polynom, notation, ordering)

    return LinearAlgebra.dot(w, f)

end

# ------------------------------------------------------------------------------
#                  fdiff_interpolation_expansion_weights(polynom)
# ------------------------------------------------------------------------------

function fwd_interpolation_expansion_weights(σ::T, ordering=reg; k=3) where T<:Real

    α = fdiff_interpolation_expansion_polynom(σ, k, fwd)
    o = fdiff_expansion_weights(α, fwd, ordering)

    return o

end
#...............................................................................
function bwd_interpolation_expansion_weights(σ::T, ordering=rev; k=3) where T<:Real

    β = fdiff_interpolation_expansion_polynom(σ, k, bwd)
    o = fdiff_expansion_weights(β, bwd, ordering)

    return o

end
#...............................................................................
@doc raw"""
    fdiff_interpolation_expansion_weights(σ::T, [, notation=bwd [, ordering=rev [, k=3]]]) where T<:Real
    fdiff_interpolation_expansion_weights(polynom [, notation=bwd [, ordering=rev]])

Finite-difference expansion weights vector defining the ``k^{th}``-order
(default *third* order) Lagrange-polynomial interpolation of a tabulated
analytic function ``f[n]`` at offset ``σ`` with respect to index
position ``n``, which is positive for increasing index and negative for
decreasing index.

**Forward difference notation** (`notation = fwd`)

In this case we consider the tabulated interval ``f[n:n+k]``. The interpolated
value ``f[n-σ]`` is given by the forward-difference expansion

```math
f[n-σ] = \sum_{p=0}^k α_p(σ) Δ^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_polynom(σ, k, fwd)`](@ref)
`` → α(σ) ≡ [α_0(σ),⋯\ α_k(σ)]``. In this notation the range
``-k ≤ σ ≤ 1`` corresponds to interpolation and the ranges ``σ < -k`` and
``σ > 1k`` to extrapolation.

**Backward difference notation** (`notation = bwd`)

In this case we consider the tabulated interval ``f[n-k:n]``. The interpolated
value ``f[n+σ]`` is given by the backward-difference expansion

```math
f[n+σ] = \sum_{p=0}^k β_p(σ) ∇^p f[n] + ⋯,
```
where the expansion coefficients are given by

[`fdiff_interpolation_expansion_polynom(σ, bwd; k=3)`](@ref)
`` → β(σ) ≡ [β_0(σ),⋯\ β_k(σ)]``. In this notation the range
``-k ≤ σ ≤ 1`` corresponds to interpolation and the ranges ``σ < -k`` and
``σ > 1k`` to extrapolation.

#### Examples:
```
julia> k=4;

julia> α = fdiff_interpolation_expansion_polynom(1, k, fwd); println("α = $α")
α = [1, -1, 1, -1, 1]

julia> Fk1 = fdiff_interpolation_expansion_weights(α, fwd, reg); println("Fk1 = $(Fk1)")
Fk1 = [5, -10, 10, -5, 1]

julia> β = fdiff_interpolation_expansion_polynom(1, k, bwd); println("β = $β")
β = [1, 1, 1, 1, 1]

julia> revBk1 = fdiff_interpolation_expansion_weights(β, bwd, rev); println("revBk1 = $(revBk1)")
revBk1 = [1, -5, 10, -10, 5]
```
"""
function fdiff_interpolation_expansion_weights(σ::T, notation=bwd, ordering=rev; k=3) where T<:Real

    o = CamiMath.isforward(notation) ? fwd_interpolation_expansion_weights(σ, k, ordering) :
                                       bwd_interpolation_expansion_weights(σ, k, ordering)
    return o

end
function fdiff_interpolation_expansion_weights(polynom, notation=bwd, ordering=rev)

    if CamiMath.isforward(notation)

        o = CamiMath.isregular(ordering) ? reg_fwd_expansion_weights(polynom) :
                                           rev_fwd_expansion_weights(polynom)
    else

        o = CamiMath.isregular(ordering) ? reg_bwd_expansion_weights(polynom) :
                                           rev_bwd_expansion_weights(polynom)
    end

    return o

end

# ------------------------------------------------------------------------------
#                          fdiff_interpolation(f, x; k=3)
# ------------------------------------------------------------------------------

@doc raw"""
    fdiff_interpolation(f::Vector{T}, v::V; k=3) where {T<:Real, V<:Real}

Finite difference lagrangian interpolation (by default *third* order) at real
position v (in index units) with respect to the elements of the uniformly
tabulated analytic function `f[1:N]`. The interpolation points are situated on a
Lagrange polynomial of degree ``k`` (by default *third* degree) running through
``k+1`` subsequenct points of the tabulated function. Outside the tabulated
range, the method represents extrapolation on the lagrangian polynomial defined
by the first/last ``k+1`` tabulated points.

Beware that the interpolation becomes inaccurate if the tabulated function
cannot be approximated by a polynomial of degree ``k``.
#### Examples:
```
f = [1,2,3,4,5,6,7]
[fdiff_interpolation(f, v; k=3) for v=1:0.5:7]
  [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]

f = [1,4,9,16,25,36,49]
[fdiff_interpolation(f, v; k=3) for v=1:0.5:7]
 [1.0, 2.25, 4.0, 6.25, 9.0, 12.25, 16.0, 20.25, 25.0, 30.25, 36.0, 42.25, 49.0]

 f = [x^3 for x=-4:2]
 f1(v) = fdiff_interpolation(f, v; k=1)
 f2(v) = fdiff_interpolation(f, v; k=2)
 f3(v) = fdiff_interpolation(f, v; k=3)
 [[f1(v),f2(v),f3(v)] for v=1:0.5:9]
   [[-64.0, -64.0, -64.0], [-45.5, -43.25, -42.875], [-27.0, -27.0, -27.0],
   [-17.5, -16.0, -15.625], [-8.0, -8.0, -8.0], [-4.5, -3.75, -3.375],
   [-1.0, -1.0, -1.0], [-0.5, -0.5, -0.125], [0.0, 0.0, 0.0],
   [0.5, -0.25, 0.125], [1.0, 1.0, 1.0], [4.5, 3.75, 3.375], [8.0, 8.0, 8.0],
   [11.5, 13.75, 15.625], [15.0, 21.0, 27.0], [18.5, 29.75, 42.875],
   [22.0, 40.0, 64.0]]
```
The result for f3(v) is exact because the function is cubic and
the expansion is third order - see Figure below. The tabulated function is
given by the black points. The interpolation and extrapolation points are red.

![Image](../assets/lagrangian_interpolation.png)
"""
function fdiff_interpolation(f::Vector{T}, v::V; k=4) where {T<:Real, V<:Real}

    l = length(f)
    k = min(k,l-1)
    k > 0 || error("Error: k ≥ 1 required for lagrangian interpolation")
    n = v < 1 ? 1 : v < l-k ? floor(Int,v) : l-k
    α = fdiff_interpolation_expansion_polynom(n-v, k, fwd)
    o = fdiff_expansion(α, f[n:n+k], fwd)

    return o

end





function sdot(x::Vector{T}, y::Vector{T}) where T<:Real

    n = length(x)

    n == length(y) || error("Error: vectors must have the same length")

    s = T(0)

    for i = 1:n
        s += x[i] * y[i]
    end

    return s

end
function sdot1(x::Vector{T}, y::Vector{T}) where T<:Real

    n = length(x)

    n == length(y) || error("Error: vectors must have the same length")

    s1 = s2 = T(0)

    for i = 1:n
        if (x[i] < 0) ⊻ (y[i] < 0)
            s1 += x[i] * y[i]
        else
            s2 += x[i] * y[i]
        end
    end

    return s1+s2

end