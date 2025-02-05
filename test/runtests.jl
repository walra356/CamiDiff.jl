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
using LinearAlgebra

import CamiMath
fwd = CamiMath.fwd
bwd = CamiMath.bwd
reg = CamiMath.reg
rev = CamiMath.rev

@testset "CamiDiff.jl" begin 

    println("CamiDiff.jl  | 153 runtests | runtime 11.2s (estimated) | start")

    @test _gridspecs(1, 1000, Float64, 0.01, 2.0, 43612.6) == "Grid: exponential, Float64, rmax = 43612.6, Ntot = 1000, h = 0.01, r0 = 2.0"
    @test _gridspecs(2, 1000, Float64, 0.01, 2.0, 2940.47) == "Grid: truncated-exponential, Float64, rmax = 2940.47, Ntot = 1000, p = 5, h = 0.01, r0 = 2.0"
    @test _gridspecs(3, 1000, Float64, 0.1, 2.0, 199.8) == "Grid: linear (uniform), Float64, rmax = 199.8, Ntot = 1000, p = 1, h = 0.1, r0 = 2.0"
    @test _gridspecs(4, 1000, Float64, 0.1, 2.0, 199.8) == "Grid: polynomial of degree 1, Float64, rmax = 199.8, Ntot = 1000, polynom = [0, 1], h = 0.1, r0 = 2.0"
    
    grid1 = castGrid(1, 1000, Float64; h = 0.01, rmax = 2.0, msg=false);
    grid2 = castGrid(2, 1000, Float64; h = 0.01, rmax = 2.0, p=5, msg=false);
    grid3 = castGrid(3, 1000, Float64; h = 0.1, rmax = 2.0, msg=false);
    grid4 = castGrid("polynomial", 1000, Float64; h = 0.1, rmax = 2.0, polynom=[0, 1], msg=false);

    grid1 = castGrid(1, 4, Float64; h = 0.1, rmax = 2.0, msg=false);
    grid2 = castGrid(2, 4, Float64; h = 0.1, rmax = 2.0);
    grid3 = castGrid(3, 4, Float64; h = 0.1, rmax = 2.0);
    grid4 = castGrid("polynomial", 4, Float64; h = 0.1, rmax = 2.0, polynom=[0, 1]);
    @test [grid3.r, grid3.r′, grid3.r′′] ≈ [grid4.r, grid4.r′, grid4.r′′]
    @test [grid1.r, grid1.r′, grid1.r′′] == [[0.0, 0.6012192107114549, 1.265669197778149, 2.0], [0.5716591827020164, 0.6317811037731619, 0.6982261024798313, 0.7716591827020164], [0.057165918270201635, 0.06317811037731619, 0.06982261024798313, 0.07716591827020164]]
    @test grid1.name == "exponential"
    n = gridPos(1.0, grid1);
    @test n == 2
    Δn = fracPos(n, 1.0, grid1);
    t = (n+Δn-1)*grid1.h;
    @test grid1.r0 * (exp(t)-1) ≈ 1.0
    
    @test gridtypename(1) == "exponential"
    @test gridtypeID("exponential") == 1
    @test_throws DomainError gridPos(100.0, grid1) == 220
    @test_throws DomainError castGrid(5, 1000, Float64)
    @test_throws DomainError gridtypename(5) 
    @test_throws DomainError gridtypeID("logarithm") 
#   =========================================================================================
    gaussian(r) = sqrt(2.0/π) * exp(-r^2/2.0);  
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, 1000, Float64; h = 0.005, rmax=9.0, msg=false);
    grid2 = castGrid(2, 1000, Float64; h = 0.005, rmax=9.0, p=5, msg=false);
    grid3 = castGrid(3, 1000, Float64; h = 0.1, rmax=9.0, msg=false);
    grid4 = castGrid(4, 1000, Float64; h = 0.1, rmax=9.0, polynom=[0,0,1], msg=false);
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
    o1 = grid_integration(f1, grid1, 1:950);
    o2 = grid_integration(f2, grid2, 1:950);
    o3 = grid_integration(f3, grid3, 1:950);
    o4 = grid_integration(f4, grid4, 1:950);
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
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);
    @test f′1 ≈ o1
    @test f′2 ≈ o2
    @test f′3 ≈ o3
    @test f′4 ≈ o4  
#   -----------------------------------------------------------------------------------------
    o1 = grid_differentiation(f1, grid1, 1:900);
    o2 = grid_differentiation(f2, grid2, 1:900);
    o3 = grid_differentiation(f3, grid3, 1:900);
    o4 = grid_differentiation(f4, grid4, 1:900);
    @test f′1[1:900] ≈ o1
    @test f′2[1:900] ≈ o2
    @test f′3[1:900] ≈ o3
    @test f′4[1:900] ≈ o4  
#   -----------------------------------------------------------------------------------------
    o1fwd = grid_differentiation(f1, grid1, 500, fwd);
    o1bwd = grid_differentiation(f1, grid1, 500, bwd);
    o2fwd = grid_differentiation(f2, grid2, 500, fwd);
    o2bwd = grid_differentiation(f2, grid2, 500, bwd);
    o3fwd = grid_differentiation(f3, grid3, 500, fwd);
    o3bwd = grid_differentiation(f3, grid3, 500, bwd); 
    o4fwd = grid_differentiation(f4, grid4, 500, fwd);
    o4bwd = grid_differentiation(f4, grid4, 500, bwd);
    @test f′1[500] ≈ o1fwd ≈ o1bwd
    @test f′2[500] ≈ o2fwd ≈ o2bwd
    @test f′3[500] ≈ o3fwd ≈ o3bwd
    @test f′4[500] ≈ o4fwd ≈ o4bwd
#   -----------------------------------------------------------------------------------------
    r = 1.0;
    f = gaussian(r);
    f′1 = -r * gaussian(r)
    @test grid_interpolation(f1, grid1, r, fwd) ≈ f
    @test grid_interpolation(f2, grid2, r, bwd) ≈ f
    @test grid_interpolation(f3, grid3, r, fwd) ≈ f
    @test grid_interpolation(f4, grid4, r, bwd) ≈ f
    @test grid_differentiation(f1, grid1, r, fwd) ≈ f′1
    @test grid_differentiation(f2, grid2, r, bwd) ≈ f′1
    @test grid_differentiation(f3, grid3, r, fwd) ≈ f′1
    @test grid_differentiation(f4, grid4, r, bwd) ≈ f′1
#   ========================================================================================= 
#   f = exp(-r)  
#   -----------------------------------------------------------------------------------------
        grid1 = castGrid(1, 9, Float64; h = 0.01, rmax=25, msg=false);
        grid2 = castGrid(2, 9, Float64; h = 0.01, rmax=25, p=6, msg=false);
        grid3 = castGrid(3, 9, Float64; h = 0.01, rmax=25, msg=false);
        grid4 = castGrid(4, 9, Float64; h = 0.01, rmax=25, polynom=[0,0,1], msg=false);
        r1 = grid1.r;
        r2 = grid2.r;
        r3 = grid3.r;
        r4 = grid4.r;
        f1 = [exp(-r1[n]) for n=1:grid1.N];
        f2 = [exp(-r2[n]) for n=1:grid2.N];
        f3 = [exp(-r3[n]) for n=1:grid3.N];
        f4 = [exp(-r4[n]) for n=1:grid4.N];
    #   -----------------------------------------------------------------------------------------
        o1 = grid_integration(f1, grid1); println("o1: ", o1)
        o2 = grid_integration(f2, grid2); println("o2: ", o2)
        o3 = grid_integration(f3, grid3); println("o3: ", o3)
        o4 = grid_integration(f4, grid4); println("o4: ", o4) 
    #    @test o1 ≈ 1.0
    #    @test o2 ≈ 1.0
    #    @test o3 ≈ 1.0
    #    @test o4 ≈ 1.0
#   ========================================================================================= 
#   f = exp(-r)  
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, 1000, Float64; h = 0.01, rmax=25, msg=false);
    grid2 = castGrid(2, 1000, Float64; h = 0.01, rmax=25, p=6, msg=false);
    grid3 = castGrid(3, 1000, Float64; h = 0.01, rmax=25, msg=false);
    grid4 = castGrid(4, 1000, Float64; h = 0.01, rmax=25, polynom=[0,0,1], msg=false);
    r1 = grid1.r;
    r2 = grid2.r;
    r3 = grid3.r;
    r4 = grid4.r;
    f1 = [exp(-r1[n]) for n=1:grid1.N];
    f2 = [exp(-r2[n]) for n=1:grid2.N];
    f3 = [exp(-r3[n]) for n=1:grid3.N];
    f4 = [exp(-r4[n]) for n=1:grid4.N];
#   -----------------------------------------------------------------------------------------
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
#   -----------------------------------------------------------------------------------------
    o1 = grid_differentiation(f1, grid1, 1:900);
    o2 = grid_differentiation(f2, grid2, 1:900);
    o3 = grid_differentiation(f3, grid3, 1:900);
    o4 = grid_differentiation(f4, grid4, 1:900);
    @test f′1[1:900] ≈ o1
    @test f′2[1:900] ≈ o2
    @test f′3[1:900] ≈ o3
    @test f′4[1:900] ≈ o4  
#   -----------------------------------------------------------------------------------------
    o1fwd = grid_differentiation(f1, grid1, 500, fwd);
    o1bwd = grid_differentiation(f1, grid1, 500, bwd);
    o2fwd = grid_differentiation(f2, grid2, 500, fwd);
    o2bwd = grid_differentiation(f2, grid2, 500, bwd);
    o3fwd = grid_differentiation(f3, grid3, 500, fwd);
    o3bwd = grid_differentiation(f3, grid3, 500, bwd); 
    o4fwd = grid_differentiation(f4, grid4, 500, fwd);
    o4bwd = grid_differentiation(f4, grid4, 500, bwd);
    @test f′1[500] ≈ o1fwd ≈ o1bwd
    @test f′2[500] ≈ o2fwd ≈ o2bwd
    @test f′3[500] ≈ o3fwd ≈ o3bwd
    @test f′4[500] ≈ o4fwd ≈ o4bwd
#   -----------------------------------------------------------------------------------------
    r = 1.0;
    @test grid_interpolation(f1, grid1, r, fwd) ≈ exp(-r)
    @test grid_interpolation(f2, grid2, r, bwd) ≈ exp(-r)
    @test grid_interpolation(f3, grid3, r, fwd) ≈ exp(-r)
    @test grid_interpolation(f4, grid4, r, bwd) ≈ exp(-r)
    @test grid_differentiation(f1, grid1, r, fwd) ≈ -exp(-r)
    @test grid_differentiation(f2, grid2, r, bwd) ≈ -exp(-r)
    @test grid_differentiation(f3, grid3, r, fwd) ≈ -exp(-r)
    @test grid_differentiation(f4, grid4, r, bwd) ≈ -exp(-r)
#   ========================================================================================= 
#   f = exp(-r)
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, 250, Float64; h = 0.01, rmax=2.5, msg=false);
    grid2 = castGrid(2, 250, Float64; h = 0.01, rmax=2.5, p=6, msg=false);
    grid3 = castGrid(3, 250, Float64; h = 0.01, rmax=2.5, msg=false);
    grid4 = castGrid(4, 250, Float64; h = 0.01, rmax=2.5, polynom=[0,0,1], msg=false);
    r1 = grid1.r;
    r2 = grid2.r;
    r3 = grid3.r;
    r4 = grid4.r;
    f1 = [exp(-r1[n]) for n=1:grid1.N];
    f2 = [exp(-r2[n]) for n=1:grid2.N];
    f3 = [exp(-r3[n]) for n=1:grid3.N];
    f4 = [exp(-r4[n]) for n=1:grid4.N];
    f′1 = -f1;
    f′2 = -f2;
    f′3 = -f3;
    f′4 = -f4;
#   -----------------------------------------------------------------------------------------
    n1 = 1;
    n2 = 245;
    o1 = grid_differentiation(f1, grid1, n1:n2);
    o2 = grid_differentiation(f2, grid2, n1:n2);
    o3 = grid_differentiation(f3, grid3, n1:n2);
    o4 = grid_differentiation(f4, grid4, n1:n2);
    @test f′1[n1:n2] ≈ o1
    @test f′2[n1:n2] ≈ o2
    @test f′3[n1:n2] ≈ o3
    @test f′4[n1:n2] ≈ o4 
#   -----------------------------------------------------------------------------------------
    n1 = 6;
    n2 = 250;
    o1 = grid_differentiation(f1, grid1, n1:n2);
    o2 = grid_differentiation(f2, grid2, n1:n2);
    o3 = grid_differentiation(f3, grid3, n1:n2);
    o4 = grid_differentiation(f4, grid4, n1:n2);
    @test f′1[n1:n2] ≈ o1
    @test f′2[n1:n2] ≈ o2
    @test f′3[n1:n2] ≈ o3
    @test f′4[n1:n2] ≈ o4
#   -----------------------------------------------------------------------------------------
    n1 = 6;
    n2 = 245;
    o1 = grid_differentiation(f1, grid1, n1:n2);
    o2 = grid_differentiation(f2, grid2, n1:n2);
    o3 = grid_differentiation(f3, grid3, n1:n2);
    o4 = grid_differentiation(f4, grid4, n1:n2);
    @test f′1[n1:n2] ≈ o1
    @test f′2[n1:n2] ≈ o2
    @test f′3[n1:n2] ≈ o3
    @test f′4[n1:n2] ≈ o4
#   ========================================================================================= 
    linear(r) = 1.00r;
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, 1000, Float64; h = 0.001, rmax=2, msg=false);
    grid2 = castGrid(2, 1000, Float64; h = 0.001, rmax=2, p=5, msg=false);
    grid3 = castGrid(3, 1000, Float64; h = 0.01, rmax=2, msg=false);
    grid4 = castGrid(4, 1000, Float64; h = 0.01, rmax=2, polynom=[0,0,1], msg=false);
    r1 = grid1.r
    r2 = grid2.r
    r3 = grid3.r
    r4 = grid4.r
    f1 = [linear(r1[n]) for n ∈ eachindex(r1)];
    f2 = [linear(r2[n]) for n ∈ eachindex(r2)]; 
    f3 = [linear(r3[n]) for n ∈ eachindex(r3)]; 
    f4 = [linear(r4[n]) for n ∈ eachindex(r4)]; 
#   -----------------------------------------------------------------------------------------
    o1 = grid_integration(f1, grid1);
    o2 = grid_integration(f2, grid2);
    o3 = grid_integration(f3, grid3);
    o4 = grid_integration(f4, grid4);
    @test o1 ≈ 2.0
    @test o2 ≈ 2.0
    @test o3 ≈ 2.0
    @test o4 ≈ 2.0
#   -----------------------------------------------------------------------------------------
    f′ = [1.00 for n ∈ eachindex(r1)];
    o1 = grid_differentiation(f1, grid1);
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);
    @test o1 ≈ f′
    @test o2 ≈ f′
    @test o3 ≈ f′
    @test o4 ≈ f′
#   -----------------------------------------------------------------------------------------
    r = 1.0;
    @test grid_interpolation(f1, grid1, r, fwd) ≈ r
    @test grid_interpolation(f2, grid2, r, bwd) ≈ r
    @test grid_interpolation(f3, grid3, r, fwd) ≈ r
    @test grid_interpolation(f4, grid4, r, bwd) ≈ r
    @test grid_differentiation(f1, grid1, r, fwd) ≈ 1.0
    @test grid_differentiation(f2, grid2, r, bwd) ≈ 1.0
    @test grid_differentiation(f3, grid3, r, fwd) ≈ 1.0
    @test grid_differentiation(f4, grid4, r, bwd) ≈ 1.0
#   ========================================================================================= 
    T = Float64
    N=1000
    k = 5
    a = 1.001
    a = T == BigFloat ? T(rationalize(a)) : a
    line(r) = T(a*r);
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, N, T; h = 0.01, rmax=2, msg=false);
    grid2 = castGrid(2, N, T; h = 0.01, rmax=2, p=4, msg=false);
    grid3 = castGrid(3, N, T; h = 0.01, rmax=2, msg=false);
    grid4 = castGrid(4, N, T; h = 0.01, rmax=2, polynom=[0,0,1], msg=false);
    r1 = grid1.r
    r2 = grid2.r
    r3 = grid3.r
    r4 = grid4.r
    f1 = [line(r1[n]) for n ∈ eachindex(r1)];
    f2 = [line(r2[n]) for n ∈ eachindex(r2)]; 
    f3 = [line(r3[n]) for n ∈ eachindex(r3)]; 
    f4 = [line(r4[n]) for n ∈ eachindex(r4)]; 
#   -----------------------------------------------------------------------------------------
    int_f1 = a* T(1//2)*r1[N]^2
    int_f2 = a* T(1//2)*r2[N]^2
    int_f3 = a* T(1//2)*r3[N]^2
    int_f4 = a* T(1//2)*r4[N]^2
    o1 = grid_integration(f1, grid1);
    o2 = grid_integration(f2, grid2);
    o3 = grid_integration(f3, grid3);
    o4 = grid_integration(f4, grid4);
    @test o1 ≈ int_f1; # println("o1 ", o1)
    @test o2 ≈ int_f2; # println("o2 ", o2)
    @test o3 ≈ int_f3; # println("o3 ", o3)
    @test o4 ≈ int_f4; # println("o4 ", o4)
#   -----------------------------------------------------------------------------------------
    f′ = [1.001 for n=1:N];
    o1 = grid_differentiation(f1, grid1; k); # println("o1 ", o1[N÷2])
    o2 = grid_differentiation(f2, grid2; k); # println("o2 ", o2[N÷2])
    o3 = grid_differentiation(f3, grid3; k); # println("o3 ", o3[N÷2])
    o4 = grid_differentiation(f4, grid4; k); # println("o4 ", o4[N÷2])
    @test o1 ≈ f′
    @test o2 ≈ f′
    @test o3 ≈ f′
    @test o4 ≈ f′
#   ========================================================================================= 
    T = Float64
    N=1000
    k = 5
    a = 1.001
    a = T == BigFloat ? T(rationalize(a)) : a
    curve(r) = T(1//2)*T(a)*T(r)*T(r);
#   -----------------------------------------------------------------------------------------
    grid1 = castGrid(1, N, T; h = 0.01, rmax=2, msg=false);
    grid2 = castGrid(2, N, T; h = 0.01, rmax=2, p=4, msg=false);
    grid3 = castGrid(3, N, T; h = 0.01, rmax=2, msg=false);
    grid4 = castGrid(4, N, T; h = 0.01, rmax=2, polynom=[0,0,1], msg=false);
    r1 = grid1.r
    r2 = grid2.r
    r3 = grid3.r
    r4 = grid4.r
    f1 = [curve(r1[n]) for n ∈ eachindex(r1)];
    f2 = [curve(r2[n]) for n ∈ eachindex(r2)]; 
    f3 = [curve(r3[n]) for n ∈ eachindex(r3)]; 
    f4 = [curve(r4[n]) for n ∈ eachindex(r4)]; 
#   -----------------------------------------------------------------------------------------
    int_f1 = a* T(1//6)*r1[N]^3
    int_f2 = a* T(1//6)*r2[N]^3
    int_f3 = a* T(1//6)*r3[N]^3
    int_f4 = a* T(1//6)*r4[N]^3
    o1 = grid_integration(f1, grid1);
    o2 = grid_integration(f2, grid2);
    o3 = grid_integration(f3, grid3);
    o4 = grid_integration(f4, grid4);
    @test o1 ≈ int_f1; # println("o1 ", o1)
    @test o2 ≈ int_f2; # println("o2 ", o2)
    @test o3 ≈ int_f3; # println("o3 ", o3)
    @test o4 ≈ int_f4; # println("o4 ", o4)
#   -----------------------------------------------------------------------------------------        
    f′1 = [line(r1[n]) for n ∈ eachindex(r1)];
    f′2 = [line(r2[n]) for n ∈ eachindex(r2)]; 
    f′3 = [line(r3[n]) for n ∈ eachindex(r3)]; 
    f′4 = [line(r4[n]) for n ∈ eachindex(r4)]; 
    o1 = grid_differentiation(f1, grid1; k); # println("o1 ", o1[N÷2])
    o2 = grid_differentiation(f2, grid2; k); # println("o2 ", o2[N÷2])
    o3 = grid_differentiation(f3, grid3; k); # println("o3 ", o3[N÷2])
    o4 = grid_differentiation(f4, grid4; k); # println("o4 ", o4[N÷2])
    @test o1 ≈ f′1
    @test o2 ≈ f′2
    @test o3 ≈ f′3
    @test o4 ≈ f′4
#   -----------------------------------------------------------------------------------------
    f = [v^3 for v=1:10];
#   -----------------------------------------------------------------------------------------
    n = 5; k=5; 
    v = 4; σ = n-v # fwd offset
    polynom = fdiff_interpolation_expansion_polynom(σ, k, fwd)
    Fk = fdiff_expansion_weights(polynom, fwd, reg) 
    @test polynom == [1, -1, 1, -1, 1, -1] 
    @test Fk ⋅ f[n:n+k] ≈ v^3

    v = 5; σ = n-v # fwd offset
    polynom = fdiff_interpolation_expansion_polynom(σ, k, fwd) 
    Fk = fdiff_expansion_weights(polynom, fwd, reg) 
    @test polynom == [1, 0, 0, 0, 0, 0]
    @test Fk ⋅ f[n:n+k] ≈ v^3 

    v = 7.5; σ = n-v # fwd offset
    polynom = fdiff_interpolation_expansion_polynom(σ, k, fwd)
    Fk = fdiff_expansion_weights(polynom, fwd, reg)  
    @test polynom == [1.0, 2.5, 1.875, 0.3125, -0.0390625, 0.01171875]
    @test Fk ⋅ f[n:n+k] ≈ v^3
#   -----------------------------------------------------------------------------------------
    n = 8; k=5; 
    v = 9; σ = -(n-v) # bwd offset
    polynom = fdiff_interpolation_expansion_polynom(σ, k, bwd)
    revBk = fdiff_expansion_weights(polynom, bwd, rev)  
    @test polynom == [1, 1, 1, 1, 1, 1]
    @test revBk ⋅ f[n-k:n] ≈ v^3
     
    v = 8; σ = -(n-v) # bwd offset
    polynom = fdiff_interpolation_expansion_polynom(σ, k, bwd)
    revBk = fdiff_expansion_weights(polynom, bwd, rev)   
    @test polynom == [1, 0, 0, 0, 0, 0]
    @test revBk ⋅ f[n-k:n] ≈ v^3
     
    v = 5.5; σ = -(n-v) # bwd offset
    polynom = fdiff_interpolation_expansion_polynom(σ, k, bwd)
    revBk = fdiff_expansion_weights(polynom, bwd, rev)   
    @test polynom == [1.0, -2.5, 1.875, -0.3125, -0.0390625, -0.01171875]
    @test revBk ⋅ f[n-k:n] ≈ v^3
    #   -----------------------------------------------------------------------------------------
    n = 5; k=5; 
    v = 4; σ = n-v # fwd offset
    polynom = fdiff_differentiation_expansion_polynom(σ, k, fwd) 
    Fk = fdiff_expansion_weights(polynom, fwd, reg)  
    @test polynom == [0, 1, -3//2, 11//6, -25//12, 137//60]
    @test Fk ⋅ f[n:n+k] ≈ 3v^2
 
    v = 5; σ = n-v # fwd offset
    polynom = fdiff_differentiation_expansion_polynom(σ, k, fwd)
    Fk = fdiff_expansion_weights(polynom, fwd, reg)   
    @test polynom == [0, 1, -1//2, 1//3, -1//4, 1//5]
    @test Fk ⋅ f[n:n+k] ≈ 3v^2
 
    v = 7.5; σ = n-v # fwd offset
    polynom = fdiff_differentiation_expansion_polynom(σ, k, fwd) 
    Fk = fdiff_expansion_weights(polynom, fwd, reg)  
    @test polynom == [0.0, 1.0, 2.0, 0.9583333333333333, -0.04166666666666674, 0.004687500000000011]
    @test Fk ⋅ f[n:n+k] ≈ 3v^2
    #   -----------------------------------------------------------------------------------------n = 8; k=5; 
    n = 8; k=5; 
    v = 9; σ = -(n-v) # bwd offset
    polynom = fdiff_differentiation_expansion_polynom(σ, k, bwd)
    revBk = fdiff_expansion_weights(polynom, bwd, rev)   
    @test polynom == [0, 1, 3//2, 11//6, 25//12, 137//60]
    @test revBk ⋅ f[n-k:n] ≈ 3v^2
 
    v = 8; σ = -(n-v) # bwd offset
    polynom = fdiff_differentiation_expansion_polynom(σ, k, bwd)
    revBk = fdiff_expansion_weights(polynom, bwd, rev)   
    @test polynom == [0, 1, 1//2, 1//3, 1//4, 1//5]
    @test revBk ⋅ f[n-k:n] ≈ 3v^2
 
    v = 5.5; σ = -(n-v) # bwd offset
    polynom = fdiff_differentiation_expansion_polynom(σ, k, bwd) 
    revBk = fdiff_expansion_weights(polynom, bwd, rev)  
    @test polynom == [0.0, 1.0, -2.0, 0.9583333333333333, 0.04166666666666674, 0.004687500000000011]
    @test revBk ⋅ f[n-k:n] ≈ 3v^2
#   -----------------------------------------------------------------------------------------
    
  #   -----------------------------------------------------------------------------------------
  #  @test [fdiff_interpolation([ν^3 for ν = -5:2], ν; k=3) for ν = 1:0.5:8] == [-125.0, -91.125, -64.0, -42.875, -27.0, -15.625, -8.0, -3.375, -1.0, -0.125, 0.0, 0.125, 1.0, 3.375, 8.0]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], fwd, reg) == [-3, 15, -33, 37, -21, 5]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], fwd, rev) == [5, -21, 37, -33, 15, -3]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], bwd, rev) == [-5, 29, -69, 85, -55, 15]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5], bwd, reg) == [15, -55, 85, -69, 29, -5]
    @test fdiff_expansion_weights([0, 1, 2, 3, 4, 5]) == [-5, 29, -69, 85, -55, 15]
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