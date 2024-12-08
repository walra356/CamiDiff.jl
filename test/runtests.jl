# SPDX-License-Identifier: MIT

using CamiDiff
using CamiMath

# using BenchmarkTools
using LinearAlgebra
using Test

@testset "CamiDiff.jl" begin 

    grid1 = castGrid(1, 4, Float64; h = 0.1, r0 = 2.0);
    grid2 = castGrid(2, 4, Float64; h = 0.1, r0 = 2.0, p=1);
    grid3 = castGrid(3, 4, Float64; h = 0.1, r0 = 2.0);
    grid4 = castGrid("polynomial", 4, Float64; h = 0.1, r0 = 2.0, polynom=[0, 1], msg=true);
    @test [grid2.r, grid2.r′, grid2.r′′] ≈ [grid3.r, grid3.r′, grid3.r′′] ≈ [grid4.r, grid4.r′, grid4.r′′]
    @test [grid1.r, grid1.r′, grid1.r′′] == [[0.0, 0.21034183615129542, 0.4428055163203397, 0.6997176151520064], [0.2, 0.22103418361512955, 0.244280551632034, 0.26997176151520064], [0.020000000000000004, 0.022103418361512958, 0.024428055163203403, 0.02699717615152007]]
#    @test [grid1.r, grid1.r′, grid1.r′′] == [[2.220446049250313e-16, 0.21034183615129542, 0.4428055163203397, 0.6997176151520064], [0.2, 0.22103418361512955, 0.244280551632034, 0.26997176151520064], [0.020000000000000004, 0.022103418361512958, 0.024428055163203403, 0.02699717615152007]]
    @test grid1.name == "exponential"
    @test findIndex(0.25, grid1) == 2
    @test_throws DomainError findIndex(100.0, grid1) == 220
    @test_throws DomainError castGrid(5, 1000, Float64)
    @test_throws DomainError gridname(5) 
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
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);
    @test f′1 ≈ o1
    @test f′2 ≈ o2
    @test f′3 ≈ o3
    @test f′4 ≈ o4
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
    grid2 = castGrid(2, 1000, Float64; h = 0.01, r0 = 0.02, p=5, msg=true);
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
    o1 = grid_differentiation(f1, grid1, 1:900);
    o2 = grid_differentiation(f2, grid2, 1:900);
    o3 = grid_differentiation(f3, grid3, 1:900);
    o4 = grid_differentiation(f4, grid4, 1:900);
    @test f′1[1:900] ≈ o1
    @test f′2[1:900] ≈ o2
    @test f′3[1:900] ≈ o3
    @test f′4[1:900] ≈ o4
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
    @test fdiff_differentiation_expansion_polynom(0, 3) == [0 // 1, 1 // 1, 1 // 2, 1 // 3]
    @test fdiff_differentiation_expansion_polynom(1, 3) == [0 // 1, 1 // 1, -1 // 2, -1 // 6]
    @test [fdiff_differentiation([16, 9, 4, 1, 0, 1, 4, 9, 16], v) for v = 1:9] == [-8 // 1, -6 // 1, -4 // 1, -2 // 1, 0 // 1, 2 // 1, 4 // 1, 6 // 1, 8 // 1]
    @test fdiff_differentiation(Real[16, 9, 4, 1, 0, 1, 4, 9, 16], 5.5) == 1.0
    @test create_lagrange_differentiation_matrix(3) == [-11//6 3//1 -3//2 1//3; -1//3 -1//2 1//1 -1//6; 1//6 -1//1 1//2 1//3; -1//3 3//2 -3//1 11//6]
    @test trapezoidal_epw(5; rationalize=true) == [95 // 288, 317 // 240, 23 // 30, 793 // 720, 157 // 160]
    @test trapezoidal_integration([1.0, 4.0, 15.0, 40.0, 85.0, 156.0], 0.0, 5.0, [3 // 8, 7 // 6, 23 // 24]) ≈ 215.4166666
    @test create_adams_moulton_weights(3; rationalize=true) == [1 // 24, -5 // 24, 19 // 24, 3 // 8]
    @test fdiff_adams_moulton_expansion_coeff(20) == -12365722323469980029//4817145976189747200000
    @test fdiff_adams_bashford_expansion_coeff(20) == 8136836498467582599787//33720021833328230400000
    
end