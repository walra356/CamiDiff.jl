# SPDX-License-Identifier: MIT

using CamiDiff
using CamiMath

# using IntervalSets
#using BenchmarkTools
using LinearAlgebra
using Test

@testset "CamiXon.jl" begin 
    @test fdiff_interpolation_expansion_coeffs(-1, 5) == [1, 1, 1, 1, 1, 1]
    coeffs = fdiff_interpolation_expansion_coeffs(-1, 5);
    @test fdiff_interpolation_expansion_weights(coeffs) ==  [-1, 6, -15, 20, -15, 6]
    @test fdiff_interpolation_expansion_coeffs(1, 5, fwd) == [1, -1, 1, -1, 1, -1]
    coeffs = fdiff_interpolation_expansion_coeffs(1, 5, fwd);
    fdiff_interpolation_expansion_weights(coeffs, fwd, reg) == [6, -15, 20, -15, 6, -1]
    fdiff_interpolation_expansion_weights(coeffs, fwd, rev) == [-1, 6, -15, 20, -15, 6]
    fdiff_interpolation_expansion_weights(coeffs, bwd, reg) == [0, 3, -6, 7, -4, 1]
    fdiff_interpolation_expansion_weights(coeffs, bwd, rev) == [1, -4, 7, -6, 3, 0]
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
    @test fdiff_differentiation_expansion_coeffs(0, 3) == [0 // 1, 1 // 1, 1 // 2, 1 // 3]
    @test fdiff_differentiation_expansion_coeffs(1, 3) == [0 // 1, 1 // 1, -1 // 2, -1 // 6]
    @test [fdiff_differentiation([16, 9, 4, 1, 0, 1, 4, 9, 16], v) for v = 1:9] == [-8 // 1, -6 // 1, -4 // 1, -2 // 1, 0 // 1, 2 // 1, 4 // 1, 6 // 1, 8 // 1]
    @test fdiff_differentiation([16, 9, 4, 1, 0, 1, 4, 9, 16], 5.5) == 1.0
    @test create_lagrange_differentiation_matrix(3) == [-11//6 3//1 -3//2 1//3; -1//3 -1//2 1//1 -1//6; 1//6 -1//1 1//2 1//3; -1//3 3//2 -3//1 11//6]
    @test trapezoidal_epw(5; rationalize=true) == [95 // 288, 317 // 240, 23 // 30, 793 // 720, 157 // 160]
    @test trapezoidal_integration([1.0, 4.0, 15.0, 40.0, 85.0, 156.0], 0.0, 5.0, [3 // 8, 7 // 6, 23 // 24]) ≈ 215.4166666
    @test create_adams_moulton_weights(3; rationalize=true) == [1 // 24, -5 // 24, 19 // 24, 3 // 8]
end
