# SPDX-License-Identifier: MIT

#using CamiDiff
using CamiMath

# using IntervalSets
#using BenchmarkTools
using LinearAlgebra
using Test

@testset "CamiXon.jl" begin 

    grid2 = castGrid1(2, 4, Float64; h = 0.1, r0 = 2.0, p=1);
    grid3 = castGrid1(3, 4, Float64; h = 0.1, r0 = 2.0);
    grid4 = castGrid1(4, 4, Float64; h = 0.1, r0 = 2.0, polynom=[0, 1]);
    @test [grid2.r, grid2.r′, grid2.r′′] ≈ [grid3.r, grid3.r′, grid3.r′′] ≈ [grid4.r, grid4.r′, grid4.r′′]
    
end