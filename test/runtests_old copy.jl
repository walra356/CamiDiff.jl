using CamiDiff
using Test # to use runtests with @testset 
using LinearAlgebra

import CamiMath
fwd = CamiMath.fwd
bwd = CamiMath.bwd
reg = CamiMath.reg
rev = CamiMath.rev

grid = castGrid(1, 1000, Float64; h = 0.005, rmax = 10, msg=true);
grid = castGrid("exponential", 1000, Float64; h = 0.005, rmax = 10, msg=true);
grid = castGrid(2, 1000, Float64; h = 0.005, rmax = 10, p=5, msg=true);
grid = castGrid(3, 1000, Float64; h = 0.1, rmax = 10, msg=true);
grid = castGrid(4, 1000, Float64; h = 0.1, rmax = 10, polynom=[0,0,1], msg=true);
r = grid.r[1:4]; println("r = ", r)
r′= grid.r′[1:4]; println("r′ = ", r′)
r′′= grid.r′′[1:4]; println("r′′ = ", r′′)