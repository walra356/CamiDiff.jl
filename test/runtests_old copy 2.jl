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



 
#   -----------------------------------------------------------------------------------------
    gaussian(r) = sqrt(2.0/π) * exp(-r^2/2.0);
    #grid1 = castGrid(1, 1000, Float64; h = 0.005, r0 = 0.1, msg=false);
    grid1 = castGrid(1, 1000, Float64; h = 0.005, r0 = 0.06137, msg=true);
    r1 = grid1.r;
    f1 = [gaussian(r1[n]) for n=1:grid1.N];
#   -----------------------------------------------------------------------------------------
    o1 = grid_integration(f1, grid1);
    @test o1 ≈ 1.0
    o1 = grid_integration(f1, grid1, 1:950);
    @test o1 ≈ 1.0
#   -----------------------------------------------------------------------------------------
    f′1 = -r1 .* f1;
    o1 = grid_differentiation(f1, grid1);
    @test f′1 ≈ o1
    a1 = f′1 ./ o1
    println("a1: ", a1[993:1000])
    println("f′1: ", f′1[993:1000])
    println("o1: ", o1[993:1000])