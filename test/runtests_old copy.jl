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
    grid1 = castGrid(1, 1000, Float64; h = 0.005, r0 = 0.06137, msg=false);
    grid2 = castGrid(2, 1000, Float64; h = 0.005, r0 = 0.0999, p=5, msg=false);
    grid3 = castGrid(3, 1000, Float64; h = 0.1, r0 = 0.0901, msg=false);
    grid4 = castGrid(4, 1000, Float64; h = 0.1, r0 = 0.000902, polynom=[0,0,1], msg=false);
    r1 = grid1.r;
    r2 = grid2.r;
    r3 = grid3.r;
    r4 = grid4.r;
    f1 = [gaussian(r1[n]) for n=1:grid1.N];
    f2 = [gaussian(r2[n]) for n=1:grid2.N];
    f3 = [gaussian(r3[n]) for n=1:grid3.N];
    f4 = [gaussian(r4[n]) for n=1:grid4.N];
#   -----------------------------------------------------------------------------------------
    f′1 = -r1 .* f1;
    f′2 = -r2 .* f2;
    f′3 = -r3 .* f3;
    f′4 = -r4 .* f4;
    o1 = grid_differentiation(f1, grid1);
    o2 = grid_differentiation(f2, grid2);
    o3 = grid_differentiation(f3, grid3);
    o4 = grid_differentiation(f4, grid4);

    println("r1max = ", r1[1000])
    println("r2max = ", r2[1000])
    println("r3max = ", r3[1000])
    println("r3max = ", r4[1000])

    println("f′1: ", f′1[992:1000])
    println("o1: ", o1[992:1000])
    println("-------")
    println(f′1[900])
    f′900 = grid_differentiation(f1, grid1, 900);
    println("-------")
    println(f′1[990])
    f′990 = grid_differentiation(f1, grid1, 990);
    println("-------")
    println(f′1[993])
    f′993 = grid_differentiation(f1, grid1, 993);

