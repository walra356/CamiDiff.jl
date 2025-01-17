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

module CamiDiff

import CamiMath
fwd = CamiMath.fwd
bwd = CamiMath.bwd
reg = CamiMath.reg
rev = CamiMath.rev

#import LinearAlgebra

export Grid
export gridname
export gridfunction
export castGrid
export findIndex
export findΔn
export grid_interpolation
export grid_differentiation
export grid_integration

export fdiff_weight
export fdiff_expansion_weights
export fdiff_expansion
export fwd_diff_expansion_weights
export fdiff_interpolation_expansion_polynom
export fdiff_interpolation_expansion_weights
export fdiff_interpolation
export fdiff_differentiation_expansion_polynom
export create_lagrange_differentiation_matrix
export fdiff_adams_moulton_expansion_coeff
export fdiff_adams_moulton_expansion_polynom
export create_adams_moulton_weights
export fdiff_adams_bashford_expansion_coeff
export fdiff_adams_bashford_expansion_polynom
export create_adams_bashford_weights
export trapezoidal_epw
export trapezoidal_integration


include("finite_differences.jl")
include("finite_difference_adams.jl")
include("grid.jl")


end
