module CamiDiff

import CamiMath

fwd = CamiMath.fwd
bwd = CamiMath.bwd
reg = CamiMath.reg
rev = CamiMath.rev

using Printf
using LinearAlgebra

export Grid
export gridname
export gridfunction
export castGrid
export findIndex
export grid_differentiation
export grid_integration

export fdiff_weight
export fdiff_expansion_weights
export fdiff_expansion
export fwd_diff_expansion_weights
export fdiff_interpolation_expansion_coeffs
export fdiff_interpolation_expansion_weights
export fdiff_interpolation
export fdiff_differentiation_expansion_coeffs
export fdiff_differentiation
export create_lagrange_differentiation_matrix
export fdiff_adams_moulton_expansion_coeff
export fdiff_adams_moulton_expansion_coeffs
export create_adams_moulton_weights
export fdiff_adams_bashford_expansion_coeff
export fdiff_adams_bashford_expansion_coeffs
export create_adams_bashford_weights
export trapezoidal_epw
export trapezoidal_integration


include("finite_differences.jl")
include("finite_difference_adams.jl")
include("grid.jl")


end
