# SPDX-License-Identifier: MIT

# Copyright (c) 2023 Jook Walraven <69215586+walra356@users.noreply.github.com> and contributors

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

# ==============================================================================
#                               dicts.jl
# ==============================================================================


# ============================ dictBernoulliNumbers ==================================

dictBernoulliNumbers = Dict{Int,Rational{Int}}(
  
0 => 1 // 1, 1 => -1 // 2, 2 => 1 // 6, 3 => 0 // 1,
  4 => -1 // 30, 5 => 0 // 1, 6 => 1 // 42, 7 => 0 // 1,
  8 => -1 // 30, 9 => 0 // 1, 10 => 5 // 66, 11 => 0 // 1,
  12 => -691 // 2730, 13 => 0 // 1, 14 => 7 // 6, 15 => 0 // 1,
  16 => -3617 // 510, 17 => 0 // 1, 18 => 43867 // 798, 19 => 0 // 1,
  20 => -174611 // 330, 21 => 0 // 1, 22 => 854513 // 138, 23 => 0 // 1,
  24 => -236364091 // 2730, 25 => 0 // 1, 26 => 8553103 // 6, 27 => 0 // 1,
  28 => -23749461029 // 870, 29 => 0 // 1, 30 => 8615841276005 // 14322, 31 => 0 // 1,
  32 => -7709321041217 // 510, 33 => 0 // 1, 34 => 2577687858367 // 6, 35 => 0 // 1
  )

# ============================ dictBigConversion ======================================

@doc raw"""
    dictBigConversion

Dictionary for conversion from Int-based types to BigInt-based types
#### Example:
```
julia> dictBigConversion
Dict{DataType, DataType} with 12 entries:
  Vector{Rational{Int64}}         => Vector{Rational{BigInt}}
  Int64                           => BigInt
  Vector{Float64}                 => Vector{BigFloat}
  ComplexF64                      => Complex{BigFloat}
  Vector{Vector{ComplexF64}}      => Vector{Vector{Complex{BigFloat}}}
  Vector{Int64}                   => Vector{BigInt}
  Vector{Vector{Int64}}           => Vector{Vector{BigInt}}
  Vector{Vector{Float64}}         => Vector{Vector{BigFloat}}
  Vector{Vector{Rational{Int64}}} => Vector{Vector{Rational{Int64}}}
  Float64                         => BigFloat
  Rational{Int64}                 => Rational{BigInt}
  Vector{ComplexF64}              => Vector{Complex{BigFloat}}
```
"""
dictBigConversion = Dict(
  
  Int => BigInt,
  Vector{Int} => Vector{BigInt},
  Vector{Vector{Int}} => Vector{Vector{BigInt}},
  Rational{Int} => Rational{BigInt},
  Vector{Rational{Int}} => Vector{Rational{BigInt}},
  Vector{Vector{Rational{Int}}} => Vector{Vector{Rational{BigInt}}},
  Float64 => BigFloat,
  Vector{Float64} => Vector{BigFloat},
  Vector{Vector{Float64}} => Vector{Vector{BigFloat}},
  Complex{Float64} => Complex{BigFloat},
  Vector{Complex{Float64}} => Vector{Complex{BigFloat}},
  Vector{Vector{Complex{Float64}}} => Vector{Vector{Complex{BigFloat}}}

  )
