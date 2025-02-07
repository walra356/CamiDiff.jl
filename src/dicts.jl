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
