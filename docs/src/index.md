```@meta
CurrentModule = CamiDiff
```

# Home

`CamiDiff.jl` is a [Julia](http://julialang.org) package for one-dimensional finite-difference analysis. 

---

### Install

The package is installed using the Julia package manager

```
julia> using Pkg; Pkg.add("CamiDiff")

julia> using CamiDiff
```

### Introduction

`CamiDiff` has been developed to analyze mathematical functions, provided by the user 
*in tabulated form*. The package is based on finite-difference analysis methods in one dimension. 

Throughout the documentation ``f(x)`` will be the 'function of interest' under investigation. In tabulated 
form we write ``f[n]``. This form can be regarded as the result of the discretization of ``f(x)`` 
on a [`Grid`](@ref) of ``N`` points, addressable by the *grid index* ``n = 1, ⋯ N``. The [`Grid`](@ref) 
can be linear (uniform) or non-linear as specified by a [`gridfunction`](@ref) - see [Discretization](@ref).

The current implementation of `CamiDiff` was developed for *real functions of a single variable*. 
A set of four predefined [`gridfunction`](@ref)`s` is included in the package: 'exponential', 
'truncated-exponential', 'linear' and 'polynomial'. These are restricted to the domain ``[0, ∞)``. 
When applicable, we shall underline this restriction to the non-negative domain by using the variable ``r`` 
rather than ``x``, writing ``f(r)`` rather than ``f(x)``, with the implicit condition ``r ≥ 0``.