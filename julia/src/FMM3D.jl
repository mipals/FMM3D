"""
    module FMM3D
Wrappers for the Flatiron Institute's FMM3D
library. 
# wrappers
All N-body codes return output in an 
`FMMVals` structure. See documentation
of N-body codes for details.
N-body interactions with the Helmholtz kernel
- [`hfmm3d`](@ref): ``O(N)`` fast mutlipole code
- [`h3ddir`](@ref): ``O(N^2)`` direct code
N-body interactions with the Laplace kernel
- [`lfmm3d`](@ref): ``O(N)`` fast mutlipole code
- [`l3ddir`](@ref): ``O(N^2)`` direct code
# lower level routines
For a list of lower level routines see the 
documentation for [`lower_level_routs`](@ref)
"""
module FMM3D

using FMM3D_jll
import Base: show, propertynames

export FMMVals, hfmm3d, lfmm3d, h3ddir, l3ddir, lower_level_routs
export besseljs3d

# fortran input/return types

Fd = Ref{Float64}
Fi = Ref{Int32}
Fc = Ref{ComplexF64}

# common input types

TFN = Union{Array{Float64},Nothing}
TCN = Union{Array{ComplexF64},Nothing}


"""
    FMMVals
Abstract return type for FMM and direct computation 
function calls. Fields are `nothing` on return 
unless requested. See individual FMM/direct 
computation function documentation for specifics.
"""
abstract type FMMVals end

include("helmholtz.jl")
include("laplace.jl")
include("stokes.jl")
include("maxwell.jl")
include("lower_level_routines.jl")



# Probably not the prettiest print. Better than nothing
function propertynames(output::T) where {T<:Union{HelmholtzOutput,LaplaceOutput}}
    Base.fieldnames(typeof(output))
    pg  = output.pg
    pgt = output.pgt

    # Possible stages
    fieldstarg = [:pottarg,:gradtarg,:hesstarg]
    fields = [:pot,:grad,:hess]

    # Only using the computed fields
    fieldstarg = fieldstarg[1:pgt]

    if pg == nothing
        return tuple(fieldstarg...,)
    else
        fields = fields[1:pg]
    end

    return tuple(cat(fields,fieldstarg,dims=1)...)
    
end
function show(io::IO, ::MIME"text/plain", output::T) where {T<:FMMVals}
    println("$T. Accesible fields:")
    for field in propertynames(output)
        println("."*string(field))
    end
end


end # module