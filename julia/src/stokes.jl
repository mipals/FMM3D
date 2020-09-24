"""
    StokesOutput
Return type for Stokes computations.
Fields are `nothing` on return unless requested.
See individual FMM/direct computation function 
documentation for specifics.
"""
mutable struct StokesOutput <: FMMVals
    pot
    pre
    grad

    pottarg
    pretarg
    gradtarg
end
function StokesOutput()
    StokesOutput(nothing,nothing,nothing,
                 nothing,nothing,nothing)
end
