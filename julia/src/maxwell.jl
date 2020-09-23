"""
    MaxwellOutput
Return type for Maxwell computations.
Fields are `nothing` on return unless requested.
See individual FMM/direct computation function 
documentation for specifics.
"""
mutable struct MaxwellOutput <: FMMVals
    e
    h
    grade
    gradh

    etarg
    htarg
    gradetarg
    gradhtarg
end
function MaxwellOutput()
    MaxwellOutput(nothing,nothing,nothing,nothing,
                  nothing,nothing,nothing,nothing)
end