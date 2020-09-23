"""
# Available lower level routines
* [`besseljs3d`](@ref) spherical Bessel function eval
"""
function lower_level_routs() end

"""
    function fj, fjder = besseljs3d(nterms,z;scale=1.0,ifder=0)
This subroutine evaluates the first `nterms` spherical Bessel 
functions, and if requested, their derivatives.
It incorporates a scaling parameter `scale` so that
      
      	fjs_n(z)=j_n(z)/SCALE^n
      	fjder_n(z)=\\frac{\\partial fjs_n(z)}{\\partial z}
# Input
* `nterms::Integer` order of expansion of output array `fjs` 
* `z::ComplexF64` argument of the spherical Bessel functions
* `scale::Float64` scaling factor
* `ifder::Integer1` flag indicating whether to calculate `fjder`
      	          0	NO
      	          1	YES
OUTPUT:
* `fjs::Array{ComplexF64}` array of length `nterms+1` of scaled Bessel functions.
* `fjder::Array{ComplexF64}` array of derivatives of scaled Bessel functions, if requested.
"""
function besseljs3d(nterms::Integer,z::ComplexF64;scale::Float64=1.0,
                    ifder::Integer=0)

    @assert (nterms >= 0)
    if (ifder !=0 && ifder != 1)
        @warn "unexpected value in ifder, no ders computed"
    end

    fjs = Array{ComplexF64}(undef,nterms+1)
    fjder = ComplexF64(0)

    
    
    if ifder == 1
        fjder = Array{ComplexF64}(undef,nterms+1)
    end
    
    # fortran interface
    # subroutine besseljs3d(nterms,z,scale,fjs,ifder,fjder)
    
    ccall((:besseljs3d_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,Fc),
          nterms,z,scale,fjs,ifder,fjder)

    if ifder != 1
        fjder = nothing
    end

    return fjs, fjder
    
end