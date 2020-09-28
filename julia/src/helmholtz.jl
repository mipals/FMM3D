"""
    HelmholtzOutput
Return type for Helmholtz computations.
Fields are `nothing` on return unless requested.
See individual FMM/direct computation function 
documentation for specifics.
"""
mutable struct HelmholtzOutput <: FMMVals
    pot
    grad

    pottarg
    gradtarg

    pgt
    pg
end
function HelmholtzOutput(pgt,pg=nothing)
    HelmholtzOutput(nothing,nothing,
                    nothing,nothing,
                    pgt,pg)
end

"""
```julia
    vals = hfmm3d(eps,zk,sources;charges=nothing,dipvecs=nothing,
                  targets=nothing,pg=0,pgt=0,nd=1)
```
This function computes the N-body Helmholtz interactions
in three dimensions where the interaction kernel is given 
by ``e^{ikr}/r`` and its gradients. This is the 
``O(N)`` fast multipole code which computes the interactions
to the requested precision.
```math
 u(x) = \\sum_{j=1}^{N} c_{j} \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
- v_{j} \\cdot \\nabla \\left( \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
\\right)  \\, , 
```
 
where ``c_{j}`` are the charge densities,  
      ``v_{j}`` are the dipole orientation vectors, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum
# Input
* `eps::Float64` precision requested
* `zk::ComplexF64` Helmholtz parameter 
* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)
* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `pg::Integer` source eval flag. 
    + Potential (``u``) at sources evaluated if `pg == 1`. 
    + Potential and gradient (``\\nabla u``) at sources evaluated if `pg == 2`
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities
# Output
        
`vals::HelmholtzOutput` with the fields
* `vals.pot::Array{ComplexF64}` size (nd,n) or (n) potential at source locations if requested
* `vals.grad::Array{ComplexF64}` size (nd,3,n) or (3,n) gradient at source locations if requested
* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
"""
function hfmm3d(eps::Float64,zk::Union{Float64,ComplexF64},
                sources::Array{Float64};
                charges::TCN=nothing,dipvecs::TCN=nothing,
                targets::TFN=nothing,pg::Integer=0,pgt::Integer=0,
                nd::Integer=1)

    # check inputs
    
    @assert size(sources,1) == 3
    @assert nd >= 0
    @assert (0 <= pg  && pg  <= 2) "flag pg not in expected range"
    @assert (0 <= pgt && pgt <= 2) "flag pgt not in expected range"   

    # default values

    vals = HelmholtzOutput(pgt,pg)
    
    ifcharge = 0
    ifdipole = 0

    zero = ComplexF64(0)
    
    pot = zero
    grad = zero
    hess = zero
    pottarg = zero
    gradtarg = zero
    hesstarg = zero

    n = div(length(sources),3)
    nt = 0
    zk = complex(zk)

    if targets !== nothing
        @assert size(targets,1) == 3
        nt = div(length(targets),3)
    else
        targets = 0.0
    end
    
    if charges !== nothing
        @assert div(length(charges),nd) == n
        ifcharge = 1
    else
        charges = zero
        ifcharge = 0
    end

    if dipvecs !== nothing
        @assert div(length(dipvecs),nd) == n*3
        ifdipole = 1
    else
        dipvecs = zero
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        return vals
    end

    if (pg != 1 && pg != 2) && (pgt != 1 && pgt != 2)
        @warn "no output requested, doing nothing"
        return vals
    end

    if (pgt == 1 || pgt == 2) && targets == nothing
        @warn "target values requested but no targets provided"
    end

    # allocate memory for return values
    
    if pg == 1 || pg == 2
        if nd > 1
            pot = zeros(ComplexF64,nd,n)
        else
            pot = zeros(ComplexF64,n)
        end
    end

    if pg == 2
        if nd > 1
            grad = zeros(ComplexF64,nd,3,n)
        else
            grad = zeros(ComplexF64,3,n)
        end
    end

    if pgt == 1 || pgt == 2
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,3,nt)
        else
            gradtarg = zeros(ComplexF64,3,nt)
        end
    end

    # actually call the function
    #
    # fortran calling sequence:
    # subroutine hfmm3d(nd,eps,zk,nsource,source,ifcharge,
    #     $    charge,ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,
    #     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg)
    #

    
    ccall((:hfmm3d_,libfmm3d),Cvoid,(Fi,Fd,Fc,Fi,Fd,Fi,Fc,
                                    Fi,Fc,Fi,Fc,Fc,Fc,Fi,
                                    Fd,Fi,Fc,Fc,Fc),
          nd,eps,zk,n,sources,ifcharge,charges,ifdipole,
          dipvecs,pg,pot,grad,hess,nt,targets,pgt,
          pottarg,gradtarg,hesstarg)

    # load requested values

    if pg == 1 || pg == 2; vals.pot = pot end
    if pg == 2; vals.grad = grad end
    if pgt == 1 || pgt == 2; vals.pottarg = pottarg end
    if pgt == 2; vals.gradtarg = gradtarg end
    
    return vals

end

"""
```julia
    vals = h3ddir(zk,sources,targets;charges=nothing,
                    dipvecs=nothing,pgt=0,nd=1,
                    thresh=1e-16)
```
This function computes the N-body Helmholtz interactions
in three dimensions where the interaction kernel is given 
by ``e^{ikr}/r`` and its gradients. This is the 
``O(N^2)`` direct evaluation code.
```math
 u(x) = \\sum_{j=1}^{N} c_{j} \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
- v_{j} \\cdot \\nabla \\left( \\frac{e^{ik \\|x-x_{j}\\|}}{\\|x-x_{j}\\|} 
\\right)  \\, , 
```
 
where ``c_{j}`` are the charge densities,  
      ``v_{j}`` are the dipole orientation vectors, and 
      ``x_{j}`` are the source locations.
      When ``x=x_{m}``, the term corresponding to 
      ``x_{m}`` is dropped from the sum
# Input
* `zk::ComplexF64` Helmholtz parameter 
* `sources::Array{Float64}` size (3,n) source locations (``x_{j}``)
* `targets::Array{Float64}` size (3,nt) target locations (``x``)
* `charges::Array{ComplexF64}` size (nd,n) or (n) charge densities (c_{j})
* `dipvecs::Array{ComplexF64}` size (nd,3,n) or (3,n)) dipole orientation vectors (``v_{j}``)
* `pgt::Integer` target eval flag. 
    + Potential at targets evaluated if `pgt == 1`. 
    + Potential and gradient at targets evaluated if `pgt == 2`
* `nd::Integer` number of densities
* `thresh::Float64` threshold for ignoring interactions when ``\\|x-x_{j}\\| \\leq thresh``
# Output
        
`vals::HelmholtzOutput` with the fields
* `vals.pottarg::Array{ComplexF64}` size (nd,nt) or (nt) potential at target locations if requested
* `vals.gradtarg::Array{ComplexF64}` size (nd,3,nt) or (3,nt) gradient at target locations if requested
"""
function h3ddir(zk::Union{ComplexF64,Float64},sources::Array{Float64},
                targets::Array{Float64};
                charges::TCN=nothing,dipvecs::TCN=nothing,
                pgt::Integer=0,nd::Integer=1,
                thresh::Float64=1e-16)

    # check inputs
    
    @assert size(sources,1) == 3
    @assert size(targets,1) == 3    
    @assert nd >= 0
    @assert (0 <= pgt || pgt <= 2) "flag pgt not in expected range"
    
    # default values

    vals = HelmholtzOutput(pgt)
    
    ifcharge = 0
    ifdipole = 0

    zero = ComplexF64(0)
    
    pottarg = zero
    gradtarg = zero
    hesstarg = zero  

    n = div(length(sources),3)
    nt = div(length(targets),3)
    zk = complex(zk)
    
    if charges !== nothing
        @assert div(length(charges),nd) == n
        ifcharge = 1
    else
        charges = zero
        ifcharge = 0
    end

    if dipvecs !== nothing
        @assert div(length(dipvecs),nd) == n*3
        ifdipole = 1
    else
        dipvecs = zero
        ifdipole = 0
    end

    if ifcharge == 0 && ifdipole == 0
        @warn "no charges or dipoles provided, doing nothing"
        return vals
    end

    if (pgt != 1 && pgt != 2)
        @warn "no output requested, doing nothing"
        return vals
    end

    # allocate memory for return values
    
    if pgt == 1 || pgt == 2
        if nd > 1
            pottarg = zeros(ComplexF64,nd,nt)
        else
            pottarg = zeros(ComplexF64,nt)
        end
    end

    if pgt == 2
        if nd > 1
            gradtarg = zeros(ComplexF64,nd,3,nt)
        else
            gradtarg = zeros(ComplexF64,3,nt)
        end
    end

    # dispatch to appropriate wrapper
    
    if pgt == 1
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdp!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,thresh)
                vals.pottarg = pottarg
            else
                h3ddirectcp!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             thresh)
                vals.pottarg = pottarg
            end
        else
            if ifdipole == 1
                h3ddirectdp!(nd,zk,sources,
                             dipvecs,n,targets,nt,
                             pottarg,thresh)
                vals.pottarg = pottarg
            end
        end
    elseif pgt == 2
        if ifcharge == 1
            if ifdipole == 1
                h3ddirectcdg!(nd,zk,sources,charges,
                              dipvecs,n,targets,nt,
                              pottarg,gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                
            else
                h3ddirectcg!(nd,zk,sources,charges,
                             n,targets,nt,pottarg,
                             gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                                
            end
        else
            if ifdipole == 1
                h3ddirectdg!(nd,zk,sources,
                              dipvecs,n,targets,nt,
                             pottarg,gradtarg,thresh)
                vals.pottarg = pottarg
                vals.gradtarg = gradtarg                                
            end
        end
    end                
            
    return vals

end

function h3ddirectcp!(nd,zk,sources,charges,n,targets,nt,
                      pot,thresh)
    ccall((:h3ddirectcp_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fd),
          nd,zk,sources,charges,n,targets,nt,pot,thresh)
    return
end

function h3ddirectdp!(nd,zk,sources,dipvecs,n,targets,nt,
                      pot,thresh)
    ccall((:h3ddirectdp_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fd),
          nd,zk,sources,dipvecs,n,targets,nt,pot,thresh)
    return
end

function h3ddirectcdp!(nd,zk,sources,charges,dipvecs,n,targets,
                       nt,pot,thresh)
    ccall((:h3ddirectcdp_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fc,Fi,
                                          Fd,Fi,Fc,Fd),
          nd,zk,sources,charges,dipvecs,n,targets,nt,pot,thresh)
    return
end

function h3ddirectcg!(nd,zk,sources,charges,n,targets,nt,
                      pot,grad,thresh)
    ccall((:h3ddirectcg_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fd),
          nd,zk,sources,charges,n,targets,nt,pot,grad,thresh)
    return
end

function h3ddirectdg!(nd,zk,sources,dipvecs,n,targets,nt,
                      pot,grad,thresh)
    ccall((:h3ddirectdg_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fd),
          nd,zk,sources,dipvecs,n,targets,nt,pot,grad,thresh)
    return
end

function h3ddirectcdg!(nd,zk,sources,charges,dipvecs,n,targets,
                       nt,pot,grad,thresh)
    ccall((:h3ddirectcdg_,libfmm3d),Cvoid,(Fi,Fc,Fd,Fc,Fc,Fi,
                                          Fd,Fi,Fc,Fc,Fd),
          nd,zk,sources,charges,dipvecs,n,targets,nt,pot,grad,thresh)
    return
end
