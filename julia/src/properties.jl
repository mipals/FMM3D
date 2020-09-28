function propertynames(output::T) where {T<:Union{HelmholtzOutput,LaplaceOutput}}
    fieldnames(typeof(output))
    pg  = output.pg
    pgt = output.pgt

    fieldstarg = [:pottarg,:gradtarg,:hesstarg]
    fieldstarg = fieldstarg[1:pgt]

    if pg === nothing
        return tuple(fieldstarg...,)
    end

    fields = [:pot,:grad,:hess]
    fields = fields[1:pg]

    return tuple(cat(fields,fieldstarg,dims=1)...)
    
end


function propertynames(output::StokesOutput)
    fieldnames(typeof(output))
    pg  = output.pg
    pgt = output.pgt

    fieldstarg = [:pottarg,pretarg,:gradtarg]
    fieldstarg = fieldstarg[1:pgt]

    if pg === nothing
        return tuple(fieldstarg...,)
    end

    fields = [:pot,:pre,:grad]
    fields = fields[1:pg]

    return tuple(cat(fields,fieldstarg,dims=1)...)
    
end


function propertynames(output::MaxwellOutput)
    fieldnames(typeof(output))
    ifeh  = output.ifeh
    ifehtarg = output.ifehtarg

    fieldstarg = [:etarg,:gradetarg,:htarg,:gradhtarg]
    if ifehtarg <= 2
        fieldstarg = fieldstarg[1:ifehtarg]
    elseif ifehtarg <= 4
        fieldstarg = fieldstarg[3:ifehtarg]
    elseif ifehtarg == 5 
        fieldstarg = [:etarg,:htarg]
    end

    fields = [:e,:grade,:h,:gradh]
    if ifeh === nothing
        return tuple(fieldstarg...,)
    elseif ifeh <= 2
        fields = fields[1:ifeh]
    elseif ifeh <= 4
        fields = fields[3:ifeh]
    elseif ifeh == 5
        fields = [:e,:h]
    end

    return tuple(cat(fields,fieldstarg,dims=1)...)
    
end

# Probably not the prettiest print. Better than nothing
function show(io::IO, ::MIME"text/plain", output::T) where {T<:FMMVals}
    for field in propertynames(output)
        strfield = string(field)
        println("."*strfield*" "^(max(9-length(strfield),0)), 
                '\t', "$(typeof(getfield(output,field)))")
    end
end