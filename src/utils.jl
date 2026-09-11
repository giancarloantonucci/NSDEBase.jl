# NSDEBase/src/utils.jl

"""
    zero!(v::AbstractVector)

fills `v` with zeros in-place. Works on flat vectors of numbers and on nested
vectors of vectors (e.g. stage arrays), recursing one level.
"""
zero!(v::AbstractVector{<:Number}) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector{<:Number}})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end
