# define BlockVectors, i.e. Vectors of Vectors
abstract type AbstractBlockVector{T} <: AbstractVector{T} end
struct BlockVector{T} <: AbstractBlockVector{T} end
function BlockVector{T}(::UndefInitializer, n_chunks, l_chunks) where T
    return Vector{T}[Vector{T}(undef, l_chunks) for i = 1:n_chunks]
end

# fill Vectors and BlockVectors with zeros
zero!(v::AbstractVector) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector})
    for i = 1:length(v)
        zero!(v[i])
    end
    return v
end

# from BlockVectors to Vectors
function flatview!(v_flat, v)
    n_chunks = length(v)
    l_chunks = length(v[1])
    for i = 1:n_chunks
        for j = 1:l_chunks
            v_flat[(i-1)*l_chunks+j] = v[i][j]
        end
    end
    return v_flat
end

# from Vectors to BlockVectors
function nestedview!(v, v_flat)
    n_chunks = length(v)
    l_chunks = length(v[1])
    for i = 1:n_chunks
        for j = 1:l_chunks
            v[i][j] = v_flat[(i-1)*l_chunks+j]
        end
    end
    return v
end
