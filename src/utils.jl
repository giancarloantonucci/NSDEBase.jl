# struct NSDEVector{data_T, idxs_T} <: AbstractVector
#     data::data_T
#     idxs::idxs_T
# end
#
# v = NSDEVector{u_T}(undef, n, d)
# v[i] == @view v.data[(i-1)*d+1:i*d]
# v[mask(i)] = mask(v, i) = v.data[i]
# extract(v, i) = [v.data[(j-1)*d+i] for j = 1:n]
#
# v = NSDEVector{u_T}(undef, n₁, n₂, d)
# v = NSDEVector{u_T}(undef, n₁, n₂, n₃, d)

# Vector{u_T}(undef, n, d) where u_T = Vector{u_T}[Vector{u_T}(undef, d) for i = 1:n]
# Vector{u_T}(undef, n₂, n₁, d) where u_T = Vector{Vector{u_T}}[Vector{u_T}(undef, n₁, d) for i = 1:n₂]

zero!(v::AbstractVector{<:Number}) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector{<:Number}})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end
