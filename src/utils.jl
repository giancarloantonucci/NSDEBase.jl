# ToDo: NSDEVector
Vector{u_T}(undef, n, d) where u_T = Vector{u_T}[Vector{u_T}(undef, d) for i = 1:n]
Vector{u_T}(undef, n₂, n₁, d) where u_T = Vector{Vector{u_T}}[Vector{u_T}(undef, n₁, d) for i = 1:n₂]

# ToDo: zero!(v::NSDEVector)
zero!(v::AbstractVector{<:Number}) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector{<:Number}})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end
