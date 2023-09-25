# 📌 Write NSDEVector
Vector{u_T}(undef, N, d) where u_T = Vector{u_T}[Vector{u_T}(undef, d) for i = 1:N]
Vector{u_T}(undef, N2, N1, d) where u_T = Vector{Vector{u_T}}[Vector{u_T}(undef, N1, d) for i = 1:N2]

# 📌 Change to zero!(v::NSDEVector)
zero!(v::AbstractVector{<:Number}) = fill!(v, zero(eltype(v)))
function zero!(v::AbstractVector{<:AbstractVector{<:Number}})
    for i in eachindex(v)
        zero!(v[i])
    end
    return v
end
