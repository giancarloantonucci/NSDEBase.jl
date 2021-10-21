"""
    show(io::IO, object::AbstractNSDEObject)

prints a full description of an [`AbstractNSDEObject`](@ref) and all its contents to a stream `io`.
"""
function Base.show(io::IO, object::AbstractNSDEObject)
    print(io, nameof(typeof(object)), ":")
    pad = get(io, :pad, "")
    lvl = get(io, :lvl, 1)
    names = propertynames(object)
    N = length(names)
    for (n, name) in enumerate(names)
        lvl_ = lvl
        field = getproperty(object, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            if lvl < 2
                lvl_ += 1
                show(IOContext(io, :pad => string(pad, "   "), :lvl => lvl_), field)
            else
                summary(io, field)
            end
        end
    end
end

"""
    summary(io::IO, problem::AbstractNSDEObject)

prints a brief description of an [`AbstractNSDEObject`](@ref) and all its contents to a stream `io`.
"""
Base.summary(io::IO, object::AbstractNSDEObject) = print(io, nameof(typeof(object)))
