"""
    show(io::IO, object::AbstractNSDEObject)

prints a description of an [`AbstractNSDEObject`](@ref) and all its contents to a stream `io`.
"""
function Base.show(io::IO, object::AbstractNSDEObject)
    print(io, nameof(typeof(object)))
    # print(io, nameof(typeof(object)), ":")
    padding = get(io, :padding, "")
    this_lvl = get(io, :this_lvl, 1)
    names = propertynames(object)
    N = length(names)
    for (n, name) in enumerate(names)
        next_lvl = this_lvl
        field = getproperty(object, name)
        if field !== nothing
            # print(io, "\n", padding, "   ‣ " * string(name) * " ≔ ")
            print(io, "\n", padding, (n < N ? "├── " : "└── ") * string(name) * " ≔ ")
            if this_lvl < 2
                next_lvl += 1
                iocontext = IOContext(
                    io,
                    # :padding => string(padding, "   "),
                    :padding => string(padding, "│   "),
                    :this_lvl => next_lvl,
                )
                show(iocontext, field)
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
