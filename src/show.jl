function lastnothingsof(object)
    names = propertynames(object)
    N = length(names)
    M = 0
    for n = N:-1:1
        field = getproperty(object, names[n])
        if field === nothing
            M += 1
        else
            break
        end
    end
    return N - M
end

"""
    show(io::IO, object::AbstractNSDEObject)

prints a description of an [`AbstractNSDEObject`](@ref) and all its contents to a stream `io`.
"""
function Base.show(io::IO, object::AbstractNSDEObject)
    print(io, nameof(typeof(object)))
    padding = get(io, :padding, "")
    thislevel = get(io, :thislevel, 1)
    names = propertynames(object)
    N = lastnothingsof(object)
    for (n, name) in enumerate(names)
        level = thislevel
        field = getproperty(object, name)
        if field !== nothing
            string₁ = (n < N ? "├── " : "└── ")
            print(io, "\n", padding, string₁, string(name))
            if thislevel < 3
                level += 1
                string₂ = (n < N ? "│   " : "    ")
                tmpio = IOContext(io, :padding => string(padding, string₂), :thislevel => level)
                string₃ = (field isa AbstractNSDEObject ? " :: " : " = ")
                print(tmpio, string₃)
                show(tmpio, field)
            else
                print(io, " :: ")
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
