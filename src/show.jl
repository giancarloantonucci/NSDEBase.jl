"""
    show(io::IO, object::AbstractNSDEObject)

prints the full description of an `object` and its contents to a stream `io`.
"""
function Base.show(io::IO, object::AbstractNSDEObject)
    print(io, nameof(typeof(object)), ":")
    padding = get(io, :padding, "")
    thislevel = get(io, :thislevel, 1)
    names = propertynames(object)
    N = length(names)
    M = 0
    for n = N:-1:1
        field = getproperty(object, names[n])
        field === nothing ? M += 1 : break
    end
    N -= M
    for (n, name) in enumerate(names)
        level = thislevel
        field = getproperty(object, name)
        if field !== nothing
            print(io, "\n", padding, n < N ? "├── " : "└── ", string(name))
            if thislevel < 3 && !(field isa AbstractArray)
                level += 1
                tmpio = IOContext(io, :padding => string(padding, n < N ? "│   " : "    "), :thislevel => level)
                print(tmpio, " = ")
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

prints the short description of an `object` to a stream `io`.
"""
Base.summary(io::IO, object::AbstractNSDEObject) = print(io, nameof(typeof(object)))
