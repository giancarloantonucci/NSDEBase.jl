"""
    show(io::IO, object::AbstractObject)

prints the full description of an `object` and its contents to a stream `io`.
"""
function Base.show(io::IO, object::AbstractObject)
    print(io, nameof(typeof(object)))
    padding = get(io, :padding, "")
    thislevel = get(io, :thislevel, 1)
    namesoffields = propertynames(object)
    numberofnames = length(namesoffields)
    numberoffieldstoignore = 0
    for n = numberofnames:-1:1
        field = getproperty(object, namesoffields[n])
        field isa Nothing ? numberoffieldstoignore += 1 : break
    end
    numberofnames -= numberoffieldstoignore
    for (n, nameoffield) in enumerate(namesoffields)
        level = thislevel
        field = getproperty(object, nameoffield)
        if field !== nothing
            print(io, "\n", padding, n < numberofnames ? "├─── " : "└─── ", string(nameoffield))
            if thislevel < 5 && !(field isa AbstractArray)
                level += 1
                tmpio = IOContext(io, :padding => string(padding, n < numberofnames ? "│    " : "     "), :thislevel => level)
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
    summary(io::IO, object::AbstractObject)

prints the short description of an `object` to a stream `io`.
"""
Base.summary(io::IO, object::AbstractObject) = print(io, nameof(typeof(object)))
