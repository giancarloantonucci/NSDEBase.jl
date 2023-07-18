"""
    show(io::IO, object::AbstractObject)

prints the full description of an `object` and its contents to a stream `io`.
"""
function Base.show(io::IO, object::AbstractObject; depth=4)
    summary(io, object)
    names = propertynames(object)
    num_properties = length(names)
    nothing_properties = 0
    for n = num_properties:-1:1
        if getproperty(object, names[n]) isa Nothing
            nothing_properties += 1
        else
            break
        end
    end
    num_properties -= nothing_properties
    padding = get(io, :padding, "")
    current_level = get(io, :current_level, 1)
    for (n, name) in enumerate(names)
        property = getproperty(object, name)
        if property isa Nothing
            continue
        end
        is_property_last = n == num_properties
        level = current_level
        print(io, "\n", padding, (is_property_last ? "└─ " : "├─ "), string(name))
        if current_level < depth
            level += 1
            io2 = IOContext(io, :padding => string(padding, is_property_last ? "   " : "│  "), :current_level => level)
            print(io, " :: ")
            summary(io, property)
            if property isa Number || property isa AbstractArray
                print(io, " = ")
                show(IOContext(io2, :limit => true), property)
            else
                show(io2, property)
            end
        else
            print(io, " :: ")
            summary(io, property)
        end
    end
end

"""
    summary(io::IO, object::AbstractObject)

prints the short description of an `object` to a stream `io`.
"""
Base.summary(io::IO, object::AbstractObject) = print(io, nameof(typeof(object)))
