_summary(io::IO, object) = print(io, nameof(typeof(object)))

function _show(io::IO, object)
    print(io, nameof(typeof(object)), ":")
    pad = get(io, :pad, "")
    names = propertynames(object)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(object, name)
        if field !== nothing
            print(io, "\n", pad, "   ‣ " * string(name) * " ≔ ")
            show(IOContext(io, :pad => string(pad, "   ")), field)
        end
    end
end
