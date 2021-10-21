function _show(io::IO, object)
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

_summary(io::IO, object) = print(io, nameof(typeof(object)))
