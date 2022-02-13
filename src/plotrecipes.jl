mutable struct _PhasePlot{object_T}
    object::object_T
end
@userplot PHASEPLOT
@recipe function f(h::PHASEPLOT)
    return _PhasePlot(h.args[1])
end

mutable struct _Convergence{object_T}
    object::object_T
end
@userplot CONVERGENCE
@recipe function f(h::CONVERGENCE)
    return _Convergence(h.args[1])
end
