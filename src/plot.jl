# NSDERungeKutta.jl
mutable struct _PhasePlot{object_T}
    object::object_T
end
@userplot PHASEPLOT
@recipe f(h::PHASEPLOT) = _PhasePlot(h.args[1])

# NSDETimeParallel.jl, NSDEMovingWindow.jl
mutable struct _Convergence{object_T}
    object::object_T
end
@userplot CONVERGENCE
@recipe f(h::CONVERGENCE) = _Convergence(h.args[1])
