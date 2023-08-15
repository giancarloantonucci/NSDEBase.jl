# NSDERungeKutta.jl
mutable struct _PhasePlot{plottable_T}
    plottable :: plottable_T
end
@userplot PHASEPLOT
@recipe f(h::PHASEPLOT) = _PhasePlot(h.args[1])

# NSDETimeParallel.jl, NSDEMovingWindow.jl
mutable struct _Convergence{plottable_T}
    plottable :: plottable_T
end
@userplot CONVERGENCE
@recipe f(h::CONVERGENCE) = _Convergence(h.args[1])
