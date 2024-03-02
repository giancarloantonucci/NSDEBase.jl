# Used by NSDERungeKutta.jl
mutable struct _PhasePlot{plottable_T}
    plottable::plottable_T
end
@userplot PHASEPLOT
RecipesBase.recipetype(::Val{:phaseplot}, args...) = PHASEPLOT(args)
@recipe f(h::PHASEPLOT) = _PhasePlot(h.args[1])

# Used by NSDETimeParallel.jl and NSDEMovingWindow.jl
mutable struct _Convergence{plottable_T}
    plottable::plottable_T
end
@userplot CONVERGENCE
RecipesBase.recipetype(::Val{:convergence}, args...) = CONVERGENCE(args)
@recipe f(h::CONVERGENCE) = _Convergence(h.args[1])
