@doc raw"""
    SplitRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t) = f_\text{s}(u, t) + f_\text{ns}(u, t)$.

# Constructors
```julia
SplitRightHandSide(f_s, f_ns)
SRHS(args...; kwargs...)
```

# Arguments
- `f_s :: Union{LinearRightHandSide, NonlinearRightHandSide}` : $f_\text{s}$, the stiff part of the right-hand side function $f$.
- `f_ns :: NonlinearRightHandSide}` : $f_\text{ns}$, the non-stiff part of the right-hand side function $f$.
"""
struct SplitRightHandSide{
    f_s_T <: Union{LinearRightHandSide, NonlinearRightHandSide},
    f_ns_T <: NonlinearRightHandSide,
    } <: AbstractRightHandSide
    f_s :: f_s_T
    f_ns :: f_ns_T
end

SplitRightHandSide(f_s::Function, f_ns::NonlinearRightHandSide) = SplitRightHandSide(NRHS(f_s), f_ns)
SplitRightHandSide(L::Union{Number, AbstractMatrix{<:Number}}, f_ns::NonlinearRightHandSide) = SplitRightHandSide(LRHS(L), f_ns)
SplitRightHandSide(f_s::Union{Number, AbstractMatrix{<:Number}, Function, LinearRightHandSide, NonlinearRightHandSide}, f_ns::Function) = SplitRightHandSide(f_s, NRHS(f_ns))
@doc (@doc SplitRightHandSide) SRHS(args...; kwargs...) = SplitRightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::SplitRightHandSide)(u, t)
    (rhs::SplitRightHandSide)(du, u, t)

returns the derivative `du` from the solution `u` and time `t`.
"""
function (rhs::SplitRightHandSide)(u, t)
    @↓ f_s, f_ns = rhs
    du = f_s(u, t) + f_ns(u, t)
    return du
end

function (rhs::SplitRightHandSide)(du, u, t)
    @↓ f_s, f_ns = rhs
    # @! du = f_s(u, t) + f_ns(u, t)
    f_s(du, u, t)
    v = similar(du)
    f_ns(v, u, t)
    du .+= v
    return du
end
