@doc """
    InitialValueProblem(rhs, u0, tspan) <: NSDEProblem
    InitialValueProblem(rhs, u0, t0::Real, tN::Real) <: NSDEProblem
    IVP(args...; kwargs...) <: NSDEProblem

returns a constructor for an initial value problem.

# Arguments
- `rhs   :: Union{Function, RightHandSideFunction}` : right-hand side derivative.
- `u0    :: Union{Number, AbstractVector{Number}}`  : initial condition.
- `tspan :: Tuple{Real, Real}`                      : time domain.
"""
struct InitialValueProblem{rhs_T, u0_T, tspan_T} <: NSDEProblem
    rhs::rhs_T
    u0::u0_T
    tspan::tspan_T
end

InitialValueProblem(rhs::RightHandSideFunction, u0::Number, tspan) = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(f::Function, u0, tspan) = InitialValueProblem(RHS(f), u0, tspan)
InitialValueProblem(rhs::RightHandSideFunction, u0, t0::Real, tN::Real) = InitialValueProblem(rhs, u0, (t0, tN))
InitialValueProblem(f::Function, u0, t0, tN) = InitialValueProblem(RHS(f), u0, t0, tN)
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

function Base.copy(problem::InitialValueProblem)
    @↓ rhs, u0, tspan = problem
    return IVP(rhs, u0, tspan)
end
