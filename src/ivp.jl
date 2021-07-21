@doc """
    InitialValueProblem(rhs, u₀, tspan) <: NSDEProblem
    InitialValueProblem(rhs, u₀, t₀::Real, T::Real) <: NSDEProblem
    IVP(args...; kwargs...) <: NSDEProblem

returns a constructor for an initial value problem.

# Arguments
- `rhs   :: Union{Function, RightHandSideFunction}` : right-hand side derivative.
- `u₀    :: Union{Number, AbstractVector{Number}}`  : initial condition.
- `tspan :: Tuple{Real, Real}`                      : time domain.
"""
struct InitialValueProblem{rhs_T, u₀_T, tspan_T} <: NSDEProblem
    rhs::rhs_T
    u₀::u₀_T
    tspan::tspan_T
end

InitialValueProblem(rhs::RightHandSideFunction, u₀::Number, tspan) = InitialValueProblem(rhs, [u₀], tspan)
InitialValueProblem(f::Function, u₀, tspan) = InitialValueProblem(RHS(f), u₀, tspan)
InitialValueProblem(rhs::RightHandSideFunction, u₀, t₀::Real, T::Real) = InitialValueProblem(rhs, u₀, (t₀, T))
InitialValueProblem(f::Function, u₀, t₀, T) = InitialValueProblem(RHS(f), u₀, t₀, T)
@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

function Base.copy(problem::InitialValueProblem)
    @↓ rhs, u₀, tspan = problem
    return IVP(rhs, u₀, tspan)
end
