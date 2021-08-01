"""
    InitialValueProblem{rhs_T, u0_T, tspan_T} <: NSDEProblem

defines a constructor for an initial value problem.

$(TYPEDFIELDS)
"""
struct InitialValueProblem{rhs_T, u0_T, tspan_T} <: NSDEProblem
    "right-hand side derivative"
    rhs::rhs_T
    "initial condition"
    u0::u0_T
    "time domain"
    tspan::tspan_T
end

# """
#     InitialValueProblem(rhs, u0, tspan) :: InitialValueProblem
#     IVP(args...; kwargs...) :: InitialValueProblem
#
# returns an `InitialValueProblem` from:
# - `rhs   :: Union{Function, RightHandSideFunction}` : right-hand side derivative.
# - `u0    :: Union{Number, AbstractVector{Number}}`  : initial condition.
# - `tspan :: Tuple{Real, Real}`                      : time domain.
# """
InitialValueProblem(rhs::RightHandSideFunction, u0::Number, tspan) = InitialValueProblem(rhs, [u0], tspan)
InitialValueProblem(f::Function, u0, tspan) = InitialValueProblem(RHS(f), u0, tspan)

"""
    InitialValueProblem(rhs, u0, t0::Real, tN::Real) :: InitialValueProblem
    IVP(args...; kwargs...) :: InitialValueProblem

returns an `InitialValueProblem` with `tspan = (t0, tN)`.
"""
InitialValueProblem(rhs::RightHandSideFunction, u0, t0::Real, tN::Real) = InitialValueProblem(rhs, u0, (t0, tN))
InitialValueProblem(f::Function, u0, t0, tN) = InitialValueProblem(RHS(f), u0, t0, tN)

@doc (@doc InitialValueProblem) IVP(args...; kwargs...) = InitialValueProblem(args...; kwargs...)

"""
    copy(problem::InitialValueProblem) :: InitialValueProblem

returns a copy of `problem`.
"""
function Base.copy(problem::InitialValueProblem)
    @↓ rhs, u0, tspan = problem
    return IVP(rhs, u0, tspan)
end
