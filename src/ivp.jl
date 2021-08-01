"""
    InitialValueProblem{rhs_T, u0_T, tspan_T} <: NSDEProblem

returns a constructor for an initial value problem.

---

    InitialValueProblem(rhs, u0, tspan)
    IVP(args...; kwargs...)

returns an `InitialValueProblem` with:
- `rhs   :: Union{Function, RightHandSideFunction}` : right-hand side derivative.
- `u0    :: Union{Number, AbstractVector{Number}}`  : initial condition.
- `tspan :: Tuple{Real, Real}`                      : time domain.

---

    InitialValueProblem(rhs, u0, t0::Real, tN::Real)
    IVP(args...; kwargs...)

returns an `InitialValueProblem` with `tspan = (t0, tN)`.
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

"""
    copy(problem::InitialValueProblem)

returns a copy of `problem`.
"""
function Base.copy(problem::InitialValueProblem)
    @↓ rhs, u0, tspan = problem
    return IVP(rhs, u0, tspan)
end

Base.summary(io::IO, problem::InitialValueProblem) = print(io, "InitialValueProblem")

function Base.show(io::IO, problem::InitialValueProblem)
    print(io, "InitialValueProblem:\n")
    pad = get(io, :pad, "")
    newline = get(io, :newline, "\n")
    names = propertynames(problem)
    N = length(names)
    for (n, name) in enumerate(names)
        field = getproperty(problem, name)
        print(io, pad, "   ‣ " * string(name) * " := ")
        show(IOContext(io, :pad => "   ", :newline => ""), field)
        n == N ? print(io, newline) : print(io, "\n")
    end
end
