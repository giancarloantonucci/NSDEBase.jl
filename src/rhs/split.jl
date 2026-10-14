@doc raw"""
    SplitRightHandSide <: AbstractRightHandSide

A composite type for the right-hand side of an [`InitialValueProblem`](@ref) in the form $f(u, t) = f_\text{s}(u, t) + f_\text{ns}(u, t)$.

# Constructors
```julia
SplitRightHandSide(fₛ, fₙₛ)
SRHS(args...; kwargs...)
```

Raw functions accept `iscomplex` and `inplace`. The latter is a shared choice
or a `(stiff_choice, nonstiff_choice)` pair, each `true`, `false`, or `nothing`.
Prebuilt RHS objects retain their own calling convention.

# Arguments
- `fₛ::Union{LinearRightHandSide,NonlinearRightHandSide}` : $f_\text{s}$, the stiff part of the right-hand side function $f$
- `fₙₛ::NonlinearRightHandSide` : $f_\text{ns}$, the non-stiff part of the right-hand side function $f$
"""
struct SplitRightHandSide{
    fₛ_T<:Union{LinearRightHandSide,NonlinearRightHandSide},fₙₛ_T<:NonlinearRightHandSide
} <: AbstractRightHandSide
    fₛ::fₛ_T
    fₙₛ::fₙₛ_T
end

function split_inplace(inplace)
    choices = inplace isa Tuple && length(inplace) == 2 ? inplace : (inplace, inplace)
    all(x -> x === nothing || x isa Bool, choices) || throw(
        ArgumentError("inplace must be true, false, nothing, or a (stiff, nonstiff) pair of these"),
    )
    return choices
end

function SplitRightHandSide(
    fₛ::Function, fₙₛ::NonlinearRightHandSide; iscomplex=false, inplace=nothing
)
    a, _ = split_inplace(inplace)
    return SplitRightHandSide(RHS(fₛ; iscomplex, inplace=a), fₙₛ)
end
function SplitRightHandSide(
    L::Union{Number,AbstractMatrix{<:Number}},
    fₙₛ::NonlinearRightHandSide;
    iscomplex=false,
    inplace=nothing,
)
    split_inplace(inplace)
    return SplitRightHandSide(LRHS(L), fₙₛ)
end
function SplitRightHandSide(
    fₛ::Union{Number,AbstractMatrix{<:Number},Function,LinearRightHandSide,NonlinearRightHandSide},
    fₙₛ::Function;
    iscomplex=false,
    inplace=nothing,
)
    a, b = split_inplace(inplace)
    stiff = if fₛ isa Function
        RHS(fₛ; iscomplex, inplace=a)
    elseif fₛ isa Union{Number,AbstractMatrix}
        LRHS(fₛ)
    else
        fₛ
    end
    return SplitRightHandSide(stiff, RHS(fₙₛ; iscomplex, inplace=b))
end
@doc (@doc SplitRightHandSide) SRHS(args...; kwargs...) = SplitRightHandSide(args...; kwargs...)

#----------------------------------- METHODS -----------------------------------

"""
    (rhs::SplitRightHandSide)(u, t)
    (rhs::SplitRightHandSide)(du, u, t)
    (rhs::SplitRightHandSide)(du, v, u, t)

returns the derivative `du` from the solution `u` and time `t`. The 3-argument
in-place form allocates a temporary; the 4-argument form takes caller-owned
scratch `v` (same shape as `du`) instead and is allocation-free. `v` is used
first as scratch for `fₛ` and then to hold `fₙₛ(u, t)`, so a single buffer
covers the whole evaluation. Hot loops should use the 4-argument form.
"""
function (rhs::SplitRightHandSide)(u, t)
    @↓ fₛ, fₙₛ = rhs
    du = fₛ(u, t) + fₙₛ(u, t)
    return du
end

function (rhs::SplitRightHandSide)(du, u, t)
    @↓ fₛ, fₙₛ = rhs
    fₛ(du, u, t)
    v = similar(du)
    fₙₛ(v, u, t)
    du .+= v
    return du
end

function (rhs::SplitRightHandSide)(du, v, u, t)
    @↓ fₛ, fₙₛ = rhs
    fₛ(du, v, u, t) # fₛ is done with the scratch `v` once it returns
    fₙₛ(v, u, t)
    du .+= v
    return du
end
