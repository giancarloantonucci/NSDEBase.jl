# NSDEBase.jl

This is the documentation of [NSDEBase.jl](https://github.com/giancarloantonucci/NSDEBase.jl), the root package of the NSDE ecosystem: the shared abstract types, the right-hand-side objects, `InitialValueProblem`, the `solve`/`solve!` contract that every solver package extends, and a set of well-known test problems.

## Installation

From the Julia REPL,

```
]add https://github.com/giancarloantonucci/NSDEBase.jl
```

## Right-hand sides and problems

An `InitialValueProblem` (alias `IVP`) pairs a right-hand side with an initial value and a time span. The right-hand side comes in three kinds, and the solver packages dispatch on them:

- [`NonlinearRightHandSide`](@ref) (`RHS`) — generic ``f(u, t)``. Build it from a function, in-place or not: the missing form and the Jacobian (ForwardDiff, or a Wirtinger derivative for complex problems via `iscomplex = true`) are derived.
- [`LinearRightHandSide`](@ref) (`LRHS`) — ``f = Lu + g(t)``, enabling direct linear solves in implicit methods.
- [`SplitRightHandSide`](@ref) (`SRHS`) — ``f = f_s + f_{ns}``, for IMEX (and, on the roadmap, exponential) solvers.

```julia
using NSDEBase
problem = IVP((du, u, t) -> du .= u .* (1 .- u), [0.5], (0.0, 1.0))
problem = IVP(reshape([-10.0], 1, 1), [0.1], (0.0, 1.0))                  # linear
problem = IVP(reshape([-10.0], 1, 1), (u, t) -> sin.(u), [0.1], (0.0, 1.0)) # split
```

`copy(problem, u0, tspan)` re-windows a problem cheaply — the right-hand side is shared, not rebuilt — which is the seam the chunked and windowed solvers stand on. The in-place right-hand-side calls have a 4-argument form taking caller-owned scratch, so hot loops stay allocation-free.

Bundled test problems: `Dahlquist`, `Logistic`, `SimplePendulum`, `DoublePendulum`, `VanDerPol`, `Rössler`, `Lorenz`, `Lorenz96`.

See the [API](api.md) for the full reference.
