# Defining problems

NSDEBase describes the equation; a solver package supplies the integration method.
The four main layers are:

| Package | Role |
|---|---|
| NSDEBase | Problems, right-hand sides, and shared interfaces |
| NSDERungeKutta | Serial explicit, implicit, IMEX, and exponential integration |
| NSDETimeParallel | Parareal on a grid of time chunks |
| NSDEMovingWindow | A sequence of overlapping Parareal windows |

ArrowMacros supplies convenience macros used by the implementation. You do not
need those macros to write a problem or inspect a solution.

## Start with a function

Run this page from top to bottom, or run `examples/guided_tour.jl` from the
repository root in an environment containing NSDEBase.

```@example problems
using NSDEBase
pendulum = IVP((u, t) -> [u[2], -sin(u[1])], [π/4, 0.0], (0.0, 10.0))
function lorenz!(du, u, t)
    du[1] = 10 * (u[2] - u[1])
    du[2] = u[1] * (28 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8/3) * u[3]
    return nothing
end
lorenz = IVP(lorenz!, [2.0, 3.0, -14.0], (0.0, 20.0))
@assert lorenz.rhs(lorenz.u0, 0.0) ≈ [10.0, 81.0, 130/3]
nothing # hide
```

An out-of-place `f(u, t)` returns a state-shaped vector. An in-place
`f!(du, u, t)` fills `du`; its return value is ignored. Avoid shared scratch
storage in a closure if the problem will be solved on multiple threads.

The shorter signature wins when both forms are present, including functions
with optional arguments. Resolve ambiguity explicitly with `inplace=true` or
`inplace=false`. An explicit choice is checked for compatible arity, not for all
possible state types. For a split problem, `inplace=(true, false)` chooses the
stiff and nonstiff conventions separately.

## Preserve the structure of the equation

```@example problems
# u' = L*u: implicit solvers can use a direct linear solve.
L = [-1.0 0.0; 0.0 -100.0]
decay = IVP(L, [1.0, 1.0], (0, 1))

# u' = L*u + g(t): forcing takes time alone.
forced = IVP(LRHS(L, t -> [sin(t), 0.0]), [0.0, 0.0], (0, 1))

# u' = L*u + f_nonstiff(u, t): for IMEX or exponential methods.
split = IVP(L, (u, t) -> [sin(u[2]), cos(u[1])], [1.0, 0.0], (0, 1))
@assert decay.tspan === (0.0, 1.0)
@assert split.rhs isa SplitRightHandSide
nothing # hide
```

`IVP(L, f, u0, tspan)` is a **split** equation: `f` takes the state and time.
Use `LRHS(L, g)` when the additive forcing takes time alone. Prebuilt `RHS`,
`LRHS`, and `SRHS` objects retain their calling conventions.

## Types, Jacobians, and remaking a problem

```@example problems
integer_data = IVP((u, t) -> -u, [1, 2], (0, 1))
rotation = IVP(-1im, [1.0], (0, 1))
new_start = [1.0, 1.0, 1.0]
short = remake(lorenz; u0=new_start, tspan=(0.0, 2.0))
@assert integer_data.u0 isa Vector{Float64}
@assert eltype(rotation.u0) <: Complex
@assert short.rhs === lorenz.rhs
@assert short.u0 === new_start
nothing # hide
```

Integer/rational data become floating-point data. Already-floating state buffers
are retained when no promotion is needed: `remake` is not a deep copy. Use
`copy(new_start)` if you need separate ownership. A complex linear operator
promotes a real state; for complex forcing with a real operator, supply a complex
initial state yourself.

Function-based IVPs infer complex differentiation from the state. Jacobian-based
complex solves require a holomorphic right-hand side; represent non-holomorphic
dynamics as real and imaginary equations. Keeping a Float32 or BigFloat time grid
does not automatically increase the precision of a solver's tableau coefficients.

Built-in problems include `Dahlquist`, `Logistic`, `SimplePendulum`,
`DoublePendulum`, `VanDerPol`, `Rössler`, `Lorenz`, and `Lorenz96`.
For example, `Logistic(0.3, (0.0, 2.0))` specifies a start and a span explicitly.

## Solve and reuse

All solver layers support `solve(problem, solver)`. Reuse is deliberately shown
with the concrete interfaces implemented by each package:

| Solver | Cached call |
|---|---|
| Runge–Kutta | `solve!(cache, solution, problem, solver)` |
| Parareal | `parareal(cache, solution, problem)` |
| MoWi | `mowi(cache, solution, problem)` |

Runge–Kutta and Parareal allocate through `initialize_cache` and
`initialize_solution`. MoWi currently uses `NSDEMovingWindow.MoWiCache` and
`MoWiSolution` instead; see its guide for the complete allocation example.

Reuse avoids rebuilding some storage; it is not a zero-allocation promise.
Follow each solver's restrictions on state shape, type, time span, and operator.

## Extending the interfaces

A custom problem for Parareal or MoWi needs writable `rhs`, `u0`, and `tspan`
fields and `remake(problem; u0=problem.u0, tspan=problem.tspan)`. The remade problem
must retain the supplied floating-point state buffer. A serial solver implements
`solve`, `initialize_cache`, `initialize_solution`, and `solve!`; its returned solution
must provide the time/state storage and interpolation needed by the caller.
See NSDETimeParallel's extension guide for a complete wrapper.

Use a subtype of `SolverFailure` for a numerical failure to produce a usable
step. Do not classify a programming error as a numerical failure. The serial
solver throws; a Parareal backend can turn a numerical failure into a diagnosed
`:nonfinite` return so a window strategy can retry.
