# NSDEBase.jl

A Julia package containing shared types and functions of [NSDE.jl](https://github.com/giancarloantonucci/NSDE.jl) and its sub-packages.

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDEBase.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDEBase.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDEBase.jl)

## Installation

NSDEBase is a [registered package](https://juliahub.com/ui/Search?q=NSDEBase&type=packages) compatible with Julia v1.6 and above. From the Julia REPL,

```
]add NSDEBase
```

## Usage

NSDEBase holds what every solver package shares: the abstract type tree, the right-hand-side objects, `InitialValueProblem`, the `solve`/`solve!` contract, and a set of well-known test problems.

```julia
using NSDEBase

# from a function — in-place or not; the missing form and the Jacobian are derived
problem = IVP((du, u, t) -> du .= u .* (1 .- u), [0.5], (0.0, 1.0))

# linear, f = L·u + g(t)
problem = IVP(reshape([-10.0], 1, 1), [0.1], (0.0, 1.0))

# split, f = L·u + fₙₛ(u, t), for IMEX and exponential solvers
problem = IVP(reshape([-10.0], 1, 1), (u, t) -> sin.(u), [0.1], (0.0, 1.0))

# cheap re-windowing for chunked and windowed solvers: the RHS is shared, not rebuilt
subproblem = copy(problem, [0.2], 0.5, 1.0)
```

Complex-valued problems get their Jacobian as a Wirtinger derivative (`RHS(f; iscomplex = true)`).

Bundled test problems: `Dahlquist`, `Logistic`, `SimplePendulum`, `DoublePendulum`, `VanDerPol`, `Rössler`, `Lorenz`, `Lorenz96`.

Read the [documentation](https://giancarloantonucci.github.io/NSDEBase.jl/dev) for a complete overview of this package.
