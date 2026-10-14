# NSDEBase.jl

Problem definitions and shared interfaces for the NSDE ecosystem.

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://giancarloantonucci.github.io/NSDEBase.jl/dev) ![Build Status](https://img.shields.io/github/actions/workflow/status/giancarloantonucci/NSDEBase.jl/CI.yml) ![Coverage Status](https://img.shields.io/codecov/c/github/giancarloantonucci/NSDEBase.jl)

## Start here

```julia
using NSDEBase
problem = IVP((u,t) -> u .* (1 .- u), [0.3], (0.0, 2.0))
short = remake(problem; u0=[0.4], tspan=(0.0, 1.0))
@assert short.rhs === problem.rhs
```

It supplies `IVP`, `RHS`, `LRHS`, `SRHS`, built-in problems, and the shared
solver interfaces. It does not integrate a problem by itself.

Integer data are converted to floating point. Complex linear operators promote
real states; function-based problems infer complex differentiation from their
state. Jacobian-based complex solves require a holomorphic right-hand side.
`remake` shares an already-compatible state buffer: copy it when independent
ownership is needed.

## Learn by running

The [guided tour](docs/src/guide.md) explains each step. Its companion
[`examples/guided_tour.jl`](examples/guided_tour.jl) runs the numerical examples
with assertions. Run it in a project containing this package and its example
dependencies; use Julia 1.10 or later.

The ecosystem is NSDEBase → NSDERungeKutta → NSDETimeParallel → NSDEMovingWindow.
ArrowMacros provides supporting macros; they are not required in user examples.

## Installation and versions

For a registered release, use `Pkg.add("NSDEBase")`. These docs describe the
current development API; a registered version may differ. To use a local
checkout, activate your working environment and `Pkg.develop(path="/path/to/NSDEBase")`.
Develop matching local dependency checkouts in that environment too, then check
`Pkg.status()`. The serial and parallel guides also require NSDERungeKutta;
plotting requires Plots.

See the [documentation](https://giancarloantonucci.github.io/NSDEBase.jl/dev/)
and [API reference](docs/src/api.md) for details.

## Build the manual locally

From the package root, with matching dependency checkouts already available:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs --threads=4 docs/make.jl
```

Documenter executes the guided-tour examples during the build. For local
unreleased dependencies, develop their paths into `docs` before instantiating.
The deployment step only publishes in an appropriately configured CI environment.
