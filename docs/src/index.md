# NSDEBase.jl

Problem definitions and shared interfaces for the NSDE ecosystem.

## First example

```julia
using NSDEBase
problem = IVP((u,t) -> u .* (1 .- u), [0.3], (0.0, 2.0))
short = remake(problem; u0=[0.4], tspan=(0.0, 1.0))
@assert short.rhs === problem.rhs
```

Start with the [guided tour](guide.md), then consult:

- [API](api.md)

## Environment

These pages describe the current development checkout (Julia 1.10 or later).
Use compatible versions of the NSDE packages in one project. For local
checkouts, develop each package with `Pkg.develop(path="/path/to/package")`;
`Pkg.status()` reports what will actually load. A registered release may expose
an older API. Plotting examples additionally need Plots.
