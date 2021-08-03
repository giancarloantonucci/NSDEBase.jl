# NSDEBase.jl

This is the documentation of [NSDEBase.jl](https://github.com/antonuccig/NSDEBase.jl), a Julia package containing shared types and functions of [NSDE.jl](https://github.com/antonuccig/NSDE.jl) and its sub-packages.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API and thus documented in this manual.

### Types

```@docs
NSDEProblem
NSDESolver
NSDESolution
RightHandSideFunction
InitialValueProblem
InitialValueSolver
InitialValueSolution
```

### ODEs

```@docs
Dahlquist
Logistic
Riccati
SimplePendulum
VanderPol
Rössler
Lorenz
Lorenz96
```

### Utilities

```@docs
show
summary
```

## Index

```@index
```
