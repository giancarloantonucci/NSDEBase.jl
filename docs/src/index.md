# NSDEBase.jl

This is the documentation of [NSDEBase.jl](https://github.com/antonuccig/NSDEBase.jl), a Julia package containing shared types and functions of [NSDE.jl](https://github.com/antonuccig/NSDE.jl) and its sub-packages.

## Manual

```@contents
Depth = 3
```

## API

Only exported types and functions are considered a part of the public API of NSDEBase. All such objects should be documented in this manual.

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

## Index

```@index
```
