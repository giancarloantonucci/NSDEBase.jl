# NSDEBase.jl

This is the documentation of [NSDEBase.jl](https://github.com/giancarloantonucci/NSDEBase.jl), a Julia package containing common types and utility functions shared by [NSDE.jl](https://github.com/giancarloantonucci/NSDE.jl) and its subpackages.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API, and thus documented in this manual.

### Abstract types

```@docs
AbstractNSDEType
AbstractNSDEProblem
AbstractNSDESolver
AbstractNSDESolution
AbstractInitialValueProblem
AbstractInitialValueSolver
AbstractInitialValueSolution
AbstractRightHandSide
```

### Composite types

```@docs
InitialValueProblem
subproblemof
RightHandSide
SplitRightHandSide
```

### ODEs

```@docs
Dahlquist
Logistic
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
