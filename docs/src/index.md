# NSDEBase.jl

This is the documentation of [NSDEBase.jl](https://github.com/giancarloantonucci/NSDEBase.jl), a Julia package containing shared types and functions of [NSDE.jl](https://github.com/giancarloantonucci/NSDE.jl) and its sub-packages.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API and thus documented in this manual.

### Abstract types

```@docs
AbstractNSDEObject
AbstractNSDEProblem
AbstractNSDESolver
AbstractNSDESolution
AbstractInitialValueProblem
AbstractInitialValueSolver
AbstractInitialValueSolution
AbstractRightHandSideFunction
```

### Composite types

```@docs
InitialValueProblem
subproblemof
RightHandSideFunction
SplitRightHandSideFunction
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
