# NSDEBase.jl

This is the documentation of [NSDEBase.jl](https://github.com/giancarloantonucci/NSDEBase.jl), a Julia package containing shared types and functions of [NSDE.jl](https://github.com/giancarloantonucci/NSDE.jl) and its sub-packages.

## Manual

```@contents
Depth = 3
```

## API

All exported types and functions are considered part of the public API and thus documented in this manual.

### Abstract types

```@autodocs
Modules = [NSDEBase]
Pages = ["abstract.jl"]
```

### Composite types

```@autodocs
Modules = [NSDEBase]
Pages = ["ivp.jl", "rhs.jl", "srhs.jl"]
```

### ODEs

```@autodocs
Modules = [NSDEBase]
Pages = ["odes.jl"]
```

### Utilities

```@docs
show
summary
```

## Index

```@index
```
