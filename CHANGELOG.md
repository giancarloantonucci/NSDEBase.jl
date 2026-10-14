# Changelog

## 0.4.0 (unreleased)

### Breaking

- When a function has both an out-of-place `f(u, t)` and an in-place
  `f!(du, u, t)` method, the out-of-place form now wins by default. Pass
  `inplace=true` to keep the former choice.
- The multi-argument problem `copy` is removed; use `remake`.
- Supported Julia: `1.10` and later (was `1.6`).

### Added

- `remake`, `chunk_solver` and `withstepsize` extension interfaces.
- `SolverFailure`, the abstract type for numerical solver failures, and the
  documented interface a custom IVP must provide.
- An `inplace` keyword (`true`, `false` or `nothing`) on `RHS`, `LRHS`,
  `SRHS` and every function-based `IVP` constructor; split constructors also
  take a `(stiff, nonstiff)` pair.

### Changed

- Integer and rational initial data are converted to floating point, and
  mixed time endpoints are promoted, both at construction and in `remake`.
  Floating-point states are kept without copying.
- Complex linear operators promote real states, and function-based IVPs infer
  complex differentiation from complex initial data.
- Linear operators and states may have independent element types.
- The Wirtinger Jacobian is extracted directly, without a dense linear solve.
- Sources, tests and examples follow the repository JuliaFormatter
  configuration (Blue style), checked in CI.

### Fixed

- Forcing-function arities are validated, and optional arguments no longer
  cause a shorter signature to be misclassified.
- Derived in-place forcing that returns `nothing` works; standalone forcing
  uses dense buffers, and RHS calls use state-compatible buffers.

## 0.3.1

### Fixed

- `NonlinearRightHandSide(f; iscomplex=true)`: the Wirtinger Jacobian returned
  the `∂conj(f)/∂conj(u)` block instead of `∂f/∂u`.
- `NonlinearRightHandSide(f!)`: the derived out-of-place `f(u, t)` returns the
  filled buffer, not `f!`'s return value, so an `f!` returning `nothing` works
  with `rhs(u, t)` and the ForwardDiff Jacobian.
- `Lorenz96` rejects `N < 4` with an `ArgumentError`.

### Added

- Allocation-free four-argument RHS calls `rhs(du, v, u, t)` with caller-owned
  scratch `v`, on linear, nonlinear and split right-hand sides. Solver
  packages from NSDERungeKutta 0.2 use this form, which is why they require
  this release.
- `IVP(fₛ, fₙₛ, u0, tspan)` builds a split problem directly.
- Defaults for the three-positional ODE constructors (`Lorenz(u0, t0, tN)`
  etc.); `Lorenz96`'s default state places the perturbation first.
- Docstrings for the solver interface `solve`, `solve!`, `initialize_cache`,
  `initialize_solution`; API page; Aqua in the test suite.

### Changed

- Supported Julia: `1.6` and later (was `1.10`). Verified on 1.6–1.13.

### Removed

- The `Vector{T}(undef, N, d)` and `Vector{T}(undef, N2, N1, d)` methods
  (type piracy on `Base`).
