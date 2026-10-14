# Changelog

## 0.4.0 development revision (unreleased)

- Add SolverFailure and document the mutable custom-IVP interface. Extract the Wirtinger block without a dense solve. Validate forced RHS arities and forward inplace through all function IVP forms and split constructors. When both signatures exist, out-of-place now wins by default; use inplace=true to retain the former preference.
- Promote real states for complex linear operators and infer complex differentiation in function-based IVPs. Add explicit RHS `inplace` overrides and handle optional arguments without misclassifying their shorter signatures.
- Normalize integer/rational initial data and mixed time endpoints at problem construction, preserving floating-state aliases. Allow independently typed linear operators and states.
- Fix derived in-place forcing that returns nothing; use dense standalone forcing buffers and state-compatible buffers in RHS calls.
- Promote mixed time-span endpoints in `remake`.
- Added `remake`, `chunk_solver`, and `withstepsize` extension interfaces. Removed multi-argument problem `copy`; use `remake` directly. Raised the Julia minimum to 1.10.
- All Julia sources, tests, and examples follow the repository JuliaFormatter configuration. Historical implementation notes have been removed from current API descriptions.

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
