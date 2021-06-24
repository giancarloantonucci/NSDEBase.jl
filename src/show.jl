function Base.show(io::IO, rhs::RightHandSideFunction)
    if get(io, :compact, false)
        print(io, "RightHandSideFunction")
    else
        @↓ f, f!, Df, Df! = rhs
        print(io,
            "RightHandSideFunction:\n",
            "  f: ", typeof(f), "\n",
            "  f!: ", typeof(f!), "\n",
            "  Df: ", typeof(Df), "\n",
            "  Df!: ", typeof(Df!), "\n",
        )
    end
end

function Base.show(io::IO, ivp::InitialValueProblem)
    if get(io, :compact, false)
        print(io, "InitialValueProblem")
    else
        @↓ rhs, u0, tspan = ivp
        print(io,
            "InitialValueProblem:\n",
            "  rhs: ",  typeof(rhs), "\n",
            "  u0: ", u0, "\n",
            "  tspan: ", tspan, "\n",
        )
    end
end
