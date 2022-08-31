
"""
f_bundle(x)

Tangent bundle vector field
```
"""
function f_bundle!(f, v, p)

    @unpack E, A, s = p   
    v₁, v₂ = eachcomponent(v)
    
    ### Cos2t Fourier Series
    cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

    project!(component(f, 1),  - v₂ )
    project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
    return f
end


function Df_bundle!(Df, v, p)

    @unpack E, A, s = p

    v₁, v₂, v₃ = eachcomponent(v)
    project!(component(Df, 1, 1), Multiplication(-σ*one(v₁)))
    project!(component(Df, 1, 2), Multiplication(σ*one(v₂)))
    project!(component(Df, 2, 1), Multiplication(ρ-v₃))
    project!(component(Df, 2, 2), Multiplication(-one(v₂)))
    return Df
end


function F_DF!(F, DF, x, σ, ρ, β)
    γ, u = x[1], component(x, 2)
    DF .= 0

    F[1] =
        (sum(component(u, 1)) - 10.205222700615433) * 24.600655549587863 +
        (sum(component(u, 2)) - 11.899530531689562) * (-2.4927169722923335) +
        (sum(component(u, 3)) - 27.000586375896557) * 71.81142025024573
        component(component(DF, 1, 2), 1)[1,:] .= 24.600655549587863
        component(component(DF, 1, 2), 2)[1,:] .= -2.4927169722923335
        component(component(DF, 1, 2), 3)[1,:] .= 71.81142025024573

    project!(component(F, 2), γ * f!(component(F, 2), u, σ, ρ, β) - differentiate(u))
    f!(component(DF, 2, 1), u, σ, ρ, β)
    project!(component(DF, 2, 2), γ * Df!(component(DF, 2, 2), u, σ, ρ, β) - Derivative(1))

    return F, DF
end

n = 400
x̄ = Sequence(ParameterSpace() × Fourier(n, 1.0)^3, zeros(ComplexF64, 1+3*(2n+1)))
newton!((F, DF, x) -> F_DF!(F, DF, x, σ, ρ, β), x̄)