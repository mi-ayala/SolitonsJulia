using RadiiPolynomial

function f!(f, P, E, A, s)
    P₁, P₂, P₃, P₄ = eachcomponent(P)
    project!(component(f, 1), P₂ )
    project!(component(f, 2), - E*P₁ + A*P₃*P₁ - s*P₁*P₁*P₁ )
    project!(component(f, 3), u₁*u₂ + P₄)
    project!(component(f, 4), u₁*u₂ - 4*P₃)
    return f
end

function Df!(Df, u, σ, ρ, β)
    P₁, P₂, P₃, P₄ = eachcomponent(P)
    project!(component(Df, 1, 1), Multiplication(-σ*one(u₁)))
    project!(component(Df, 1, 2), Multiplication(σ*one(u₂)))
    project!(component(Df, 1, 3), Multiplication(zero(u₃)))
    project!(component(Df, 1, 4), Multiplication(zero(u₃)))
    project!(component(Df, 2, 1), Multiplication(ρ-u₃))
    project!(component(Df, 2, 2), Multiplication(-one(u₂)))
    project!(component(Df, 2, 3), Multiplication(-u₁))
    project!(component(Df, 2, 4), Multiplication(zero(u₃)))
    project!(component(Df, 3, 1), Multiplication(u₂))
    project!(component(Df, 3, 2), Multiplication(u₁))
    project!(component(Df, 3, 3), Multiplication(-β*one(u₃)))
    project!(component(Df, 3, 4), Multiplication(zero(u₃)))
    project!(component(Df, 4, 1), Multiplication(u₂))
    project!(component(Df, 4, 2), Multiplication(u₁))
    project!(component(Df, 4, 3), Multiplication(-β*one(u₃)))
    project!(component(Df, 4, 4), Multiplication(zero(u₃)))
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




γ = project(
    Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
    space(v))


    (Derivative(1) - Df + λ*α*I)

function jac_vectorfield(N, A, E)
    Df_ = LinearOperator(Fourier(N,1.0)^4, Fourier(N,1.0)^4, zeros(ComplexF64, 4*(2*N+1), 4*(2*N+1)))
    cos2t = Sequence(Fourier(2, 1.0), [.5, 0, 0, 0, .5])

    project!(component(Df_, 1, 2), I)
    project!(component(Df_, 2, 1), Multiplication(A*cos2t - E))
    project!(component(Df_, 3, 4), I)
    project!(component(Df_, 4, 3), -4*I)
    return Df_
end


function F!(F, x, Df)

    λ = component(x,1)[1]
    v =  component(x,2)

    F[1] = component(v, 1)(0)[0] - 1
    project!(component(F, 2), (Derivative(1) - Df + λ*I)*v  )

    return F
end

function DF!(DF, x, Df)
    DF .= 0

    λ = component(x,1)[1]
    v = component(x,2)


    project!(component(component(DF,1,2), 1), Evaluation(0))
    project!(component(DF,2,1), v)
    project!(component(DF,2,2), Derivative(1) - Df + λ*I)
   
   
    return DF
end
# Complex{Float64} == ComplexF64

N = 20
A = -3
E = 1

s = ParameterSpace() × Fourier(N, 1.0)^4
x = Sequence(s, rand(Complex{Float64}, dimension(s)))

Df = jac_vectorfield(N, A, E)

newton!((F, DF, x) -> (F!(F, x, Df), DF!(DF, x, Df)), x)

λ = real(x[1])
v = component(x, 2)
newton!((F, DF, x) -> (F!(F, x, Df), DF!(DF, x, Df)), x)
