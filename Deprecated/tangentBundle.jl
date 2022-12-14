using .SolitonsJulia
using RadiiPolynomial


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

# # Complex{Float64} == ComplexF64

# N = 20
# A = -3
# E = 1

# s = ParameterSpace() × Fourier(N, 1.0)^4
# x = Sequence(s, rand(Complex{Float64}, dimension(s)))

# Df = jac_vectorfield(N, A, E)

# newton!((F, DF, x) -> (F!(F, x, Df), DF!(DF, x, Df)), x)

# λ = real(x[1])
# v = component(x, 2)
# newton!((F, DF, x) -> (F!(F, x, Df), DF!(DF, x, Df)), x)


# x̄ = x
# ### Bound Y
# x̄_interval = Interval.(x̄)
# F_interval = F!(x̄_interval)
# DF_interval = DF(x̄_interval, space(x̄_interval), space(x̄_interval), eltype(x̄_interval))