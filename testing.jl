
using RadiiPolynomial, UnPack

# function jac_vectorfield(N, A, E)
#     Df_ = LinearOperator(Fourier(N,1.0)^4, Fourier(N,1.0)^4, zeros(ComplexF64, 4*(2*N+1), 4*(2*N+1)))
#     cos2t = Sequence(Fourier(2, 1.0), [.5, 0, 0, 0, .5])

#     project!(component(Df_, 1, 2), I)
#     project!(component(Df_, 2, 1), Multiplication(A*cos2t - E))
#     project!(component(Df_, 3, 4), I)
#     project!(component(Df_, 4, 3), -4*I)
#     return Df_
# end

# function F!(F, x, Df)

#     λ = component(x,1)[1]
#     v =  component(x,2)

#     F[1] = component(v, 1)(0)[0] - 1
#     project!(component(F, 2), (Derivative(1) - Df + λ*I)*v  )

#     return F
# end

# function DF!(DF, x, Df)
#     DF .= 0

#     λ = component(x,1)[1]
#     v = component(x,2)


#     project!(component(component(DF,1,2), 1), Evaluation(0))
#     project!(component(DF,2,1), v)
#     project!(component(DF,2,2), Derivative(1) - Df + λ*I)
   
   
#     return DF
# end

# N = 20
# A = -3
# E = 1

# s = ParameterSpace() × Fourier(N, 1.0)^4
# x = Sequence(s, rand(Complex{Float64}, dimension(s)))

# Df = jac_vectorfield(N, A, E)
# newton!((F, DF, x) -> (F!(F, x, Df), DF!(DF, x, Df)), x)

# λ = real(x[1])
# v = component(x, 2)

N_T = 1
N_F = 1
s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4
DF = LinearOperator(s_mani, s_mani, zeros(ComplexF64, 4*(2*N_F+1)*(N_T+1), 4*(2*N_F+1)*(N_T+1)))
P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))
F = P
f = F
bundle = P

para = soliton_parameters()

    ### Function
    project!(F, -f!(f, P, para) + Derivative(1,0)*P)
    
    γ = project(
        Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
        space(bundle))
    
    for i = 1:4
    component(F, i) .= component(F, i) +  λ*component(Derivative(0,1)*P,i)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) 
    component(F, i)[(:,0)] .= component(P, i)[(:,0)] .- component(γ, i)
    component(F, i)[(:,1)] .= component(P, i)[(:,1)] .- component(bundle, i)
    end
        

    F
    # ### Derivative
    # DF .= 0

    # project!(component(DF, 1, 1), component(Df, 1, 1) + Derivative(1,0) + λ*Derivative(0,1))
    # project!(component(DF, 1, 2), component(Df, 1,2))
    # project!(component(DF, 1, 3), component(Df, 1,3))
    # project!(component(DF, 1, 4), component(Df, 1,4))
    # project!(component(DF, 2, 1), component(Df, 2,1))
    # project!(component(DF, 2, 2), component(Df, 2,2) + Derivative(1,0) + λ*Derivative(0,1))
    # project!(component(DF, 2, 3), component(Df, 2,3))
    # project!(component(DF, 2, 4), component(Df, 2,4))
    # project!(component(DF, 3, 1), component(Df, 3,1))
    # project!(component(DF, 3, 2), component(Df, 3,2))
    # project!(component(DF, 3, 3), component(Df, 3,3) + Derivative(1,0) + λ*Derivative(0,1))
    # project!(component(DF, 3, 4), component(Df, 3,4))
    # project!(component(DF, 4, 1), component(Df, 4,1))
    # project!(component(DF, 4, 2), component(Df, 4,2))
    # project!(component(DF, 4, 3), component(Df, 4,3))
    # project!(component(DF, 4, 4), component(Df, 4,4) + Derivative(1,0) + λ*Derivative(0,1))

    # for i = 1:4
    #     component(Df_,i,i)[(i,0),(i,0)] = I(N_T)
    # end


