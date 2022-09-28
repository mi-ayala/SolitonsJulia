using .SolitonsJulia
using RadiiPolynomial, UnPack

function f!(f, P, parameters)

    @unpack E, A, s = parameters   
    P₁, P₂, P₃, P₄ = eachcomponent(P)

    project!(component(f, 1), P₂ )
    project!(component(f, 2), - E*P₁ + A*P₃*P₁ + s*P₁*P₁*P₁ )
    project!(component(f, 3), + P₄)
    project!(component(f, 4), - 4*P₃)

    return f

end

function Df!(Df, P, parameters)

    @unpack E, A, s = parameters
    P₁, P₂, P₃, P₄ = eachcomponent(P)

    project!(component(Df, 1, 1), Multiplication(zero(P₁)))
    project!(component(Df, 1, 2), Multiplication(one(P₂)))
    project!(component(Df, 1, 3), Multiplication(zero(P₃)))
    project!(component(Df, 1, 4), Multiplication(zero(P₄)))

    project!(component(Df, 2, 1), Multiplication(A*P₃ + s*P₁*P₁))
    project!(component(Df, 2, 2), Multiplication(zero(P₂)))
    project!(component(Df, 2, 3), Multiplication(A*P₁))
    project!(component(Df, 2, 4), Multiplication(zero(P₄)))

    project!(component(Df, 3, 1), Multiplication(zero(P₁)))
    project!(component(Df, 3, 2), Multiplication(zero(P₂)))
    project!(component(Df, 3, 3), Multiplication(zero(P₃)))
    project!(component(Df, 3, 4), Multiplication(one(P₄)))

    project!(component(Df, 4, 1), Multiplication(zero(P₁)))
    project!(component(Df, 4, 2), Multiplication(zero(P₂)))
    project!(component(Df, 4, 3), Multiplication(-4*one(P₃)))
    project!(component(Df, 4, 4), Multiplication(zero(P₄)))

    return Df
end

# function F_DF!(F, DF, P, f, para, bundle,N_T,N_F)

#     ### Function
#     project!(F, -f!(f, P, para) + Derivative(1,0)*P)
    
#     γ = project(
#         Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#         space(bundle))
    
#     for i = 1:4
#     component(F, i) .= component(F, i) +  λ*component(Derivative(0,1)*P,i)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) 
#     component(F, i)[(:,0)] .= component(P, i)[(:,0)] .- component(γ, i)
#     component(F, i)[(:,1)] .= component(P, i)[(:,1)] .- component(bundle, i)
#     end
        
#     ### Derivative
#     DF .= 0

#     project!(component(DF, 1, 1), component(Df, 1, 1) + Derivative(1,0) + λ*Derivative(0,1))
#     project!(component(DF, 1, 2), component(Df, 1,2))
#     project!(component(DF, 1, 3), component(Df, 1,3))
#     project!(component(DF, 1, 4), component(Df, 1,4))
#     project!(component(DF, 2, 1), component(Df, 2,1))
#     project!(component(DF, 2, 2), component(Df, 2,2) + Derivative(1,0) + λ*Derivative(0,1))
#     project!(component(DF, 2, 3), component(Df, 2,3))
#     project!(component(DF, 2, 4), component(Df, 2,4))
#     project!(component(DF, 3, 1), component(Df, 3,1))
#     project!(component(DF, 3, 2), component(Df, 3,2))
#     project!(component(DF, 3, 3), component(Df, 3,3) + Derivative(1,0) + λ*Derivative(0,1))
#     project!(component(DF, 3, 4), component(Df, 3,4))
#     project!(component(DF, 4, 1), component(Df, 4,1))
#     project!(component(DF, 4, 2), component(Df, 4,2))
#     project!(component(DF, 4, 3), component(Df, 4,3))
#     project!(component(DF, 4, 4), component(Df, 4,4) + Derivative(1,0) + λ*Derivative(0,1))

#     for i = 1:4
#         component(Df_,i,i)[(i,0),(i,0)] = I(N_T)
#     end

#     return F, DF
# end


# N = 20
# A = -3
# E = 1

