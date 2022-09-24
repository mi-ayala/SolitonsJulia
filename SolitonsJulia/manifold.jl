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

function F_DF!(F, DF, P, parameters, bundle)

    DF .= 0

    ### θ, σ,
    ### Prejection part
    project!(F, -f!(f, P, parameters) + (Derivative(1,0) + λ*Derivative(0,1))*P*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1])
    )
    
    #project!(component(DF, 2, 2), γ * Df!(component(DF, 2, 2), u, σ, ρ, β) - Derivative(1))


    ### Periodic orbit condition
    # γ = project(
    # Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
    # space(bundle))

    # component(component(F, 1),1) = component(component(P, 1),1) - component(γ, 1) 
    # component(component(F, 2),2) = component(component(P, 2),2) - component(γ, 2) 
    # component(component(F, 3),3) = component(component(P, 3),3) - component(γ, 3) 

    # ### Bundle condition  
    # component(component(F, 1),1) = component(component(P, 1),1) - component(bundle, 1) 
    # component(component(F, 2),2) = component(component(P, 2),2) - component(bundle, 2) 
    # component(component(F, 3),3) = component(component(P, 3),3) - component(bundle, 3) 

    return F, DF
end



N = 20
A = -3
E = 1





