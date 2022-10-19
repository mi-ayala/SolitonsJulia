function get_soliton(N_F, N_T, parameters, λ, v, γ )
        
    ### Manifold space and first guess
    s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4

    S = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))
    
    for i in 1:4
        component(P, i)[(:,0)] .= component(γ, i)
        component(P, i)[(:,1)] .= component(v, i)
    end

    ### Tapes for Derivarive
    DF_ = LinearOperator(s_mani, s_mani, zeros(ComplexF64, 4*(2*N_F+1)*(N_T+1), 4*(2*N_F+1)*(N_T+1)))
    F_ = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

    newton!((F_, DF_, x) -> (F!(F_, x, parameters, N_F, N_T,v,λ, γ), DF!(DF_, x, parameters, N_F,λ)), P)          

    return P
end 

### Vector field 
function f(S, parameters, N_F, N_T)

    s_soliton = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4
    f=Sequence(s_mani, zeros(ComplexF64, dimension(s_soliton)))
    
    @unpack E, A, s = parameters   
    S₁, S₂, S₃, S₄ = eachcomponent(S)

    project!(component(f, 1), S₂ )
    project!(component(f, 2), - E*S₁ + A*S₃*S₁ + s*S₁*S₁*S₁ )
    project!(component(f, 3), + S₄)
    project!(component(f, 4), - 4*S₃)

    return f

end

### Jacobian
function Df!(Df, S, parameters)

    @unpack E, A, s = parameters
    S₁, S₂, S₃, S₄ = eachcomponent(S)

    project!(component(Df, 1, 1), 0*I) 
    project!(component(Df, 1, 2), Multiplication(one(S₂)))
    project!(component(Df, 1, 3), 0*I)
    project!(component(Df, 1, 4), 0*I)

    project!(component(Df, 2, 1), Multiplication(A*S₃ + s*S₁*S₁-E))
    project!(component(Df, 2, 2), 0*I)
    project!(component(Df, 2, 3), Multiplication(A*S₁))
    project!(component(Df, 2, 4), 0*I)

    project!(component(Df, 3, 1), 0*I)
    project!(component(Df, 3, 2), 0*I)
    project!(component(Df, 3, 3), 0*I)
    project!(component(Df, 3, 4), Multiplication(one(S₄)))

    project!(component(Df, 4, 1), 0*I)
    project!(component(Df, 4, 2), 0*I)
    project!(component(Df, 4, 3), Multiplication(-4*one(S₃)))
    project!(component(Df, 4, 4), 0*I)

    return Df
end

### Zero Finding Problem
function F!(F, S, para, N_F, N_T,v,λ,γ)

    project!(F, f(S, para, N_F, N_T) + Derivative(1,0)*S)

      for i = 1:4
          component(F, i) .= component(F, i) +  λ*component(Derivative(0,1)*P,i)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) 
          component(F, i)[(:,0)] .= component(P, i)[(:,0)] .- component(γ, i)
          component(F, i)[(:,1)] .= component(P, i)[(:,1)] .- component(v, i)
      end
 
    return F
end

### Derivative
function DF!(DF, P, para, N_F, λ)
  
        DF .= 0
        Df = copy(DF)

        Df = Df!(Df, P, para)

        M = project(Multiplication(Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) ), domain(component(Df, 1, 1)), codomain(component(Df,1,1)))
        
        project!(component(DF, 1, 1), -component(Df, 1, 1) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 1, 1)), codomain(component(Df,1,1)), ComplexF64))
        project!(component(DF, 1, 2), -component(Df, 1,2))
        project!(component(DF, 1, 3), -component(Df, 1,3))
        project!(component(DF, 1, 4), -component(Df, 1,4))
        project!(component(DF, 2, 1), -component(Df, 2,1))
        project!(component(DF, 2, 2), -component(Df, 2,2) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 2, 2)), codomain(component(Df,2,2)), ComplexF64))
        project!(component(DF, 2, 3), -component(Df, 2,3))
        project!(component(DF, 2, 4), -component(Df, 2,4))
        project!(component(DF, 3, 1), -component(Df, 3,1))
        project!(component(DF, 3, 2), -component(Df, 3,2))
        project!(component(DF, 3, 3), -component(Df, 3,3) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 3, 3)), codomain(component(Df,3,3)), ComplexF64))
        project!(component(DF, 3, 4), -component(Df, 3,4))
        project!(component(DF, 4, 1), -component(Df, 4,1))
        project!(component(DF, 4, 2), -component(Df, 4,2))
        project!(component(DF, 4, 3), -component(Df, 4,3))
        project!(component(DF, 4, 4), -component(Df, 4,4) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 4, 4)), codomain(component(Df,4,4)), ComplexF64))

         for i = 1:4
             component(DF,i,i)[(:,0),(:,0)] = I(2*N_F+1)
             component(DF,i,i)[(:,1),(:,1)] = I(2*N_F+1)
        end

        return DF
end


