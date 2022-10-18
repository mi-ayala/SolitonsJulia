using .SolitonsJulia
using RadiiPolynomial, UnPack

### Parameters
N_F = 3
A = -3
E = 1
para = soliton_parameters()

### Bundle
λ, v = get_bundle(N_F, A, E)

### Manifold initial guess
N_T = 2
s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4

### Periodic orbit data
P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

#  γ = project(
#      Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#      space(v))

 γ = project(
     Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
     space(v))

for i in 1:4
    component(P, i)[(:,0)] .= component(γ, i)
    component(P, i)[(:,1)] .= component(v, i)
end

### Tapes 
DF_ = LinearOperator(s_mani, s_mani, zeros(ComplexF64, 4*(2*N_F+1)*(N_T+1), 4*(2*N_F+1)*(N_T+1)))
F_ = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

### Vector field manifold
function f(P, parameters, N_F, N_T)

    s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4
    f=Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))
    
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

    project!(component(Df, 2, 1), Multiplication(A*P₃ + s*P₁*P₁-E))
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


function F!(F, P, para, N_F, N_T,v,λ)

    project!(F, -f(P, para, N_F, N_T) + Derivative(1,0)*P)

    γ = project(
        Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
        space(v))

      for i = 1:4
          component(F, i) .= component(F, i) +  λ*component(Derivative(0,1)*P,i)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) 
          component(F, i)[(:,0)] .= component(P, i)[(:,0)] .- component(γ, i)
          component(F, i)[(:,1)] .= component(P, i)[(:,1)] .- component(v, i)
      end
 
    return F
end


    function DF!(DF, P, para, N_F, N_T,λ)
  
        DF .= 0
        Df = DF

        Df = Df!(Df, P, para)
        
        project!(component(DF, 1, 1), -component(Df, 1, 1) + Derivative(1,0) + λ*project(Derivative(0,1), domain(component(Df, 1, 1)), domain(component(Df,1,1)), ComplexF64))
        project!(component(DF, 1, 2), -component(Df, 1,2))
        project!(component(DF, 1, 3), -component(Df, 1,3))
        project!(component(DF, 1, 4), -component(Df, 1,4))
        project!(component(DF, 2, 1), -component(Df, 2,1))
        project!(component(DF, 2, 2), -component(Df, 2,2) + Derivative(1,0) + λ*project(Derivative(0,1), domain(component(Df, 2, 2)), domain(component(Df,2,2)), ComplexF64))
        project!(component(DF, 2, 3), -component(Df, 2,3))
        project!(component(DF, 2, 4), -component(Df, 2,4))
        project!(component(DF, 3, 1), -component(Df, 3,1))
        project!(component(DF, 3, 2), -component(Df, 3,2))
        project!(component(DF, 3, 3), -component(Df, 3,3) + Derivative(1,0) + λ*project(Derivative(0,1), domain(component(Df, 3, 3)), domain(component(Df,3,3)), ComplexF64))
        project!(component(DF, 3, 4), -component(Df, 3,4))
        project!(component(DF, 4, 1), -component(Df, 4,1))
        project!(component(DF, 4, 2), -component(Df, 4,2))
        project!(component(DF, 4, 3), -component(Df, 4,3))
        project!(component(DF, 4, 4), -component(Df, 4,4) + Derivative(1,0) + λ*project(Derivative(0,1), domain(component(Df, 4, 4)), domain(component(Df,4,4)), ComplexF64))

         for i = 1:4
             component(DF,i,i)[(:,0),(:,0)] = I(2*N_F+1)
             component(DF,i,i)[(:,1),(:,1)] = I(2*N_F+1)
        end

        return DF
    end


### Testing
      F!(F_, P, para, N_F, N_T, v,λ)

#     DF!(DF_, P, para, N_F, N_T,λ)
#    x_0 = newton!((F_, DF_, x) -> (F!(F_, x, para, N_F, N_T,v,λ), DF!(DF_, x, para, N_F, N_T,λ)), P)
