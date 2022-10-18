
# function F_bundle!(v, p)

#     @unpack E, A, s = p   
#     v₁, v₂ = eachcomponent(v)
    
#     ### Cos2t Fourier Series
#     cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

#     project!(component(f, 1),  - v₂ )
#     project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
#     return F
# end


# function G_bundle!(v, p)

#     @unpack E, A, s = p   
#     v₁, v₂ = eachcomponent(v)
    
#     ### Cos2t Fourier Series
#     cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

#     project!(component(f, 1),  - v₂ )
#     project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
#     return G
# end


# function B!(v, p)

#     @unpack E, A, s = p   
#     v₁, v₂ = eachcomponent(v)
    
#     ### Cos2t Fourier Series
#     cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

#     project!(component(f, 1),  - v₂ )
#     project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
#     return F

# end


# function Df(N, A, E)
#     Df_ = LinearOperator(Fourier(N,1.0)^4, Fourier(N,1.0)^4, zeros(ComplexF64, 4*(2*N+1), 4*(2*N+1)))
#     cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

#     project!(component(Df_, 1, 2), I)
#     project!(component(Df_, 2, 1), Multiplication(A*cos2t - E))
#     project!(component(Df_, 3, 4), I)
#     project!(component(Df_, 4, 3), -4*I)
#     return Df_
# end


# function F(x, Df)

#     λ = component(x,1)
#     v =  component(x,2)

#     F_ = similar(x)

#     F_[1] = component(v,1)(0)[0] - 1 
#     project!(component(F_, 2), (Derivative(1) - Df + λ*I)*v  )

#     return F_
# end


# function DF(x, Df, N)

#     λ = component(x,1)
#     v =  component(x,2)

#     DF_ = LinearOperator(ParameterSpace() × Fourier(N,1.0)^4, ParameterSpace() × Fourier(N,1.0)^4, zeros(ComplexF64, 4*(2*N+1) + 1, 4*(2*N+1) +1))

#    project!(componenet(component(DF_,1,2),1), Evaluation(0))
#    project!(component(DF_,2,1), v)
#    project!(component(DF_,2,2), Derivative(1) - Df + λ*I)
   
   
#     return DF_
# end


# Or more general





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

# # compute manifold

# function f_hat(P_view, A, E, α)
#     f_hat_ = Sequence(space(P_view)[1][1]^4, zeros(ComplexF64, 4dimension(space(P_view)[1][1])))
#     result = A* component(P_view, 1)*component(P_view, 3) + component(P_view, 1)^3
#     h = result[(:,α)] # all fourier mode of α-th order
#     h_ = Sequence(space(result)[1], h) # turn it into the Fourier sequence
#     project!(component(f_hat_, 2), h_)
#     return f_hat_
# end

# N_T = 10
# s_mani = (Fourier(N, 1.0) ⊗ Taylor(N_T))^4
# P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

# γ = project(
#     Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#     space(v))

# for i in 1:4
#     component(P, i)[(:,0)] .= component(γ, i)
#     component(P, i)[(:,1)] .= component(v, i)
# end

# for α ∈ 2:N_T
#     P_view = Sequence((Fourier(N, 1.0) ⊗ Taylor(α))^4, [view(component(P, 1), (:,0:α)); view(component(P, 2), (:,0:α)) ; view(component(P, 3), (:,0:α)) ; view(component(P, 4), (:,0:α))])

#     g = (Derivative(1) - Df + λ*α*I) \ f_hat(P_view, A, E, α)

#     for i in 1:4
#         component(P, i)[(:,α)] .= component(g, i)
#     end
# end

# using GLMakie
# t_range = range(0, stop=2π, length = 100)
# σ_range = -1:0.1:1
# surface([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range, σ in σ_range])



# for σ = σ_range
#     lines!([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range];
#         transparency = true, linewidth = 3)
# end


# surface([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 2)(t, σ)[(0,0)]) for t in t_range, σ in σ_range])



# for σ = σ_range
#     lines!([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range])
# end



### ANother way of computing manifold
# using SolitonsJulia
# using UnPack

# ### Parameters
# N_F = 10
# A = -3
# E = 1
# para = soliton_parameters()

# ### Bundle
# λ, v = get_bundle(N_F, A, E)

# ### Manifold initial guess
# N_T = 4
# s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4

# ### Periodic orbit data
# P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

#  γ = project(
#      Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#      space(v))
 
# for i in 1:4
#     component(P, i)[(:,0)] .= component(γ, i)
#     component(P, i)[(:,1)] .= component(v, i)
# end

# ### Tapes 
# DF_ = LinearOperator(s_mani, s_mani, zeros(ComplexF64, 4*(2*N_F+1)*(N_T+1), 4*(2*N_F+1)*(N_T+1)))
# F_ = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

# ### Vector field manifold
# function f(P, parameters, N_F, N_T)

#     s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4
#     f=Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))
    
#     @unpack E, A, s = parameters   
#     P₁, P₂, P₃, P₄ = eachcomponent(P)

#     project!(component(f, 1), P₂ )
#     project!(component(f, 2), - E*P₁ + A*P₃*P₁ + s*P₁*P₁*P₁ )
#     project!(component(f, 3), + P₄)
#     project!(component(f, 4), - 4*P₃)

#     return f

# end

# function Df!(Df, P, parameters)

#     @unpack E, A, s = parameters
#     P₁, P₂, P₃, P₄ = eachcomponent(P)

#     project!(component(Df, 1, 1), 0*I) #Multiplication(zero(P₁))
#     project!(component(Df, 1, 2), Multiplication(one(P₂)))
#     project!(component(Df, 1, 3), Multiplication(zero(P₃)))
#     project!(component(Df, 1, 4), Multiplication(zero(P₄)))

#     project!(component(Df, 2, 1), Multiplication(A*P₃ + s*P₁*P₁-E))
#     project!(component(Df, 2, 2), Multiplication(zero(P₂)))
#     project!(component(Df, 2, 3), Multiplication(A*P₁))
#     project!(component(Df, 2, 4), Multiplication(zero(P₄)))

#     project!(component(Df, 3, 1), Multiplication(zero(P₁)))
#     project!(component(Df, 3, 2), Multiplication(zero(P₂)))
#     project!(component(Df, 3, 3), Multiplication(zero(P₃)))
#     project!(component(Df, 3, 4), Multiplication(one(P₄)))

#     project!(component(Df, 4, 1), Multiplication(zero(P₁)))
#     project!(component(Df, 4, 2), Multiplication(zero(P₂)))
#     project!(component(Df, 4, 3), Multiplication(-4*one(P₃)))
#     project!(component(Df, 4, 4), Multiplication(zero(P₄)))

#     return Df
# end


# function F!(F, P, para, N_F, N_T,v,λ)

#     project!(F, -f(P, para, N_F, N_T) + Derivative(1,0)*P)

#     γ = project(
#         Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#         space(v))

#       for i = 1:4
#           component(F, i) .= component(F, i) +  λ*component(Derivative(0,1)*P,i)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) 
#           component(F, i)[(:,0)] .= component(P, i)[(:,0)] .- component(γ, i)
#           component(F, i)[(:,1)] .= component(P, i)[(:,1)] .- component(v, i)
#       end
 
#     return F
# end


#     function DF!(DF, P, para, N_F, N_T,λ)
  
#         DF .= 0
#         Df = copy(DF)

#         Df = Df!(Df, P, para)

#         M = project(Multiplication(Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) ), domain(component(Df, 1, 1)), codomain(component(Df,1,1)))
        
#         project!(component(DF, 1, 1), -component(Df, 1, 1) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 1, 1)), codomain(component(Df,1,1)), ComplexF64))
#         project!(component(DF, 1, 2), -component(Df, 1,2))
#         project!(component(DF, 1, 3), -component(Df, 1,3))
#         project!(component(DF, 1, 4), -component(Df, 1,4))
#         project!(component(DF, 2, 1), -component(Df, 2,1))
#         project!(component(DF, 2, 2), -component(Df, 2,2) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 2, 2)), codomain(component(Df,2,2)), ComplexF64))
#         project!(component(DF, 2, 3), -component(Df, 2,3))
#         project!(component(DF, 2, 4), -component(Df, 2,4))
#         project!(component(DF, 3, 1), -component(Df, 3,1))
#         project!(component(DF, 3, 2), -component(Df, 3,2))
#         project!(component(DF, 3, 3), -component(Df, 3,3) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 3, 3)), codomain(component(Df,3,3)), ComplexF64))
#         project!(component(DF, 3, 4), -component(Df, 3,4))
#         project!(component(DF, 4, 1), -component(Df, 4,1))
#         project!(component(DF, 4, 2), -component(Df, 4,2))
#         project!(component(DF, 4, 3), -component(Df, 4,3))
#         project!(component(DF, 4, 4), -component(Df, 4,4) + Derivative(1,0) + λ*M*project(Derivative(0,1), domain(component(Df, 4, 4)), codomain(component(Df,4,4)), ComplexF64))

#          for i = 1:4
#              component(DF,i,i)[(:,0),(:,0)] = I(2*N_F+1)
#              component(DF,i,i)[(:,1),(:,1)] = I(2*N_F+1)
#         end

#         return DF
#     end


# ### Testing
#       F!(F_, P, para, N_F, N_T, v,λ)

#     DF!(DF_, P, para, N_F, N_T,λ)
# newton!((F_, DF_, x) -> (F!(F_, x, para, N_F, N_T,v,λ), DF!(DF_, x, para, N_F, N_T,λ)), P)

# component(P, 1)
