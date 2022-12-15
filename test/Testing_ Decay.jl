####################################
### Testing
### Testing bundle/manifold decay
####################################

####################################
###  Things we are testing:
####################################
###  Decay of coefficients 
####################################

####################################
### Code Blueprint
####################################
### 1. Bundle with N_F=17, phase condition with .5 
### 2. Manifold inverting operator
### 3. Zero finding problem for the manifold with FD
### 4. 
### 5. 
####################################

####################################
### Observations after run
####################################
### 1 . Decay node 17 (4 components)
###    2.430815405772803e-22 - 2.063058758888647e-21im
###    3.913901216750381e-20 + 5.120173270428039e-21im
###    -1.472380965647887e-30 - 2.479032950290945e-31im
###    4.7618617874974406e-32 - 2.5979237654867842e-30im
####################################

####################################
### Notes
####################################
### 1 . Last two components are zero.
####################################


using SolitonsJulia
using RadiiPolynomial
using FiniteDifferences


### Modes
N_F = 17
N_T = 16

### Parameters
E = 1
A = -1
s = 1

### Parameters
parameters = soliton_parameters()

### Bundle
λ, v = get_bundle(N_F, parameters)

### First  bundle component decay
### component(v,1)

#################################
### Manifold inverting operator
#################################

# function jac_vectorfield(N, A, E)
#     Df_ = LinearOperator(Fourier(N,1.0)^4, Fourier(N,1.0)^4, zeros(ComplexF64, 4*(2*N+1), 4*(2*N+1)))
#     cos2t = Sequence(Fourier(2, 1.0), [.5, 0, 0, 0, .5])

#     project!(component(Df_, 1, 2), I)
#     project!(component(Df_, 2, 1), Multiplication(A*cos2t - E))
#     project!(component(Df_, 3, 4), I)
#     project!(component(Df_, 4, 3), -4*I)
#     return Df_
# end

# Df = jac_vectorfield(N_F, A, E)

# function f_hat(P_view, A, E, α)
#     f_hat_ = Sequence(space(P_view)[1][1]^4, zeros(ComplexF64, 4dimension(space(P_view)[1][1])))
#     result = A*component(P_view, 1)*component(P_view, 3) + component(P_view, 1)^3
#     h = result[(:,α)] # all fourier mode of α-th order
#     h_ = Sequence(space(result)[1], h) # turn it into the Fourier sequence
#     project!(component(f_hat_, 2), h_)
#     return f_hat_
# end


# s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4

# P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

# γ = project(
#     Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#     space(v))

# for i in 1:4
#     component(P, i)[(:,0)] .= component(γ, i)
#     component(P, i)[(:,1)] .= component(v, i)
# end

# for α ∈ 2:N_T
#     P_view = Sequence((Fourier(N_F, 1.0) ⊗ Taylor(α))^4, [view(component(P, 1), (:,0:α)); view(component(P, 2), (:,0:α)) ; view(component(P, 3), (:,0:α)) ; view(component(P, 4), (:,0:α))])

#     g = (Derivative(1) - Df + λ*α*I) \ f_hat(P_view, A, E, α)

#     for i in 1:4
#         component(P, i)[(:,α)] .= component(g, i)
#     end
# end

##############################
##### Manifold using autodiff
##############################

function F_manifold_vectorfield(P, space)
    E = 1
    A = -1
    s = 1

    f = Sequence(space, zeros(ComplexF64, dimension(space)))

    P₁, P₂, P₃, P₄ = eachcomponent(P)

    project!(component(f, 1), P₂)
    project!(component(f, 2), -E * P₁ + A * P₃ * P₁ + s * P₁ * P₁ * P₁)
    project!(component(f, 3), + P₄)
    project!(component(f, 4), -4 * P₃)

    return f

end

function F_manifold(P_input, v, λ, γ, space, f)


    P_real = P_input[1:dimension(space)]
    P_im = 1im*P_input[dimension(space)+1:end]

    P = Sequence(space, P_real + P_im)

    F_ = Sequence(space, zeros(ComplexF64, dimension(space)))

    project!(F_, -f(P, space) + Derivative(1, 0) * P)

    for i = 1:4
        component(F_, i) .= component(F_, i) + λ * component(Derivative(0, 1) * P, i) * Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1])
        component(F_, i)[(:, 0)] .= component(P, i)[(:, 0)] .- component(γ, i)
        component(F_, i)[(:, 1)] .= component(P, i)[(:, 1)] .- component(v, i)
    end

    return [real(F_[:]); imag(F_[:])]
end


s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4
P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

γ = project(
    Sequence(Fourier(2, 1.0)^4, [zeros(5); zeros(5); [0.5, 0, 0, 0, 0.5]; -2 * [0.5im, 0, 0, 0, -0.5im]]),
    space(v))
    
for i in 1:4
        component(P, i)[(:,0)] .= component(γ, i)
        component(P, i)[(:,1)] .= component(v, i)
end

P = [real(P[:]);imag(P[:])]

f_manifold = x -> F_manifold(x, v, λ, γ, s_mani, F_manifold_vectorfield)
df_manifold = x -> jacobian(central_fdm(8, 1), f_manifold, x)[1]

w =  newton(x -> (f_manifold(x), df_manifold(x)), P)
w = w[1]
w = w[1:dimension(s_mani)] + im*w[dimension(s_mani)+1:end]
P = Sequence(s_mani,w)


### Manifold first component a0 Fourier expansion 
### Last two coefficients when N_T=10
###   2.1311307173139004e-19 - 1.4806628564341287e-18im
###   1.6710254050046438e-31 + 2.0682860718537375e-31im
### Last two coefficients when N_T=16
###   2.2388971965058967e-25 + 4.145201206017607e-26im
###   2.1311319619674126e-19 + 1.4806629195431269e-18im
component(P,1)[(:,1)]

### ???
### component(P,1)[(2,:)]



 using JLD2

 @save "Manifold_Bundle_FD.jld2"  P v λ





