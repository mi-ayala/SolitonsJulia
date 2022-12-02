using .SolitonsJulia
using RadiiPolynomial
using DifferentialEquations
using JLD2


### Finding Solitons - Shooting

### Modes
N_F = 10
N_T = 5


E = 1.0
A = -1.0
s = 1.0

### Parameters
parameters = soliton_parameters()


### Bundle
λ, v = get_bundle(N_F, parameters)


 γ = project(
     Sequence(Fourier(2, 1.0)^4, [zeros(5); zeros(5); [0.5, 0, 0, 0, 0.5]; -2 * [0.5im, 0, 0, 0, -0.5im]]),
     space(v))


### Manifold

function jac_vectorfield(N, A, E)
    Df_ = LinearOperator(Fourier(N,1.0)^4, Fourier(N,1.0)^4, zeros(ComplexF64, 4*(2*N+1), 4*(2*N+1)))
    cos2t = Sequence(Fourier(2, 1.0), [.5, 0, 0, 0, .5])

    project!(component(Df_, 1, 2), I)
    project!(component(Df_, 2, 1), Multiplication(A*cos2t - E))
    project!(component(Df_, 3, 4), I)
    project!(component(Df_, 4, 3), -4*I)
    return Df_
end

Df = jac_vectorfield(N_F, A, E)

function f_hat(P_view, A, E, α)
    f_hat_ = Sequence(space(P_view)[1][1]^4, zeros(ComplexF64, 4dimension(space(P_view)[1][1])))
    result = A* component(P_view, 1)*component(P_view, 3) + component(P_view, 1)^3
    h = result[(:,α)] # all fourier mode of α-th order
    h_ = Sequence(space(result)[1], h) # turn it into the Fourier sequence
    project!(component(f_hat_, 2), h_)
    return f_hat_
end


s_mani = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^4
P = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

γ = project(
    Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
    space(v))

for i in 1:4
    component(P, i)[(:,0)] .= component(γ, i)
    component(P, i)[(:,1)] .= component(v, i)
end

for α ∈ 2:N_T
    P_view = Sequence((Fourier(N_F, 1.0) ⊗ Taylor(α))^4, [view(component(P, 1), (:,0:α)); view(component(P, 2), (:,0:α)) ; view(component(P, 3), (:,0:α)) ; view(component(P, 4), (:,0:α))])

    g = (Derivative(1) - Df + λ*α*I) \ f_hat(P_view, A, E, α)

    for i in 1:4
        component(P, i)[(:,α)] .= component(g, i)
    end
end







function checking(P, t_final, N, σ)

    t_range = range(0, stop=2π, length = N)

    u₀= []
    
    for t ∈ t_range 

        push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

    end

    ### Finding Solitons.


    integration_ending = zeros(N,4)
   
    for i ∈  1:N

        tspan = (0, t_final)
        prob = ODEProblem(vectorField!,u₀[i],tspan,parameters)
        sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14)

        integration_ending[i,:] = sol[end]
    end

    return integration_ending

end   

t_range = range(0, stop=2π, length = 50)
σ = .9
    u₀= []
    
    for t ∈ t_range 

        push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

    end