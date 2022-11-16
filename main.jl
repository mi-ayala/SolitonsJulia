using SolitonsJulia
using RadiiPolynomial
using DifferentialEquations, Plots


### Finding Solitons

### Modes
N_F = 20
N_T = 20


### Parameters
parameters = soliton_parameters()


### Bundle
λ, v = get_bundle(N_F, parameters)


### Periodic Orbit
γ = project(
    Sequence(Fourier(2, 1.0)^4, [zeros(5); zeros(5); [0.5, 0, 0, 0, 0.5]; -2 * [0.5im, 0, 0, 0, -0.5im]]),
    space(v))


### Manifold 
P = get_manifold(N_F, N_T, parameters, λ, v, γ );

N = 100
### Testing manifold
 t_range = range(0, stop=2π, length = N)


 σ_range = -1:1
 σ = .9

 u₀ = []
for t ∈ t_range 

    push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

end
    

cb = VectorContinuousCallback(condition,affect!,2)
initial_condition = []  
bump = []

for i ∈  1:N

    u0 = u₀[i]
    prob = ODEProblem(vectorField!,u0,tspan,p)  
    tspan = -(0, 10)
    prob = ODEProblem(vectorField!,u0,tspan,p)
    sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13, callback=cb)
    pues!(bump, sol[end] )

end


#plot(sol,idxs=(0,1))


function condition(out,u,t,integrator) # Event when event_f(u,t) == 0
    out[1] = u[2]
    out[2] = u[4]
  end
  
  function affect!(integrator, idx)
      push!(initial_condition, u0 )
      terminate!(integrator)
  end