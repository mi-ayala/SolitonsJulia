using .SolitonsJulia
using RadiiPolynomial
using DifferentialEquations, Plots


### Finding Solitons - Shooting

### Modes
N_F = 35
N_T = 30


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

N = 25000

### Manifold points
t_range = range(0, stop=2π, length = N)


 σ_range = -1:1
 σ = .9

 u₀= []
 
for t ∈ t_range 

    push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

end
    


### Finding
function condition(u,t,integrator)
    u[2]
 end
  
initial_condition = []  
bump = []
p=0

for i ∈  1:N

    u₀[i]

    function affect!(integrator) 
        terminate!(integrator)
    end

    cb = ContinuousCallback(condition,affect!)

    tspan = -(0, 5)
    prob = ODEProblem(vectorField!,u₀[i],tspan,p)
    sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14, callback=cb)
    push!(bump, sol[end] )

end


for i ∈  1:N

    if abs(bump[i][4]) < 1e-4
        push!(initial_condition, u₀[i] )
    end    

end



### Plots
tspan = -(0, 4)
prob = ODEProblem(vectorField!,initial_condition[1],tspan,p)
sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14)
fig = plot(sol.t[:] , sol[1,:])


# for i ∈  2:length(initial_condition)

#     tspan1 = -(0, 5)
#     prob1 = ODEProblem(vectorField!,initial_condition[i],tspan1,p)
#     sol1 = solve(prob1, VCABM(),abstol = 1e-13, reltol = 1e-13)
 
#     plot!(fig, sol1.t[:] , sol1[1,:])

# end