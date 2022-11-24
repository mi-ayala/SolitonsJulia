using SolitonsJulia
using RadiiPolynomial
using DifferentialEquations
using JLD2

### Finding Solitons - Shooting

### Modes
N_F = 35
N_T = 30


### Parameters
parameters = soliton_parameters()


# ### Bundle
# λ, v = get_bundle(N_F, parameters)


# ### Periodic Orbit
# γ = project(
#     Sequence(Fourier(2, 1.0)^4, [zeros(5); zeros(5); [0.5, 0, 0, 0, 0.5]; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#     space(v))


# ### Manifold 
# P = get_manifold(N_F, N_T, parameters, λ, v, γ );

# 	file = jldopen("Manifold.jld2");
# 	P = file["P"];

# jldsave("Manifold.jld2"; P)

file = jldopen("Manifold.jld2");
P = file["P"];

N = 10000

### Manifold points
t_range = range(0, stop=2π, length = N)


 σ_range = -1:1
 σ = 1

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


for i ∈  1:N

    u₀[i]

    function affect!(integrator) 
        terminate!(integrator)
    end

    cb = ContinuousCallback(condition,affect!)

    tspan = -(0, 10)
    prob = ODEProblem(vectorField!,u₀[i],tspan,parameters)
    sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14, callback=cb)
    push!(bump, sol[end] )

end


for i ∈  1:N

    if abs(bump[i][4]) < 1e-3
        push!(initial_condition, u₀[i] )
    end    

end



function search(x,  parameters)
	
    σ = x[1]
    L = x[2]
    t = 1
    
    u0 = [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])]

   tspan = -(0, L)
   prob = ODEProblem(vectorField!,u0,tspan,parameters)
   sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-14)
   out = [sol[end][2], sol[end][4] ]

    return out
      
end

DF = x -> ForwardDiff.jacobian(x -> search(x,parameters) , [.5 1])
















# σ = -.9

 
# for t ∈ t_range 

#     push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

# end
    
  
# p=0

# for i ∈  1:N

#     u₀[i]

#     function affect!(integrator) 
#         terminate!(integrator)
#     end

#     cb = ContinuousCallback(condition,affect!)

#     tspan = -(0, 5)
#     prob = ODEProblem(vectorField!,u₀[i],tspan,parameters)
#     sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14, callback=cb)
#     push!(bump, sol[end] )

# end


# for i ∈  1:N

#     if abs(bump[i][4]) < 1e-4
#         push!(initial_condition, u₀[i] )
#     end    

# end


# ### Plots
# tspan = -(0, 4)
# prob = ODEProblem(vectorField!,initial_condition[1],tspan,p)
# sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14)
# fig = plot(sol.t[:] , sol[1,:])


# # for i ∈  2:length(initial_condition)

# #     tspan1 = -(0, 5)
# #     prob1 = ODEProblem(vectorField!,initial_condition[i],tspan1,p)
# #     sol1 = solve(prob1, VCABM(),abstol = 1e-13, reltol = 1e-13)
 
# #     plot!(fig, sol1.t[:] , sol1[1,:])

# # end