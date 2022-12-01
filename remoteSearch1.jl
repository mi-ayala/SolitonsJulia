using SolitonsJulia
using RadiiPolynomial
using DifferentialEquations
using JLD2
using ForwardDiff

### Finding Solitons - Shooting

### Modes
N_F = 15
N_T = 10


### Parameters
parameters = soliton_parameters()


### Bundle
λ, v = get_bundle(N_F, parameters)


### Periodic Orbit
# γ = project(
#     Sequence(Fourier(2, 1.0)^4, [zeros(5); zeros(5); [0.5, 0, 0, 0, 0.5]; -2 * [0.5im, 0, 0, 0, -0.5im]]),
#     space(v))


file = jldopen("ManifoldFD.jld2");
P = file["P"];


initial_condition = []  
N = 1000
σ = .9


function shooting!(initial_condition,t_final, N, σ)


        t_range = range(0, stop=2π, length = N)

        
        u₀= []
        
        for t ∈ t_range 

            push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

        end

        ### Finding Solitons.

        function condition(u,t,integrator)
            u[2]
        end
        
       
        for i ∈  1:N

            function affect!(integrator) 
                

                if abs(integrator.u[4]) < 1e-4
                    push!(initial_condition, u₀[i] )
                end    


            end

            cb = ContinuousCallback(condition,affect!)

            tspan = -(0, t_final)
            prob = ODEProblem(vectorField!,u₀[i],tspan,parameters)
            sol = solve(prob, VCABM(),abstol = 1e-14, reltol = 1e-14, callback=cb)
            #push!(bump, sol[end])

        end



end    