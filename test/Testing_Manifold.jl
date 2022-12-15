####################################
### Testing
### Testing manifold boundary points
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
### 2. 
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
using DifferentialEquations
using JLD2

file = jldopen("Manifold_Bundle.jld2")
P = file["P"]
v = file["v"]
λ = file["λ"]


### Modes
N_F = 17
N_T = 16

### Parameters
E = 1
A = -1
s = 1

parameters = soliton_parameters()


function shooting!(initial_condition,P,t_final, N, σ)


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
                

                if abs(integrator.u[4]) < 1e-5
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


### Building boundary points
N = 1000
t_range = range(0, stop=2π, length = N)
u₀= []
σ = .9

for t ∈ t_range 

    push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

end

### Checking 
### checking(P, -10, N, σ)

### Finding
ic = []
shooting!(ic,P,20, N, σ)