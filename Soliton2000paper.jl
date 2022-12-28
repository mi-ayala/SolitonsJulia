using SolitonsJulia
using DifferentialEquations
using Plots

include("vectorField2000!.jl")


#### Testing shooting
   
#    u0 = [ 0.71; 0;  1; 0]
#    tspan = (0, 50)
#    prob = ODEProblem(vectorField2000!,u0,tspan)

#    #sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13)
#    sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13)

#    # Using the plot recipe tools defined on the plotting page, we can choose to do a 3D phase space plot between the different variables:
   
   #plot(sol,idxs=(0,1))

   N=150

   initial_conditions = range(.69, stop=.72, length=N)
   tspan = (0, 50)
   candidates=[]

   for  i in 1:N

    u0 = [initial_conditions[i]; 0;  1; 0]

    prob = ODEProblem(vectorField2000!,u0,tspan)
    sol = solve(prob, Tsit5(),abstol = 1e-13, reltol = 1e-13)

    println(i)
    println(sol[end][1])
    if abs(sol[end][1]) < 1e-5
        push!(u0 , candidates)
    end    
    
   end