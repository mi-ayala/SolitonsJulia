#The way we do this is we simply write the output to the 1st input of the function. 
# For example, our Lorenz equation problem would be defined by the function:

using SolitonsJulia
using DifferentialEquations
using Plots

function vectorField!(du,u,p,t)

    E = 1.0
    A = -3.0
    s = 1.0

    du[1] = u[2]
    du[2] = -E*u[1] + A*u[3]*u[1] + s*u[1]*u[1]*u[1]
    du[3] = u[4]
    du[4] = -4.0*u[3]
end
   
   # and then we can use this function in a problem:
   
   u0 = [ 0.909917938392809;  -1.076582247110947;  -0.216389906598672;  -1.952614051288393]
   tspan = (0,-2.5)
   prob = ODEProblem(vectorField!,u0,tspan)
   sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13)
   # Using the plot recipe tools defined on the plotting page, we can choose to do a 3D phase space plot between the different variables:
   
   plot(sol,idxs=(0,1))



   
