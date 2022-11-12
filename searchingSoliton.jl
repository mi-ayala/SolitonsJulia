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
   
# function condition(u,t,integrator) # Event when event_f(u,t) == 0
#     u[2]
# end

condition1(u,t,integrator) = u[2]
affect!(integrator) = terminate!(integrator)
cb1 = ContinuousCallback(condition1,affect!)

condition2(u,t,integrator) = u[4]
affect!(integrator) = terminate!(integrator)
cb2 = ContinuousCallback(condition2,affect!)

   #### Testing shooting
   
   u0 = [ 0.909917938392809;  -1.076582247110947;  -0.216389906598672;  -1.952614051288393]
   tspan = -(0, 10)
   prob = ODEProblem(vectorField!,u0,tspan)


   #sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13)
   sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13, callback=cb1)

   # Using the plot recipe tools defined on the plotting page, we can choose to do a 3D phase space plot between the different variables:
   
   plot(sol,idxs=(0,1))


   ### Effect function with saving value














   


##Function below computes the Chebyshev nodes in [-1,1].
 
# function cheb_nodes(N::Integer)
#     req_nodes = cos.(pi*((0:1:N))/N);
#     return req_nodes
# end
 
# ##Function below gets the Chebyshev coefficients of the usual Chebyshev interpolant. Note that f0 needs to be your function evaluated at the output of cheb_nodes, and N needs to match your input from cheb_nodes.
 
# function cheb_interp(f₀,N::Integer)
#     a = zeros(eltype(f₀),N+1);
#     e₁ = exp.(im*pi*(0:1:N)/N);
#     e₂ = exp.(-im*pi*(0:1:N)/N);
#     T = similar(a);
#     for n in 0:N
#         T[:] = real.(e₁.^n + e₂.^n)/2;
#         T[1] = T[1]/2;  T[end] = T[end]/2;
#         if n==N
#             T[:] = T[:]/2;
#         end
#         a[n+1] = (T'*f₀)/N;
#     end
#     return a
# end


   



