using DifferentialEquations, Plots
#### "Soliton" from MATLAB

function condition(out,u,t,integrator) # Event when event_f(u,t) == 0
  out[1] = u[2]
  out[2] = u[4]
end

event_idx = []  # global variable

function affect!(integrator, idx)

    push!(event_idx, u0 )
    terminate!(integrator)
end

cb = VectorContinuousCallback(condition,affect!,2)


u0 = [ 0.909917938392809;  -1.076582247110947;  -0.216389906598672;  -1.952614051288393]
tspan = -(0, 10)
p=0.0
prob = ODEProblem(vectorField!,u0,tspan,p)  
tspan = -(0, 10)
prob = ODEProblem(vectorField!,u0,tspan,p)
# sol = solve(prob,Tsit5(),callback=cb,dt=1e-3,adaptive=false)
sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13, callback=cb)
println(sol.retcode)
println(event_idx) # [1]
plot(sol,idxs=(0,1))
### sol[end] is the bump.