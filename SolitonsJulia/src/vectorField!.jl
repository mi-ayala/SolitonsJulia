
# export condition, affect!, initial_condition, u0

function vectorField!(du,u,p,t)

    E = 1.0
    A = -3.0
    s = 1.0

    du[1] = u[2]
    du[2] = -E*u[1] + A*u[3]*u[1] + s*u[1]*u[1]*u[1]
    du[3] = u[4]
    du[4] = -4.0*u[3]
end

# function condition(out,u,t,integrator) # Event when event_f(u,t) == 0
#     out[1] = u[2]
#     out[2] = u[4]
#   end
  
 
#   initial_condition = []  
#   function affect!(integrator, idx)
#       push!(initial_condition, u0 )
#       terminate!(integrator)
#   end