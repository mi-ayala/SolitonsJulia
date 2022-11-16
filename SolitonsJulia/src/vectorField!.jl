function vectorField!(du,u,p,t)

    E = 1.0
    A = -3.0
    s = 1.0

    du[1] = u[2]
    du[2] = -E*u[1] + A*u[3]*u[1] + s*u[1]*u[1]*u[1]
    du[3] = u[4]
    du[4] = -4.0*u[3]
end

