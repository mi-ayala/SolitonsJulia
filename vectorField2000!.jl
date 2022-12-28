"""
    vectorField2000!(du,u,parameters,t)


Compute the vector field given in Eq. (27) 
G. Alfimov, V.V. Konotop / Physica D 146 (2000) 307–327
"""


function vectorField2000!(du,u,parameters,t)

    omega = 1.05
    E2 =  -0.5

    E = omega^2
    A = -E2*E
    s = - (3/4)*E 

    du[1] = u[2]
    du[2] = -E*u[1] + A*u[3]*u[1] + s*u[1]*u[1]*u[1]
    du[3] = u[4]
    du[4] = -4.0*u[3]

end

