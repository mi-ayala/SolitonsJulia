using SolitonsJulia
using RadiiPolynomial

### Modes
N_F = 20
N_T = 20

### Parameters
parameters = soliton_parameters()

### Bundle
λ, v = get_bundle(N_F, parameters)

### Periodic Orbit
γ = project(
    Sequence(Fourier(2, 1.0)^4, [zeros(5); zeros(5); [0.5, 0, 0, 0, 0.5]; -2 * [0.5im, 0, 0, 0, -0.5im]]),
    space(v))

### Manifold 
P = get_manifold(N_F, N_T, parameters, λ, v, γ );

### Testing manifold
 t_range = range(0, stop=2π, length = 10)


 σ_range = .9

 u₀ = []
for t ∈ t_range 

    push!(u₀ , [real(component(P, 1)(t, σ)[(0,0)]), real(component(P, 2)(t, σ)[(0,0)]), real(component(P, 3)(t, σ)[(0,0)]), real(component(P, 4)(t, σ)[(0,0)])])

end
    
#  surface([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#      [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#      [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#      colormap = ColorSchemes.BrBG_10.colors,)


for i ∈ 1:10

    tspan = (0, 1)
    prob = ODEProblem(vectorField!,u₀[i],tspan)
    sol = solve(prob, VCABM(),abstol = 1e-13, reltol = 1e-13)
    # Using the plot recipe tools defined on the plotting page, we can choose to do a 3D phase space plot between the different variables:
    
    plot(sol,idxs=(0,1))
    

end


### Plotting manifold

# using GLMakie
# using ColorSchemes 

# t_range = range(0, stop=2π, length = 100)
# σ_range = -1:0.1:1
# surface([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     colormap = ColorSchemes.BrBG_10.colors,)



# for σ = σ_range
#     lines!([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range];
#         transparency = true, linewidth = 3)
# end


# surface([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 2)(t, σ)[(0,0)]) for t in t_range, σ in σ_range])



# for σ = σ_range
#     lines!([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range])
# end

