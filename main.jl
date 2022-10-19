using SolitonsJulia
using RadiiPolynomial

### Modes
N_F = 10
N_T = 5

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


### Plotting manifold

using GLMakie
t_range = range(0, stop=2π, length = 100)
σ_range = -1:0.1:1
surface([real(component(P, 1)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
    [real(component(P, 2)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
    [real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range])

for σ = σ_range
    lines!([real(component(P, 1)(t, σ)[(0,0)]) for t in t_range],
    [real(component(P, 2)(t, σ)[(0,0)]) for t in t_range],
    [real(component(P, 3)(t, σ)[(0,0)]) for t in t_range];
        transparency = true, linewidth = 3)
end

# surface([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range, σ in σ_range],
#     [real(component(P, 2)(t, σ)[(0,0)]) for t in t_range, σ in σ_range])

# for σ = σ_range
#     lines!([real(component(P, 3)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 4)(t, σ)[(0,0)]) for t in t_range],
#     [real(component(P, 1)(t, σ)[(0,0)]) for t in t_range])
# end



