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


