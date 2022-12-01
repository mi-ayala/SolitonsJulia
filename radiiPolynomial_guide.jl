### Parameters
using RadiiPolynomial

N_F = 1
N_T = 1
space = (Fourier(N_F, 1.0) ⊗ Taylor(N_T))^2

x = Sequence(space, zeros(ComplexF64, dimension(space)))

x₁ = component(x,1)

#### Fourier series for Taylor coefficient a0
x₁[(:,0)]

A = LinearOperator(space, space, zeros(ComplexF64, 2*(2*N_F+1)*(N_T+1), 2*(2*N_F+1)*(N_T+1)))

A₁ = component(A,1,1)

### Derivative first variable with respect to its variables
A₁.= rand(6,6)

### Taylor 0 - Taylor -0 / First fourier column
component(A,1,1)[(:,0),(-1,0)]

### Taylor 0 - Taylor -0 / First fourier row
component(A,1,1)[(-1,0),(:,0)]