using RadiiPolynomial, UnPack
using .SolitonsJulia

N=2
N_T=5

para = soliton_parameters()
λ = 1

### Tensor product components 
s_mani = (Fourier(N, 1.0) ⊗ Taylor(N_T))^4

P = Sequence(s_mani, ones(ComplexF64, dimension(s_mani)))
f = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))

#h = Sequence(space(result)[1], h) # Turn it into sequence
#Derivative(1,0)*P*Sequence(Taylor(1) ⊗ Fourier(0, 1.0), [0, 1])


function f!(f, P, parameters)

    @unpack E, A, s = parameters   
    P₁, P₂, P₃, P₄ = eachcomponent(P)

    project!(component(f, 1), P₂ )
    project!(component(f, 2), - E*P₁ + A*P₃*P₁ + s*P₁*P₁*P₁ )
    project!(component(f, 3), + P₄)
    project!(component(f, 4), - 4*P₃)

    return f

end

F = Sequence(s_mani, zeros(ComplexF64, dimension(s_mani)))
P_input = P

project!(F, -f!(f, P_input, para) + Derivative(1,0)*P_input)
 
 for i = 1:4
 component(F, i) = component(F, i) +  λ*component(Derivative(0,1)*P_input,i)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1]) 
 end
    
γ = project(
     Sequence(Fourier(2, 1.0)^4, [zeros(5) ; zeros(5) ; [.5, 0, 0, 0, .5] ; -2 * [0.5im, 0, 0, 0, -0.5im]]),
     Fourier(N, 1.0)^4)
v = γ


for i = 1:4
component(F, i)[(:,0)] .= component(P_input, i)[(:,0)] .- component(γ, i)
component(F, i)[(:,1)] .= component(P_input, i)[(:,1)] .- component(v, i)
end