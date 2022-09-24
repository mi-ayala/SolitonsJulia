using RadiiPolynomial, UnPack

N=3
N_T=3

p = parameters()
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

# project!(P, -f!(f, P, parameters) + (Derivative(1,0) + λ*Derivative(0,1))*P*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1])
#  )

# f!(f, P, p)

#(Derivative(1,0)*P + Derivative(0,1)*P*λ)*Sequence(Fourier(0, 1.0) ⊗ Taylor(1), [0, 1])


Derivative(0,1)*P*Sequence((Fourier(0, 1.0) ⊗ Taylor(1))^4, [0,1,0,1,0,1,0,1])