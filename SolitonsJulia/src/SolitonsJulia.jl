module SolitonsJulia

using RadiiPolynomial, Parameters

export soliton_parameters, cos2t

### Parameter space
@with_kw struct soliton_parameters{R}
    E::R = 1
    A::R = -3
    s::R = 1
end

### Fourier series for cosine
cos2t = Sequence(Fourier(2, 1.0), [.5, 0, 0, 0, .5])



end 