module SolitonsJulia

using RadiiPolynomial
using Parameters


export soliton_parameters,  get_bundle, get_manifold


include("get_bundle.jl")
include("get_manifold.jl")


### Parameter space
@with_kw struct soliton_parameters{R}
    E::R = 1
    A::R = -3
    s::R = 1
end







end