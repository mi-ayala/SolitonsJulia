module SolitonsJulia

using RadiiPolynomial
using Parameters


export soliton_parameters, get_bundle, get_manifold, vectorField!


include("get_bundle.jl")
include("get_manifold.jl")
include("vectorField!.jl")


### Parameter space
@with_kw struct soliton_parameters{R}
    E::R = 1
    A::R = -1
    s::R = 1
end







end