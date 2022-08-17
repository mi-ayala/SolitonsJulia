
"""
f_bundle(x)

Tangent bundle vector field
```
"""
function f_bundle!(f, v, p)

    @unpack E, A, s = p   
    v₁, v₂ = eachcomponent(v)
    
    ### Cos2t Fourier Series
    cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

    project!(component(f, 1),  - v₂ )
    project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
    return f
end


function Df_bundle!(Df, v, p)

    @unpack E, A, s = p

    v₁, v₂, v₃ = eachcomponent(v)
    project!(component(Df, 1, 1), Multiplication(-σ*one(v₁)))
    project!(component(Df, 1, 2), Multiplication(σ*one(v₂)))
    project!(component(Df, 2, 1), Multiplication(ρ-v₃))
    project!(component(Df, 2, 2), Multiplication(-one(v₂)))
    return Df
end