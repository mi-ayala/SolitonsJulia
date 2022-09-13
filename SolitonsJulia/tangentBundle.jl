
function F_bundle!(v, p)

    @unpack E, A, s = p   
    v₁, v₂ = eachcomponent(v)
    
    ### Cos2t Fourier Series
    cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

    project!(component(f, 1),  - v₂ )
    project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
    return F
end


function G_bundle!(v, p)

    @unpack E, A, s = p   
    v₁, v₂ = eachcomponent(v)
    
    ### Cos2t Fourier Series
    cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

    project!(component(f, 1),  - v₂ )
    project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
    return G
end


function B!(v, p)

    @unpack E, A, s = p   
    v₁, v₂ = eachcomponent(v)
    
    ### Cos2t Fourier Series
    cos2t = Sequence(Fourier(2, 1.0), [.5 0 0 0 .5])

    project!(component(f, 1),  - v₂ )
    project!(component(f, 2), + E*v₁ - A*(v₁*cos2t))
    
    return F

end
