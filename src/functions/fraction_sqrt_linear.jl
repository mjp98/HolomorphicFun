function singularities(f::ScalarF{ScalarL{T},Sqrt{true,S}}) where {T,S}
    l = numerator(f)
    s = denominator(f)
    cut = first(singularities(inv(s)))
    d = domain(cut)
    jsp = JacobiWeight(-0.5,-0.5,d)
    Δ = Fun(x->cut.Δ(x)*l(x),jsp)
    Σ = cut.Σ
    return [Cut(d,Δ,Σ)]
end

function singularities(f::ScalarF{Sqrt{true,S},ScalarL{T}}) where {T,S}
    s = numerator(f)
    l = denominator(f)
    cut = first(singularities(s))
    d = domain(cut)
    jsp = JacobiWeight(-0.5,-0.5,d)
    Δ = Fun(x->cut.Δ(x)/l(x),jsp)
    Σ = cut.Σ
    x = root(l)
    return Cut(d,Δ,Σ)+Pole(x,s(x)/l.a,0)
end
