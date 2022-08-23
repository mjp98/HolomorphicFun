"""
    P = ScalarPower(z,n,θ,s) defines a function

    P(α) = s*(α-z)^n

    with

    1) P(z+ε) ∼ s*ε^n for ϵ>0
    2) branch cut parallel to θ

    i.e.

    P(α) = s*((θ*(α-z))^n)/((θ)^n)

    Computed as scale*((θ*(α-z))^n) with scale = (s/((θ)^n))

"""
struct ScalarPower{T} <: AbstractScalarFunction
    scale::T
    center::T
    exponent::T
    θ::T
    function ScalarPower(z,n,θ,s)
        s /= ((θ)^n)
        z = cfloat(z)
        new{eltype(promote(z,n,θ,s))}(promote(z,n,θ,s)...)
    end
end

scale(K::ScalarPower) = K.scale
center(K::ScalarPower) = K.center
exponent(K::ScalarPower) = K.exponent
angle(K::ScalarPower) = K.θ

function (K::ScalarPower)(α)
    s = scale(K)
    z = center(K)
    n = exponent(K)
    θ = angle(K)
    return s*((θ*(α-z))^n)/((θ)^n)
end
