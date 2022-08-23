
# Pole

struct Pole{T,S,R} <: AbstractSingularity
    domain::T
    Δ::S
    Σ::R
end

Pole(d::Number,Δ::Number,Σ::Number) = Pole(ScalarPole(d,Δ,Σ))

removedsingularity(f::Pole,z) = (f.Σ(z) + im*hilbert(f.Δ,z))

function *(a::Pole,b::Fun)
    ac = coefficients(Δ(a))
    bc = coefficients(b)

    @assert length(ac) == 1

    Δab = Δ(a)*first(bc)
    Σab = Σ(a)*b

    if length(bc) > 1
        Σab += Fun(space(Σ(a)),-first(ac)*bc[2:end]/(2π*im))
    end
    return Pole(domain(a),Δab,Σab)
end
