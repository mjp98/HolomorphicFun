
# Cut

struct Cut{T,S,R} <: AbstractSingularity
    domain::T
    Δ::S
    Σ::R
end

removedsingularity(f::Cut,z) = (f.Σ(z) + im*hilbert(f.Δ,z))/2

function truncate(cut::Cut{<:Ray},cutoff)
    d = truncate(domain(cut),cutoff)
    cΔ = Fun(Δ(cut),setdomain(space(Δ(cut)),d))
    cΣ = Fun(Σ(cut),setdomain(space(Σ(cut)),d))
    Cut(d,cΔ,cΣ)
end
