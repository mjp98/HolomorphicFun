import Base: angle

function complexsqrt(x,y,θ=1)
    eiθ = cispi(1-θ)
    eiθ2 = cispi((1-θ)/2)
    return sqrt(eiθ*(x-y))/eiθ2
end

# -----------------------------
# Sqrt - shifted with branch cut at angle theta
# -----------------------------

export Sqrt

struct Sqrt{S,T} <: AbstractScalarFunction
    center::Complex{T}
    angle::T
    scale::Complex{T}
    function Sqrt{S}(center,angle::T,scale=one(complex(T))) where {S,T}
        c,a,s = promote(float(center),angle,scale)
        new{S,float(T)}(complex(c),real(a),complex(s))
    end
end
center(x::Sqrt) = x.center
scale(x::Sqrt) = x.scale
angle(x::Sqrt) = x.angle

evaluate(x::Sqrt{true},z) = scale(x)*complexsqrt(z,center(x),angle(x))
evaluate(x::Sqrt{false},z) = scale(x)/complexsqrt(z,center(x),angle(x))

*(x::Sqrt{T},y::Union{Number,ScalarC}) where T = Sqrt{T}(center(x),angle(x),scale(x)*y)
/(x::Sqrt{T},y::Union{Number,ScalarC}) where T = Sqrt{T}(center(x),angle(x),scale(x)/y)

*(x::Union{Number,ScalarC},y::Sqrt{T}) where T = Sqrt{T}(center(y),angle(y),x*scale(y))
/(x::Union{Number,ScalarC},y::Sqrt{T}) where T = Sqrt{!T}(center(y),angle(y),x/scale(y))
inv(y::Sqrt{T,S}) where {T,S} = Sqrt{!T}(center(y),angle(y),inv(scale(y)))

-(x::Sqrt{T,S}) where {T,S} = Sqrt{T}(center(x),angle(x),-scale(x))

function singularities(x::Sqrt{T}) where T
    d = Ray{angle(x)}(center(x))
    csp = Chebyshev(d)
    jsp = JacobiWeight(-0.5,-0.5,d)
    λ = T  ? -2*scale(x)*cispi(angle(x)/2) : -2*scale(x)*cispi(-angle(x)/2)
    Σ = Fun(zero, csp)
    Δ = T ? Fun(jsp,[λ; λ]) : Fun(jsp,[λ; -λ])
    return [Cut(d,Δ,Σ)]
end
