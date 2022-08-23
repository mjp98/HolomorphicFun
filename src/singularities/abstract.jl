
abstract type AbstractSingularity <: AbstractDefinedFunction end

domain(x::AbstractSingularity) = x.domain
Δ(x::AbstractSingularity) = x.Δ
Σ(x::AbstractSingularity) = x.Σ

cauchy(x::AbstractSingularity,z) = cauchy(Δ(x),z)
evaluate(x::AbstractSingularity,z) = cauchy(x,z)

# Vector of singularities

cauchy(f::AbstractVector{<:AbstractSingularity},z) = isempty(f) ? zero(z) : sum(cauchy(a,z) for a in f)
cauchy(f::AbstractVector{<:AbstractSingularity},d,z) = isempty(f) ? zero(z) :  sum(cauchy(a,d,z) for a in f)

function cauchy(f::AbstractSingularity,d,z)
    ret = zero(z)
    if contains(d,domain(f))
        ret+=cauchy(f,z)
    end
    ret
end

# Algebra

function *(a::T,b::Fun) where T<:AbstractSingularity
    T(domain(a),Δ(a)*b,Σ(a)*b)
end

function *(a::T,b::Number) where T<:AbstractSingularity
    T(domain(a),Δ(a)*b,Σ(a)*b)
end

function -(a::T) where T<:AbstractSingularity
    T(domain(a),-Δ(a),-Σ(a))
end

function +(a::Vector{T},b::S) where {T<:AbstractSingularity,S<:AbstractSingularity}
    Union{T,S}[a;[b]]
end

*(a::AbstractSingularity,b::Function) = a*Fun(b,domain(a))
*(a::Function,b::AbstractSingularity) = *(b,a)

*(a::AbstractSingularity,b::ScalarA) = a*Fun(b,domain(a))
*(a::ScalarA,b::AbstractSingularity) = *(b,a)

/(a::AbstractSingularity,b::Function) = a/Fun(b,domain(a))
/(a::AbstractSingularity,b::ScalarA) = a/Fun(b,domain(a))

function +(a::Vector{A},b::Vector{B}) where {A<:AbstractSingularity,B<:AbstractSingularity}
    Union{A,B}[a; b]
end

function +(a::A,b::B) where {A<:AbstractSingularity,B<:AbstractSingularity}
    Union{A,B}[a; b]
end

function *(a::Vector{A},b::Function,args...) where {A<:AbstractSingularity}
    [*(x,b,args...) for x in a]
end
function /(a::Vector{A},b::Function,args...) where {A<:AbstractSingularity}
    [/(x,b,args...) for x in a]
end
*(a::Function,b::Vector{A},args...) where {A<:AbstractSingularity} = *(b,a,args...)
