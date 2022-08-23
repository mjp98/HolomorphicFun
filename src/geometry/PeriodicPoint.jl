# From DomainSets
# using DomainSets
# using ApproxFun

# import LinearAlgebra: norm
# import Base: convert, Number, isempty, hash, ==,+,-, getindex, first, last
# import DomainSets: similardomain, indomain, dimension, approx_indomain, canonicaldomain, isopenset, isclosedset, mapfrom_canonical, boundary, boundingbox, infimum,supremum, interior, closure, point_in_domain, distance_to, mapped_domain, map_domain, parametric_domain, issubset1, setdiffdomain1, setdiffdomain, StaticTypes, Translation,  EmptySpace, hashrec

# import ApproxFun: checkpoints, points, fromcanonical,tocanonical, isambiguous, setdiff, reverseorientation, intersect, issubset, AnyDomain

export PeriodicPoint

"""
PeriodicPoint(x)
represents a single point at `x`.
"""
struct PeriodicPoint{T} <: PeriodicDomain{T}
x::T
end

similardomain(d::PeriodicPoint, ::Type{T}) where {T} = PeriodicPoint{T}(d.x)

convert(::Type{Number}, d::PeriodicPoint{<:Number}) = d.x
convert(::Type{N}, d::PeriodicPoint{<:Number}) where N<:Number = convert(N, convert(Number, d.x))
Number(d::PeriodicPoint) = convert(Number, d)

==(d1::PeriodicPoint,d2::PeriodicPoint) = d1.x == d2.x
hash(d::PeriodicPoint, h::UInt) = hashrec("PeriodicPoint", d.x, h)

indomain(x, d::PeriodicPoint) = x == d.x
isempty(::PeriodicPoint) = false

approx_indomain(x, d::PeriodicPoint, tolerance) = norm(x-d.x) <= tolerance

dimension(d::PeriodicPoint{Vector{T}}) where {T} = length(d.x)

canonicaldomain(d::PeriodicPoint{T}) where {T<:StaticTypes} = PeriodicPoint(zero(T))
canonicaldomain(d::PeriodicPoint{T}) where {T<:AbstractVector} =
PeriodicPoint(zeros(eltype(T),dimension(d)))

mapfrom_canonical(d::PeriodicPoint) = Translation(d.x)

isopenset(d::PeriodicPoint) = false
isclosedset(d::PeriodicPoint) = true

boundary(d::PeriodicPoint) = d
boundingbox(d::PeriodicPoint) = d.x..d.x

infimum(d::PeriodicPoint) = d.x
supremum(d::PeriodicPoint) = d.x

interior(d::PeriodicPoint{T}) where {T} = EmptySpace{T}()
closure(d::PeriodicPoint) = d

point_in_domain(d::PeriodicPoint) = d.x

distance_to(d::PeriodicPoint, x) = norm(x-d.x)

mapped_domain(invmap, p::PeriodicPoint) = PeriodicPoint(inverse(invmap, p.x))
map_domain(map, p::PeriodicPoint) = PeriodicPoint(applymap(map, p.x))
parametric_domain(map, p::PeriodicPoint) = PeriodicPoint(applymap(map, p.x))

for op in (:+,:-)
@eval $op(a::PeriodicPoint, b::PeriodicPoint) = PeriodicPoint($op(a.x,b.x))
end

# Interval minus a point:
setdiffdomain(d::Interval, p::PeriodicPoint) = setdiffdomain(promote_domains((d,p))...)
function setdiffdomain(d::Interval{L,R,T}, p::PeriodicPoint{T}) where {L,R,T}
a = leftendpoint(d)
b = rightendpoint(d)
x = p.x

a == x && return Interval{:open,R,T}(a,b)
a < x < b && return UnionDomain(Interval{L,:open,T}(a,p.x), Interval{:open,R,T}(p.x,b))
b == x && return Interval{L,:open,T}(a,b)
return d
end

issubset1(d1::PeriodicPoint, d2) = d1.x ∈ d2

setdiffdomain1(p::PeriodicPoint, d2) = issubset(p, d2) ? EmptySpace{eltype(p)}() : p

intersectdomain(d1::PeriodicPoint, d2::PeriodicPoint) = d1.x ∈ d2 ? d1 : EmptySpace{eltype(d1)}()

# From ApproxFun

PeriodicPoint(::AnyDomain) = PeriodicPoint(NaN)
PeriodicPoint{T}(x::AnyDomain) where {T} = new(T(NaN))

convert(::Type{PeriodicPoint},::AnyDomain) = PeriodicPoint(NaN)
convert(::Type{PeriodicPoint{T}},::AnyDomain) where T = PeriodicPoint{T}(NaN)


isambiguous(d::PeriodicPoint) = isnan(d.x)


norm(p::PeriodicPoint) = norm(p.x)

getindex(p::PeriodicPoint,k...) = p.x[k...]
first(p::PeriodicPoint) = p.x
last(p::PeriodicPoint) = p.x

Base.isnan(d::PeriodicPoint) = false


issubset(a::PeriodicPoint,d::UnionDomain) = a.x in d

intersect(a::PeriodicPoint,b::PeriodicPoint) = b.x in a ? b : EmptyDomain()
intersect(a::UnionDomain,b::PeriodicPoint) = b.x in a ? b : EmptyDomain()
intersect(a::Domain,b::PeriodicPoint) = b.x in a ? b : EmptyDomain()
intersect(b::PeriodicPoint,a::UnionDomain) = b.x in a ? b : EmptyDomain()
intersect(b::PeriodicPoint,a::Domain) = b.x in a ? b : EmptyDomain()

setdiff(a::PeriodicPoint,b::PeriodicPoint) = a==b ? EmptyDomain() : a
reverseorientation(a::PeriodicPoint) = a


canonicaldomain(a::PeriodicPoint) = PeriodicPoint(0.)
tocanonical(a::PeriodicPoint,x) = x-a.x
fromcanonical(a::PeriodicPoint,x) = x+a.x

points(a::PeriodicPoint,n::Integer) = eltype(a)[a.x]
checkpoints(a::PeriodicPoint) = eltype(a)[a.x]




Fun(f::Function,sp::Taylor{<:PeriodicPoint}) = Fun(sp,complex.(taylor_expand(f,domain(sp).x).coeffs))

Fun(f::Function,d::PeriodicPoint) = Fun(f,Taylor(d))

function evaluate(f::AbstractVector,S::Hardy{false,D},z) where D<:PeriodicPoint
    z = inv(z-domain(S).x)
    z*horner(f,z)
end

function evaluate(f::AbstractVector,S::Hardy{true,D},z) where D<:PeriodicPoint
    horner(f,z-domain(S).x)
end
