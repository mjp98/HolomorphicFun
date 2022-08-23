export unitnormal, unittangent
export extendtoray, extend, parallel,intersects, truncate, isunbounded
export HalfPlane


unittangent(d::Union{Line,Ray}) = cis(angle(d))
unittangent(d::Segment) = cis(angle(d.b - d.a))
unittangent(::PeriodicPoint{T}) where T = one(T)
unittangent(::Point{T}) where T = one(T)
unitnormal(d::Domain) = im*unittangent(d)
unittangent(sp::Space) = unittangent(domain(sp))
unitnormal(sp::Space) = unitnormal(domain(sp))

isunbounded(sp::Space) = isunbounded(domain(sp))
isunbounded(d::Domain) = false
isunbounded(d::Union{Line,Ray,PeriodicLine}) = true

function extend(x::Ray{S}) where S
    Line{S}(x.center)
end
function extend(x::Segment)
    θ = angle(x.b - x.a)/pi
    return Line{θ}(x.a)
end
function extendtoray(x::Segment)
    θ = angle(x.b - x.a)/pi
    return Ray{θ}(x.a)
end

function truncate(x::Ray{S},y::Real) where S
    return Segment(x.center,x.center+y*cispi(S))
end

parallel(x::Line,y::Line) = (mod(angle(x)-angle(y),pi) == 0)
intersects(x::Line,y::Line) = !parallel(x,y) || ((x.center ∈ y) || (y.center ∈ x))

abstract type AbstractComplexDomain end

struct HalfPlane{T} <: AbstractComplexDomain
end
HalfPlane(uhp::Bool) = HalfPlane{uhp}()

function contains(::HalfPlane{T},x) where T
    if T
        return imag(x)>0 || ((imag(x)==0) && (real(x)>0))
    else
        return imag(x)<0 || ((imag(x)==0) && (real(x)<0))
    end
end

contains(d::HalfPlane,p::Point) = contains(d,p.x)
contains(d::HalfPlane,p::PeriodicPoint) = contains(d,p.x)

function contains(domain::HalfPlane{S},s::Segment{T}) where {S,T}
    # True for any convex domain
    return contains(domain,Point(s.a)) && in(Point(s.b),domain)
end

function contains(domain::HalfPlane{T},r::Ray{S}) where {T,S}
    containscenter = contains(domain,Point(r.center))
    if T
        containsangle = S >= 0
    else
        containsangle = S <= 0 || isone(S)
    end
    return containscenter && containsangle
end
