struct Exponential{T} <:  AbstractScalarFunction
    exponent::T
    phase::T
end
Exponential(exponent::T,phase=zero(T)) where T = Exponential(promote(exponent,phase)...)

exponent(x::Exponential) = x.exponent
phase(x::Exponential) = x.phase

evaluate(x::Exponential,z) = exp(im*(exponent(x)*z-phase(x)))

isanalytic(f::Exponential) = true
growsbelow(f::Exponential) = real(exponent(f))>0

*(a::Exponential,b::Exponential) = Exponential(exponent(a)+exponent(b),phase(a)+phase(b))
/(a::Exponential,b::Exponential) = Exponential(exponent(a)-exponent(b),phase(a)-phase(b))

struct ExponentialShift{T<:Function,S<:Number} <: AbstractScalarFunction
    fun::T
    exponent::S
    #origin::S
    phase::S
    function ExponentialShift(f::T, e, p) where {T}
        args = promote(e,  p)
        return new{T,eltype(args)}(f, args...)
    end
end

ExponentialShift(f::ExponentialShift,e,p) = ExponentialShift(f.fun,exponent(f)+e,phase(f)+p)

fun(f::ExponentialShift) = f.fun
fun(f::ExponentialShift, z) = f.fun(z)
exponent(f::ExponentialShift) = f.exponent
phase(f::ExponentialShift) = f.phase
#origin(f::ExponentialShift) = f.origin

evaluate(f::ExponentialShift,z) = fun(f, z) * exp(im * exponent(f) * z - im * phase(f))

-(f::ExponentialShift) = ExponentialShift(-f.fun,exponent(f),phase(f))
growsbelow(f::ExponentialShift) = real(exponent(f))>0
function growsindomain(f::Union{ExponentialShift,Exponential},d)
    # @show contains(HalfPlane(true),d)
    # @show growsabove(f)
    # @show contains(HalfPlane(false),d)
    # @show growsbelow(f)
    (contains(HalfPlane(true),d) && growsabove(f)) || (contains(HalfPlane(false),d) && growsbelow(f))
end
Exponential(a::ExponentialShift) = Exponential(exponent(a),phase(a))

*(a::Exponential,b::ExponentialShift) = ExponentialShift(b.fun,exponent(a)+exponent(b),phase(a)+phase(b))
*(a::ExponentialShift,b::Exponential) = *(b,a)
*(a::ExponentialShift,b::ScalarA) = ExponentialShift(a.fun*b,exponent(a),phase(a))
*(a::ScalarA,b::ExponentialShift) = *(b,a)

*(a::Exponential,b::ScalarA) = ExponentialShift(b,exponent(a),phase(a))
*(a::Exponential,b::ScalarΣ) = ScalarΣ([ExponentialShift(c,exponent(a),phase(a)) for c in b.elements])

function *(b::Cut,a::Union{Exponential,ExponentialShift})
    d = domain(b)
    if isunbounded(d) && growsindomain(a,d)
        @error "singularity non-integrable"
        return b*NaN
    else
        return b*Fun(a,domain(b))
    end
end
*(a::Union{Exponential,ExponentialShift},b::Cut) = *(b,a)
singularities(a::ExponentialShift) = Exponential(a)*singularities(a.fun)
sumsplit(f::SFun{<:ExponentialShift}) = sumsplit(f,!growsabove(f))

growsbelow(f::ScalarA) = false
growsabove(f::ScalarA) = !growsbelow(f)

growsbelow(f::SFun) = growsbelow(f.fun)
growsabove(f::SFun) = growsabove(f.fun)
