import SingularIntegralEquations:stieltjes,setcanonicaldomain,tocanonical,hilbert

#! TYPE PIRACY
function stieltjes(S::Space{<:Ray},f,z::T) where T
    stieltjes(setcanonicaldomain(S),f,tocanonical(S,z)) +pi*hilbert(setcanonicaldomain(S),f,1)
end

function hilbert(S::Space{<:Ray},f,z::T) where T
    x = tocanonical(S,z)
    s = setcanonicaldomain(S)
    if abs(imag(x))<epsreal(T)
        x = real(x)
    end
    hilbert(s,f,x)-hilbert(s,f,1)
end

stieltjes(S::PointSpace,f,z) = first(f)/(first(domain(S).domains).x - z)
hilbert(::PointSpace,f,z) = zero(z)

stieltjes(S::Hardy{false}{<:PeriodicPoint},f,z) = -first(f)/(domain(S).x - z)
hilbert(::Hardy{false}{<:PeriodicPoint},f,z) = zero(z)
hilbert(::Hardy{true}{<:PeriodicPoint},f,z) = zero(z)
