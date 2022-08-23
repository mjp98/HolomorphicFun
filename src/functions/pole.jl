struct ScalarPole{T<:Number} <: AbstractScalarFunction
    position::T
    residue::T
    order::Int
    function ScalarPole(x,r,order=1)
        xr = promote(complex(x),r)
        new{eltype(xr)}(xr...,order)
    end
end
residue(f::ScalarPole) = f.residue
position(f::ScalarPole) = f.position
order(f::ScalarPole) = f.order
evaluate(f::ScalarPole,z) = residue(f)/(position(f)-z)^order(f)

-(f::ScalarPole) = ScalarPole(position(f),-residue(f),order(f))

function Pole(f::ScalarPole)
d = PeriodicPoint(position(f))
Δ = Fun(Hardy{false}(PeriodicPoint(position(f))),[residue(f)*(2π*im)])
Σ = Fun(Hardy{true}(PeriodicPoint(position(f))),[0.])
return Pole(d,Δ,Σ)
end

function singularities(f::ScalarPole)
    d = PeriodicPoint(position(f))
    Δ = Fun(Hardy{false}(PeriodicPoint(position(f))),[residue(f)*(2π*im)])
    Σ = Fun(Hardy{true}(PeriodicPoint(position(f))),[0.0+0.0im])
    return [Pole(d,Δ,Σ)]

    # d = Point(position(f))
    # Δ = Fun(d,[residue(f)])
    # Σ = Fun(zero,d)
    # return [Pole(d,Δ,Σ)]
end
