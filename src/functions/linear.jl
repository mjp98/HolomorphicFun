struct ScalarLinear{T} <: AbstractScalarFunction
    a::T
    b::T
    function ScalarLinear(a,b)
        args = promote(a,b)
        new{eltype(args)}(args...)
    end
end
const ScalarL = ScalarLinear
root(f::ScalarL) = -f.b/f.a
evaluate(f::ScalarLinear, z, args...) = f.a*z + f.b
