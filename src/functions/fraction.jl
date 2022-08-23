struct ScalarFraction{T,S} <: AbstractScalarFunction
    numerator::T
    denominator::S
end
const ScalarF = ScalarFraction

numerator(f::ScalarF) = f.numerator
denominator(f::ScalarF) = f.denominator
numerator(f::ScalarF, args...) = f.numerator(args...)
denominator(f::ScalarF, args...) = f.denominator(args...)
evaluate(f::ScalarF, args...) = numerator(f, args...) / denominator(f, args...)
inv(a::ScalarF) = ScalarF(denominator(a),numerator(a))
/(a::ScalarA, b::ScalarA) = ScalarF(a, b)
