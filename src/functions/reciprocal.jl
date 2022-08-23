struct ScalarReciprocal{T} <: AbstractScalarFunction
    value::T
end
const ScalarR = ScalarReciprocal
evaluate(f::ScalarR, args...) = inv(f.value(args...))
inv(a::ScalarA) = ScalarR(a)
inv(a::ScalarR) = a.value
