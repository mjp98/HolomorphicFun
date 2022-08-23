struct ScalarWrapper{T} <: AbstractScalarFunction
    fun::T
end

evaluate(f::ScalarWrapper,z) = f.fun(z)
