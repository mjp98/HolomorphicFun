struct ScalarConstant{T} <: AbstractScalarFunction
    value::T
end
const ScalarC = ScalarConstant
evaluate(f::ScalarC, args...) = f.value

*(a::AbstractSingularity,b::ScalarC) = isnan(b.value) ? a*NaN : a*Fun(b,domain(a))


*(a::ScalarC,b::ScalarC) = ScalarC(a.value*b.value)

singularities(::ScalarConstant{T}) where T = Pole{T}[]
