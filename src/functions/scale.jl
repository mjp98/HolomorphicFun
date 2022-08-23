struct ScalarScale{S,T<:Number} <: AbstractScalarFunction
    fun::S
    scale::T
end
fun(f::ScalarScale) = f.fun
fun(f::ScalarScale,z...) = f.fun(z...)
scale(f::ScalarScale) = f.scale

evaluate(f::ScalarScale,z) = scale(f)*fun(f, z)

-(a::AbstractScalarFunction) = ScalarScale(a,-1)
