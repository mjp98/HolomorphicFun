
abstract type AbstractCompositeScalarFunction <: AbstractScalarFunction end
elements(f::AbstractCompositeScalarFunction) = f.elements
getindex(f::AbstractCompositeScalarFunction, i) = getindex(elements(f), i)

eltype(f::AbstractCompositeScalarFunction) = eltype(elements(f))
ncoefficients(f::AbstractScalarFunction) = 0
ncoefficients(f::AbstractCompositeScalarFunction) = sum(ncoefficients(x) for x in f)

import Base: iterate,//, eltype
iterate(f::AbstractCompositeScalarFunction,args...) = iterate(elements(f),args...)

struct ScalarSum{T} <: AbstractCompositeScalarFunction
    elements::Vector{T}
end
struct ScalarProduct{T} <: AbstractCompositeScalarFunction
    elements::Vector{T}
end
const ScalarΠ = ScalarProduct
const ScalarΣ = ScalarSum
evaluate(f::ScalarΣ, z::Number) = sum(evaluate(x,z) for x in f.elements)
evaluate(f::ScalarΠ, z::Number) = prod(evaluate(x,z) for x in f.elements)
evaluate(f::ScalarΣ, z,args...) = sum(evaluate(x,z,args...) for x in elements(f))
evaluate(f::ScalarΠ, z,args...) = prod(evaluate(x,z,args...) for x in elements(f))

# evaluate(f::ScalarΣ, z::Number) = sum(x(z) for x in f.elements)
# evaluate(f::ScalarΠ, z::Number) = prod(x(z) for x in f.elements)
# evaluate(f::ScalarΣ, z,args...) = sum(x(z,args...) for x in f.elements)
# evaluate(f::ScalarΠ, z,args...) = prod(x(z,args...) for x in f.elements)

+(a::A, b::B) where {A<:ScalarA,B<:ScalarA} = ScalarΣ(Union{A,B}[a,b])
-(a::A, b::B) where {A<:ScalarA,B<:ScalarA} = ScalarΣ(Union{A,B}[a,-b])

-(a::A, b::B) where {A<:ScalarΣ,B<:Vector{<:AbstractSingularity}} = ScalarΣ(Union{A,B}[a,-b])


+(a::ScalarΣ, b::ScalarA) = ScalarΣ{Union{eltype(a),typeof(b)}}(vcat(elements(a), b))
+(a::ScalarA, b::ScalarΣ) = +(b, a)

+(a::ScalarΣ, b::AbstractSingularity) = ScalarΣ{Union{eltype(a),typeof(b)}}(vcat(elements(a), b))
+(a::AbstractSingularity, b::ScalarΣ) = +(b, a)

+(a::ScalarΣ, b::Vector{<:AbstractSingularity}) = ScalarΣ{Union{eltype(a),typeof(b)}}(vcat(elements(a), b))
+( b::Vector{<:AbstractSingularity},a::ScalarΣ) = ScalarΣ{Union{eltype(a),typeof(b)}}(vcat(elements(a), b))


*(a::ScalarΣ, b::ScalarA) = ScalarΣ([x*b for x in a.elements])
/(a::ScalarΣ, b::ScalarA) = ScalarΣ([x*b for x in a.elements])
//(a::ScalarA, b::ScalarA) = ScalarFraction(a,b)
Π(a::ScalarA, b::ScalarA) = ScalarΠ([a,b])
+(a::ScalarΣ, b::ScalarΣ) = ScalarΣ{Union{eltype(a),eltype(b)}}(vcat(elements(a), elements(b)))

*(a::A, b::B) where {A<:ScalarA,B<:ScalarA} = ScalarΠ{Union{A,B}}(Union{A,B}[a,b])

*(a::ScalarΠ, b::ScalarA) = ScalarΠ{Union{eltype(a),typeof(b)}}(vcat(elements(a), b))
*(a::ScalarA, b::ScalarΠ) = *(b, a)
*(a::ScalarΠ, b::ScalarΠ) = ScalarΠ{Union{eltype(a),eltype(b)}}(vcat(elements(a), elements(b)))


SFun(f::ScalarΣ) = ScalarΣ(SFun.(f.elements))

*(a::ScalarΣ{<:SFun},b::SFun) = ScalarΣ([x*b for x in a.elements])
*(a::SFun,b::ScalarΣ{<:SFun}) = *(b,a)
sumsplit(f::ScalarΣ,args...) = ScalarΣ([sumsplit(x,args...) for x in f.elements])


ScalarScale(a::ScalarΣ,b::Number) = ScalarΣ([ScalarScale(x,b) for x in elements(a)])


-(a::A, b::B) where {A<:ScalarΣ,B<:ScalarΣ} = +(a,-b)
-(a::A) where {A<:ScalarΣ} = ScalarΣ([ScalarScale(x,-1) for x in elements(a)])


export simplifysumsplit

function simplifysumsplit(f::ScalarSum{<:SumSplit},u)
    ScalarSum([RFun(x,u) for x in elements(f)])
end

function simplifysumsplit(f::ScalarSum{<:SumSplit})
    SumPair(simplifysumsplit(f,true), simplifysumsplit(f,false))
end
