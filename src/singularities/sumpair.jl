export SumPair

struct SumPair{P,M} <: AbstractDefinedFunction
    p::P
    m::M
end
getindex(K::SumPair,u::Bool) = u ? K.p : K.m
evaluate(K::SumPair,u::Bool) = K[u]
evaluate(K::SumPair,α,u) = K[u](α)
evaluate(K::SumPair,α) = evaluate(K,α,true)+evaluate(K,α,false)
factorise(K::SumPair,args...) = K
