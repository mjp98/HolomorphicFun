
struct SumSplit{T,S} <: AbstractDefinedFunction
    fun::T
    domain::S
end
evaluate(f::SumSplit,z) = evaluate(f.fun,z)

function evaluate(f::SumSplit,z,u)
    if (u && f.domain) || (!u && !f.domain)
        return regular(f.fun,z)
    else
        return singular(f.fun,z)
    end
end

sumsplit(f::SFun,u::Bool) = SumSplit(filtersingularities(f,HalfPlane(u)),u)
