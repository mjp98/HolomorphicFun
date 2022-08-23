struct SingularFunction{T,S} <: AbstractDefinedFunction
    fun::T
    singularities::S
end

const SFun = SingularFunction

evaluate(f::SFun,z::Number) = evaluate(f.fun,z)
singularities(f::SFun) = f.singularities
function domain(f::SFun)
    sing = singularities(f)
    if isempty(sing)
        return EmptySpace()
    else
        return UnionDomain([domain(s) for s in sing]...)
    end
end

getsingularity(f::SFun,i) = getindex(singularities(f),i)
function singularityindex(f::SFun,d::Domain)
    findfirst(z->domain(z)==d,singularities(f))
end
function singularityondomain(f::SFun,d::Domain)
    getsingularity(f,singularityindex(f,d))
end

hassingularities(f::ScalarA) = false
SFun(f::ScalarA) = SFun(f,singularities(f))

singular(f::SFun,args...) = cauchy(singularities(f),args...)
regular(f::SFun,z) = evaluate(f,z) - singular(f,z)
regular(f::SFun,d,z) = evaluate(f,z) - singular(f,d,args...)

removeuhp(f::SFun,args...) = cauchy(singularities(f),args...)
regular(f::SFun,args...) = evaluate(f,args...) - singular(f,args...)

removedsingularity(f::SFun,i,z) = removedsingularity(getsingularity(f,i),z)
removedsingularityapprox(f::SFun,i,z) = evaluate(f,z) - getsingularity(f,i)(z)

function filtersingularities(f::SFun,d)
    s = filter(x->contains(d,domain(x)),singularities(f))
    return SFun(f.fun,s)
end

# Algebra

function *(a::SFun,b::Function,args...)
    SFun(a.fun*b,*(singularities(a),b,args...))
end
*(a::Function,b::SFun,args...) = *(b,a,args...)
function /(a::SFun,b::Function,args...)
    SFun(a.fun/b.fun,/(singularities(a),b.fun,args...))
end
function *(a::SFun,b::SFun,args...)
    s = *(singularities(a),b.fun,args...) + *(a.fun,singularities(b),args...)
    SFun(a.fun*b.fun,s)
end

# Misc

function truncate(f::SFun,cutoff)
    S = Vector{AbstractSingularity}(undef,0)
    for s in singularities(f)
        if domain(s) isa Ray
            push!(S,truncate(s,cutoff))
        else
            push!(S,s)
        end
    end
    return SFun(f.fun,S)
end
