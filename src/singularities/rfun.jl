
struct RFun{U,T,S} <: AbstractDefinedFunction
    fun::T
    domain::S
end
function RFun{U}(fun::T,domain::S) where {U,T,S}
    RFun{U,T,S}(fun,domain)
end

*(a::RFun,b::ScalarA) = ScalarΠ([a,b])
*(b::ScalarA,a::RFun) = ScalarΠ([a,b])
/(a::RFun,b::ScalarA) = ScalarΠ([a,inv(b)])


evaluate(f::RFun{true},z) = regular(f.fun,z)
evaluate(f::RFun{false},z) = singular(f.fun,z)

function RFun(f::SumSplit,u)
    if (u && f.domain) || (!u && !f.domain)
        return RFun{true}(f.fun,f.domain)
    else
        return RFun{false}(f.fun,f.domain)
    end
end

RFun(f::SumSplit) = RFun(f,true),RFun(f,false)
