import ApproxFun: Fun, domain
import SingularIntegralEquations: cauchy
import Base: *, /,-

export AbstractSingularity, Cut, Pole

export SFun, RFun
export regular, singular
export singularities, singularityondomain, filtersingularities
export removedsingularity, removedsingularityapprox

export testremoved


export SumSplit, sumsplit

include("sumpair.jl")
include("cauchy.jl")
include("abstract.jl")
include("cut.jl")
include("pole.jl")
include("sfun.jl")
include("sumsplit.jl")
include("rfun.jl")
include("test.jl")
