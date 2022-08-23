module HolomorphicFun

using ApproxFun
using DomainSets
using SingularIntegralEquations

import Base: +,-,/,*, inv, getindex, angle, in, contains, truncate,  convert, Number, isempty, hash, ==, first, last
import LinearAlgebra: norm

import DomainSets: Point, elements, similardomain, indomain, dimension, approx_indomain, canonicaldomain, isopenset, isclosedset, mapfrom_canonical, boundary, boundingbox, infimum,supremum, interior, closure, point_in_domain, distance_to, mapped_domain, map_domain, parametric_domain, issubset1, setdiffdomain1, setdiffdomain, StaticTypes, Translation,  EmptySpace, hashrec
import ApproxFun: center, evaluate, Fun, ncoefficients, domain, checkpoints, points, fromcanonical,tocanonical, isambiguous, setdiff, reverseorientation, intersect, issubset, AnyDomain
import ApproxFunFourier: horner
import TaylorSeries: Taylor1, taylor_expand

# import DomainSets:

# import ApproxFun:


export ScalarA, AbstractScalarFunction
export ScalarΠ, ScalarProduct
export ScalarΣ, ScalarSum
export ScalarF, ScalarFraction
export ScalarR, ScalarReciprocal
export ScalarC, ScalarConstant
export ScalarL, ScalarLinear
export ScalarWrapper
export ScalarRational
export ScalarScale
export ScalarPower
export ScalarPole, residue
export ExponentialShift, Exponential, exponent, phase
export Sqrt

export AbtractSingularity, Cut, Pole

export growsabove, growsbelow

abstract type AbstractDefinedFunction <: Function end

(f::AbstractDefinedFunction)(z) = evaluate_safe(f,z)
(f::AbstractDefinedFunction)(z,u) = evaluate_safe(f,z,u)
(f::AbstractDefinedFunction)(args...) = evaluate_safe(f, args...)

abstract type AbstractScalarFunction <: AbstractDefinedFunction end
const ScalarA = AbstractScalarFunction

include("util.jl")
include("evaluate.jl")
include("geometry/PeriodicPoint.jl")
include("geometry/geometry.jl")
include("singularities/singularities.jl")
include("functions/functions.jl")

end
