struct ScalarRational{T} <: AbstractScalarFunction
    poles::T
    roots::T
end

roots(K::ScalarRational) = K.roots
poles(K::ScalarRational) = K.poles

nroots(K::ScalarRational) = length(roots(K))
npoles(K::ScalarRational) = length(poles(K))

inv(K::ScalarRational) = ScalarRational(roots(K),poles(K))

function evaluate(K::ScalarRational{T},z) where T
    if nroots(K)==0
        num = one(eltype(T))
    else
        num = prod(z-x for x in roots(K))
    end
    if npoles(K)==0
        den = one(eltype(T))
    else
        den = prod(z-x for x in poles(K))
    end
    return num/den
end

function residues(K::ScalarRational{T}) where T
    r = Vector{T}(undef,npoles(K))
    for (i,p) in enumerate(poles(K))
        # numerator
        if nroots(K)>0
            r[i] = prod(p-x for x in roots(K))
        else
            r[i] = one(T)
        end
        # denominator
        den = one(T)
        for (j,q) in enumerate(poles(K))
            if !(p ≈ q) # i!==j depending...
                den * (p-q)
            end
        end
        r[i] /= den
    end
    return r
end
