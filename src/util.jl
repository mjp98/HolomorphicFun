const complexfloat = complex ∘ float
const cfloat = complexfloat
candidatefloat(x) = float(x)
candidatefloat(::Taylor1{T}) where T = float(one(T))
epsreal(z) = eps(real(candidatefloat(one(z))))
