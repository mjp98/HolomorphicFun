function evaluate_close(f, z::T; ε=100epsreal(T)) where {T}
    evaluate(f, z + ε)
end
function evaluate_close(f, z::T, args...; ε=100epsreal(T)) where {T}
    evaluate(f, z + ε, args...)
end
function evaluate_safe(f, z::T; ε=100epsreal(T)) where {T}
    small = 100ε
    exact = evaluate_close(f, z; ε=0)
    if isnan(exact)
        close = evaluate_close(f, z; ε)
        absval = abs(close)
        if absval < small
            close=0*close
        end
        return close
    else
        return exact
    end
end


function evaluate_safe(f, z::T, args...; ε=100epsreal(T)) where {T}
    small = 100ε
    large = inv(small)
    exact = evaluate_close(f, z, args...; ε=0)
    if isnan(exact)
        close = evaluate_close(f, z, args...; ε)
        absval = abs(close)
        absval < small && return zero(close)
        absval > large && return NaN * close
        return close
    else
        return exact
    end
end
