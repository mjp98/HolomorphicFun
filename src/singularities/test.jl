function testremoved(f::SFun,N=10,ε=1e-8)
    all(testremoved.(f,singularities(f),N,ε))
end

function testremoved(f::SFun,c::AbstractSingularity,N=10,ε=1e-8)
    x = points(domain(c),N)
    n = unitnormal(domain(c))
    # @show evaluate.(f,x .+ε*n) - c.(x .+ε*n)
    # @show removedsingularity.(c,x)
    isapprox(evaluate.(f,x .+ε*n) - c.(x .+ε*n),removedsingularity.(c,x);atol=10ε ,rtol = 10ε)
end
