using HolomorphicFun
using ApproxFun
using DomainSets
using Test

@testset "HolomorphicFun.jl" begin
    @testset "sum" begin
        f = ScalarΣ([ScalarWrapper(sin),ScalarWrapper(cos)])
        z = randn(ComplexF64)
        @test f(z) == sin(z) + cos(z)

        h = ScalarΣ([ScalarWrapper(sin)])
        g = ScalarΣ([ScalarWrapper(cos)])

        hg = h+g
        @test hg(z) == h(z) + g(z)
    end
    @testset "constant" begin
        c = randn(ComplexF64)
        s = ScalarC(c)
        z = randn(ComplexF64)
        @test s(z) == c
    end

    @testset "reciprocal" begin
        c = randn(ComplexF64)
        s = ScalarR(ScalarC(c))
        z = randn(ComplexF64)
        @test s(z) == inv(c)
    end

    @testset "product" begin
        f = ScalarΠ([ScalarWrapper(sin),ScalarWrapper(cos)])
        z = randn(ComplexF64)
        @test f(z) == sin(z)*cos(z)

        h = ScalarΠ([ScalarWrapper(sin)])
        g = ScalarΠ([ScalarWrapper(cos)])

        hg = h*g
        @test hg(z) == h(z)*g(z)

        s = ScalarC(1)
        sh = s*h
        @test sh(z) == s(z)*h(z)
    end

    @testset "fraction" begin
        f = ScalarF(ScalarWrapper(sin),ScalarWrapper(cos))
        z = randn(ComplexF64)
        @test f(z) == sin(z)/cos(z)

        h = ScalarΠ([ScalarWrapper(sin)])
        g = ScalarΠ([ScalarWrapper(cos)])

        hg = h/g
        @test hg(z) == h(z)/g(z)

        hi = inv(hg)
        @test hi(z) == g(z)/h(z)
    end

    @testset "exponential" begin
    end

    @testset "pole" begin
    end

    @testset "power" begin
    end

    function testsing(f::SFun,args...)
        sing = singularities(f)
        for s in sing
            testsing(f,s,args...)
        end
    end

    function testsing(f::SFun,s::AbstractSingularity,args...)
        testΔ(f,s,args...)
        testΣ(f,s,args...)
    end

    function testΔ(f,s::Cut,N=10,ε = 1e-8)
        x = points(domain(s),N)
        n = unitnormal(domain(s))
        @test s.Δ.(x) ≈ f.(x .+ε*n) - f.(x .-ε*n)  rtol = 100*N*ε  atol = 100*N*ε
    end
    function testΣ(f,s::Cut,N=10,ε = 1e-8)
        x = points(domain(s),N)
        n = unitnormal(domain(s))
        @test s.Σ.(x) ≈ f.(x .+ε*n) + f.(x .-ε*n) rtol = 100*N*ε atol = 100*N*maximum(abs.(f.(x .+ε*n)))*ε
    end

    function testΔ(f,s::Pole,N=10,ε = 1e-8)
        x = points(domain(s),N)[1]
        n = unitnormal(domain(s))
        @test s.Δ.coefficients[1] ≈ (-n*ε*(f(x +ε*n) - f(x -ε*n))/2)*(2π*im)  rtol = 100*N*ε  atol = 100*N*ε
    end

    function testΣ(f,s::Pole,N=10,ε = 1e-8)
        x = points(domain(s),N)
        n = unitnormal(domain(s))
        @test s.Σ.(x) ≈ (f.(x .+ε*n) + f.(x .-ε*n))/2 rtol = 100*N*ε atol = 100*N*maximum(abs.(f.(x .+ε*n)))*ε
    end

    @testset "sqrt" begin


        f = Sqrt{true}(1.0+0im,0.)
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))

        f = Sqrt{true}(1.0+0im,0.5)
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))

        f = Sqrt{true}(1.0+0im,-0.5)
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))


        f = Sqrt{true}(randn(ComplexF64),randn(Float64))
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))
        @testset  "oversqrt" begin

        f = Sqrt{false}(1.0+0im,0.)
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))

        f = Sqrt{false}(1.0+0im,0.5)
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))

        f = Sqrt{false}(1.0+0im,-0.5)
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))
        f = Sqrt{false}(randn(ComplexF64),randn(Float64))
        sing = singularities(f)

        testΔ(f,first(sing))
        testΣ(f,first(sing))
        end
    end


    @testset "singularities" begin
        a = SFun(Sqrt{false}(1.0+0im,-0.5))
        b = SFun(ScalarPole(-1.0,2))
        c = a*b
        testsing(c)
        @test testremoved(c)
    end

    @testset "geometry" begin
        @testset  "extend" begin
            z = randn(ComplexF64)
            @test extend(Ray{0.5}(z)) == Line{0.5, ComplexF64}(z)
            @test extend(Segment(z,z+1im)) == Line{0.5, ComplexF64}(z)
            @test extendtoray(Segment(z,z+1im)) == Ray{0.5, ComplexF64}(z)
        end
        @testset "parallel" begin
            z = randn(ComplexF64)
            @test parallel(Line{0.5}(0.1),Line{0.5}(z))
            @test parallel(Line{0.5}(0.1),Line{-0.5}(z))
            @test !parallel(Line{0.5}(0.1),Line{0.49}(z))
        end
        @testset "intersects" begin
            z = randn(ComplexF64)
            @test intersects(Line{0.5}(z),Line{0.5}(z+1im))
            @test intersects(Line{0.5}(z),Line{0.2}(z+1im))
            @test !intersects(Line{0.5}(z),Line{0.5}(z+1))
        end

        @testset "contains" begin
            z = randn(ComplexF64)
            uhp = HalfPlane(true)
            @test contains(uhp,Point(1im))
            lhp = HalfPlane(false)
            @test contains(lhp,Point(-1im))
        end
    end

    @testset "Exponential Sqrt" begin
        e  = Exponential(1.)
        s =  Sqrt{true}(1.,0.5)
        f = e*SFun(s)
        testsing(f)
        @test testremoved(f)

        e  = Exponential(-2.)
        s =  Sqrt{false}(4.,-0.5)
        f = e*SFun(s)
        testsing(f)
        @test testremoved(f)

    end

    @testset "filter" begin
        e  = Exponential(1.)
        s =  Sqrt{true}(1.,0.5)
        f = e*SFun(s)
        @test length(singularities(filtersingularities(f,HalfPlane(true)))) == 1
        a = SFun(Sqrt{false}(1.0+0im,0.5))
        b = SFun(ScalarPole(-1.0,2))
        c = a*b
        @test length(singularities(filtersingularities(c,HalfPlane(true)))) == 1
    end


end
