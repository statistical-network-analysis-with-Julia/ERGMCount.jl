using ERGMCount
using Test

@testset "ERGMCount.jl" begin
    @testset "Module loading" begin
        @test @isdefined(ERGMCount)
    end

    @testset "Reference measures" begin
        @testset "PoissonReference" begin
            ref = PoissonReference()
            @test ref.lambda == 1.0
            ref2 = PoissonReference(2.0)
            @test ref2.lambda == 2.0
        end

        @testset "GeometricReference" begin
            ref = GeometricReference()
            @test ref.prob == 0.5
            ref2 = GeometricReference(0.3)
            @test ref2.prob == 0.3
            @test_throws ArgumentError GeometricReference(0.0)
            @test_throws ArgumentError GeometricReference(1.0)
        end

        @testset "BinomialReference" begin
            ref = BinomialReference(10)
            @test ref.n == 10
            @test ref.prob == 0.5
            ref2 = BinomialReference(5, 0.3)
            @test ref2.n == 5
            @test ref2.prob == 0.3
            @test_throws ArgumentError BinomialReference(0, 0.5)
        end

        @testset "DiscUnifReference" begin
            ref = DiscUnifReference(5)
            @test ref.max == 5
            @test_throws ArgumentError DiscUnifReference(-1)
        end

        @testset "DiscUnif2Reference" begin
            ref = DiscUnif2Reference(1, 5)
            @test ref.a == 1
            @test ref.b == 5
            @test_throws ArgumentError DiscUnif2Reference(5, 1)
        end
    end

    @testset "Count-specific terms" begin
        @test SumTerm() isa SumTerm
        @test NonzeroTerm() isa NonzeroTerm
        @test GreaterthannTerm(3) isa GreaterthannTerm
        @test CountAtleastnTerm(2) isa CountAtleastnTerm
        @test CountMutualTerm() isa CountMutualTerm
        @test TransitiveTiesTerm() isa TransitiveTiesTerm
        @test CyclicalTiesTerm() isa CyclicalTiesTerm
        @test NodeOSumTerm() isa NodeOSumTerm
        @test NodeISumTerm() isa NodeISumTerm
        @test NodeSumTerm() isa NodeSumTerm
    end

    @testset "Estimation API" begin
        @test ergm_count === fit_count_ergm
    end

    @testset "Simulation API" begin
        @test isdefined(ERGMCount, :simulate_count_ergm)
    end
end
