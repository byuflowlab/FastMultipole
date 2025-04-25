M̃test(n,m) = sqrt(Float64(factorial(big(n+abs(m)))) * Float64(factorial(big(n-abs(m)))) / (2*n+1))

L̃test(n,m) = sqrt(1 / (Float64(factorial(big(n+abs(m)))) * Float64(factorial(big(n-abs(m)))) * (2*n+1)))

@testset "multipole powers" begin
    
M̃ = FastMultipole.M̃

i = 1
for n in 0:5
    for m in 0:n
        @test isapprox(M̃[i], M̃test(n,m), atol=1e-12)
        i += 1
    end
end

end


@testset "local powers" begin
    
L̃ = FastMultipole.L̃

i = 1
for n in 0:5
    for m in 0:n
        @test isapprox(L̃[i], L̃test(n,m), atol=1e-12)
        i += 1
    end
end

end