@testset "solid harmonics" begin

p = 7
expansion_order = Val(p)
harmonics = zeros(2,2, FastMultipole.harmonic_index(p,p))
ρ = 1.7
θ = 3π/8
ϕ = 5π/3

# test regular harmonics
FastMultipole.regular_harmonics!(harmonics, ρ, θ, ϕ, expansion_order)
for n in 0:p
    for m in 0:n
        Rnm = (-1)^n * im^(m*1.0) * ρ^(n*1.0) * Plm(cos(θ),n,m) * exp(im*m*ϕ) / factorial(n+m)
        @test isapprox(real(Rnm), harmonics[1,1,FastMultipole.harmonic_index(n,m)]; atol=1e-12)
        @test isapprox(imag(Rnm), harmonics[2,1,FastMultipole.harmonic_index(n,m)]; atol=1e-12)
    end
end

# test irregular harmonics
FastMultipole.irregular_harmonics!(harmonics, ρ, θ, ϕ, expansion_order)
for n in 0:p
    for m in 0:n
        Snm = (-1)^m * im^(m*1.0) * ρ^((-n-1)*1.0) * Plm(cos(θ),n,m) * exp(im*m*ϕ) * factorial(n-m)
        @test isapprox(real(Snm), harmonics[1,1,FastMultipole.harmonic_index(n,m)]; atol=1e-12)
        @test isapprox(imag(Snm), harmonics[2,1,FastMultipole.harmonic_index(n,m)]; atol=1e-12)
    end
end

end
