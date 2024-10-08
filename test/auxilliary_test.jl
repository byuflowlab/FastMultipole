@testset "complex" begin
    z1 = rand(Complex{Float64})
    z2 = rand(Complex{Float64})
    z3 = rand(Complex{Float64})

    # addition
    @test real(z1+z2) ≈ FastMultipole.complex_add(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1+z2) ≈ FastMultipole.complex_add(real(z1), imag(z1), real(z2), imag(z2))[2]

    # subtraction
    @test real(z1-z2) ≈ FastMultipole.complex_subtract(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1-z2) ≈ FastMultipole.complex_subtract(real(z1), imag(z1), real(z2), imag(z2))[2]

    # multiplication
    @test real(z1*z2) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1*z2) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2))[2]

    @test real(z1*z2*z3) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2), real(z3), imag(z3))[1]
    @test imag(z1*z2*z3) ≈ FastMultipole.complex_multiply(real(z1), imag(z1), real(z2), imag(z2), real(z3), imag(z3))[2]

    # division
    @test real(z1/z2) ≈ FastMultipole.complex_divide(real(z1), imag(z1), real(z2), imag(z2))[1]
    @test imag(z1/z2) ≈ FastMultipole.complex_divide(real(z1), imag(z1), real(z2), imag(z2))[2]
    @test real(z1/z2) ≈ FastMultipole.complex_divide_real(real(z1), imag(z1), real(z2), imag(z2))
    @test imag(z1/z2) ≈ FastMultipole.complex_divide_imag(real(z1), imag(z1), real(z2), imag(z2))

	# cross product
	z1_vec = SVector{3}(rand(Complex{Float64}) for _ in 1:3)
	z2_vec = SVector{3}(rand(Complex{Float64}) for _ in 1:3)
	@test real.(cross(z1_vec,z2_vec)) ≈ FastMultipole.complex_cross_real(real(z1_vec[1]), imag(z1_vec[1]), real(z1_vec[2]), imag(z1_vec[2]), real(z1_vec[3]), imag(z1_vec[3]), real(z2_vec[1]), imag(z2_vec[1]), real(z2_vec[2]), imag(z2_vec[2]), real(z2_vec[3]), imag(z2_vec[3]))
end

@testset "cartesian to spherical" begin
    # cartesian to spherical
    rho = 1.0
    theta = pi/4
    phi = pi/2
    that = [rho, theta, phi]
    x = rho * sin(theta) * cos(phi)
    y = rho * sin(theta) * sin(phi)
    z = rho * cos(theta)
    this = [x,y,z]
    ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(this)
    this = [ρ,θ,ϕ]
    for i in 1:3
        @test isapprox(this[i], that[i]; atol=1e-10)
    end
end

