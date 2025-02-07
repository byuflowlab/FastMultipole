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

@testset "minimum_edge_distance" begin
    center1 = [0.0, 0.0, 0.0]
    box1 = [1.0, 1.0, 1.0]
    center2 = [2.0, 2.0, 2.0]
    box2 = [1.0, 1.0, 1.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 0.0; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    center2 = [3.0, 3.0, 3.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 1.0; atol=1e-10)
    @test isapprox(dy, 1.0; atol=1e-10)
    @test isapprox(dz, 1.0; atol=1e-10)

    center2 = [-2.0, -2.0, -2.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 0.0; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    center2 = [-3.0, -3.0, -3.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, -1.0; atol=1e-10)
    @test isapprox(dy, -1.0; atol=1e-10)
    @test isapprox(dz, -1.0; atol=1e-10)

    center2 = [2.0, 0.0, 0.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 0.0; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    center2 = [1.0, 3.0, 0.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 1.0; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    center2 = [-1.0, 3.0, 0.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 1.0; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    center2 = [1.0, 3.0, 0.0]
    box2 = [0.1, 0.1, 0.1]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 2.0-0.1; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    center2 = [0.0, 3.0, 0.0]
    box2 = [0.1, 0.1, 0.1]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 2.0-0.1; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)

    box2 = [5.0, 1.0, 1.0]
    dx, dy, dz = FastMultipole.minimum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 0.0; atol=1e-10)
    @test isapprox(dy, 1.0; atol=1e-10)
    @test isapprox(dz, 0.0; atol=1e-10)
end

@testset "maximum_edge_distance" begin
    center1 = [0.0, 0.0, 0.0]
    box1 = [1.0, 1.0, 1.0]
    center2 = [2.0, 2.0, 2.0]
    box2 = [1.0, 1.0, 1.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 4.0; atol=1e-10)
    @test isapprox(dy, 4.0; atol=1e-10)
    @test isapprox(dz, 4.0; atol=1e-10)

    center2 = [3.0, 3.0, 3.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 5.0; atol=1e-10)
    @test isapprox(dy, 5.0; atol=1e-10)
    @test isapprox(dz, 5.0; atol=1e-10)

    center2 = [-2.0, -2.0, -2.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, -4.0; atol=1e-10)
    @test isapprox(dy, -4.0; atol=1e-10)
    @test isapprox(dz, -4.0; atol=1e-10)

    center2 = [-3.0, -3.0, -3.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, -5.0; atol=1e-10)
    @test isapprox(dy, -5.0; atol=1e-10)
    @test isapprox(dz, -5.0; atol=1e-10)

    center2 = [2.0, 0.0, 0.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 4.0; atol=1e-10)
    @test isapprox(dy, 2.0; atol=1e-10)
    @test isapprox(dz, 2.0; atol=1e-10)

    center2 = [1.0, 3.0, 0.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 3.0; atol=1e-10)
    @test isapprox(dy, 5.0; atol=1e-10)
    @test isapprox(dz, 2.0; atol=1e-10)

    center2 = [-1.0, 3.0, 0.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, -3.0; atol=1e-10)
    @test isapprox(dy, 5.0; atol=1e-10)
    @test isapprox(dz, 2.0; atol=1e-10)

    center2 = [1.0, 3.0, 0.0]
    box2 = [0.1, 0.1, 0.1]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 2.1; atol=1e-10)
    @test isapprox(dy, 4.1; atol=1e-10)
    @test isapprox(dz, 1.1; atol=1e-10)

    center2 = [0.0, 3.0, 0.0]
    box2 = [0.1, 0.1, 0.1]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 1.1; atol=1e-10)
    @test isapprox(dy, 4.1; atol=1e-10)
    @test isapprox(dz, 1.1; atol=1e-10)

    box2 = [5.0, 1.0, 1.0]
    dx, dy, dz = FastMultipole.maximum_edge_distance(center1, box1, center2, box2)
    @test isapprox(dx, 6.0; atol=1e-10)
    @test isapprox(dy, 5.0; atol=1e-10)
    @test isapprox(dz, 2.0; atol=1e-10)
end