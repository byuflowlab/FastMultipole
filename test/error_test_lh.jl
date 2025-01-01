using Random
using LegendrePolynomials
using Statistics

include("../test/vortex.jl")
include("../test/evaluate_multipole.jl")

function flatten_derivatives!(jacobian)
    # velocity
    velocity = zeros(3)
    velocity[1] = jacobian[2,3] - jacobian[3,2]
    velocity[2] = jacobian[3,1] - jacobian[1,3]
    velocity[3] = jacobian[1,2] - jacobian[2,1]

    return velocity
end

function Snm(r,θ,ϕ,n,m)
    if abs(m) > n
        return 0.0 + 0im
    end
    t1 = (1.0im)^(-abs(m)) * factorial(big(n-abs(m))) / r^(n+1)
    t2 = Plm(cos(θ),n,abs(m)) * exp(im*m*ϕ)
    return Float64(t1*t2)
end

function Rnm(ρ,θ,ϕ,n,m)
    if abs(m) > n
        return 0.0 + 0im
    end
    t1 = (-1)^n * (1.0*im)^abs(m) / Float64(factorial(big(n+abs(m))))
    t2 = ρ^n
    t3 = Plm(cos(θ),n,abs(m))
    t4 = exp(im*m*ϕ)
    return t1 * t2 * t3 * t4
end

function rotate(r̂)
    ẑ = r̂
    x̂ = SVector{3}(1.0,0,0) - dot(SVector{3}(1.0,0,0), r̂) * r̂
    if iszero(x̂)
        x̂ = SVector{3}(0.0,1.0,0.0) - dot(SVector{3}(0.0,1,0), r̂) * r̂
    end
    x̂ /= norm(x̂)
    ŷ = cross(ẑ, x̂)
    R = vcat(x̂', ŷ', ẑ')
    return x̂, ŷ, ẑ, R
end

function preintegrate_multipole_χ(r̂, n_max; s=1.0)
    # rotate frame
    _, _, _, R = rotate(r̂)

    # preallocate results
    ITns = zeros(Complex{Float64},2,n_max)

    # perform integration
    nx = 100
    dx = dy = dz = s / nx
    x0 = -s*0.5 + dx * 0.5
    y0 = -s*0.5 + dy * 0.5
    z0 = -s*0.5 + dz * 0.5
    for x in range(x0, step=dx, length=nx)
        for y in range(y0, step=dy, length=nx)
            for z in range(z0, step=dz, length=nx)
                ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(R * SVector{3}(x,y,z))
                sθ, cθ = sincos(θ)
                e2imϕ = exp(-2*im*ϕ)

                # Legendre polynomial recursions
                P_nm1_0 = 1.0

                # n=1
                ITns[1,1] += -im * abs(P_nm1_0)
                ITns[2,1] += -abs(P_nm1_0)

                # recurse
                P_nm1_0 = cθ
                P_nm2_0 = 1.0

                # n=2
                ITns[1,2] += -1.5 * ρ * im * abs(P_nm1_0)
                ITns[2,2] += -1.5 * ρ * abs(P_nm1_0)

                # next polynomial of order 0
                n = 2
                P_n_0 = ( (2*n-1) * cθ * P_nm1_0 - (n-1) * P_nm2_0 ) / n

                # recurse (n=3)
                P_nm2_0 = P_nm1_0
                P_nm1_0 = P_n_0

                # next polynomial of order 1 (n=2)
                P_nm1_1 = -sθ
                P_n_1 = (2*n-1) * cθ * P_nm1_1 / (n-1)

                # recurse (n=3)
                P_nm1_1 = P_n_1

                # next polynomial of order 2 (n=3)
                P_nm2_2 = 0.0
                P_nm1_2 = 3 * sθ * sθ

                # n>2
                ρnm1 = ρ * ρ
                for n in 3:n_max
                    ITns[1,n] += ρnm1 * im * 0.5 * (-(n+1) * abs(P_nm1_0) + 1/n * abs(P_nm1_2) * e2imϕ)
                    ITns[2,n] += ρnm1 * 0.5 * (-(n+1) * abs(P_nm1_0) - 1/n * abs(P_nm1_2) * e2imϕ)

                    #--- recurse ---#

                    # ρ^(n-1)
                    ρnm1 *= ρ

                    # order 0
                    P_n_0 = ( (2*n-1) * cθ * P_nm1_0 - (n-1) * P_nm2_0 ) / n
                    P_nm2_0 = P_nm1_0
                    P_nm1_0 = P_n_0

                    # order 2
                    P_n_2 = ( (2*n-1) * cθ * P_nm1_2 - (n+1) * P_nm2_2 ) / (n-2)
                    P_nm2_2 = P_nm1_2
                    P_nm1_2 = P_n_2
                end
            end
        end
    end
    dV = dx * dy * dz
    V = dV * nx * nx * nx
    dV_V = dV / V
    ITns .*= dV_V

    return ITns
end

function preintegrate_multipole_χ(n_max, nθ, nϕ; s=1.0, θ_max=π, ϕ_max=2*π)
    ITns = zeros(Float64,2,2,n_max,nθ,nϕ)
    container = zeros(Complex{Float64},2,n_max)
    for (iϕ,ϕ) in enumerate(range(0.0, ϕ_max, nϕ))
        @show ϕ
        for (iθ,θ) in enumerate(range(0.0, θ_max, nθ))
            @show θ
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            r̂ = SVector{3}(sθ * cϕ, sθ * sϕ, cθ)
            container .= preintegrate_multipole_χ(r̂, n_max; s)
            for n in 1:n_max
                for i in 1:2 # z component is always zero
                    ITns[1,i,n,iθ,iϕ] = real(container[i,n])
                    ITns[2,i,n,iθ,iϕ] = imag(container[i,n])
                end
            end
        end
    end
    return ITns
end

const ε_Nθ_χ = 80
const ε_Δθ_χ = π / ε_Nθ_χ
const ε_Nϕ_χ = 20
const ε_Δϕ_χ = 2π / ε_Nϕ_χ

function get_iθ_χ(θr)
    θ = 0.0
    for i in 1:ε_Nθ_χ
        if θ >= θr
            θ-θr > θr-θ+ε_Δθ_χ && (return i-1)
            return i
        end
        θ += ε_Δθ_χ
    end
    π-θr > θr-π+ε_Δθ_χ && (return ε_Nθ_χ-1)
    return ε_Nθ_χ
end

function get_iϕ_χ(ϕr)
    if ϕr < 0
        ϕr += FastMultipole.π2
    end
    ϕ = 0.0
    for i in 1:ε_Nϕ_χ
        if ϕ >= ϕr
            ϕ-ϕr > ϕr-ϕ+ε_Δϕ_χ && (return i-1)
            return i
        end
        ϕ += ε_Δϕ_χ
    end
    FastMultipole.π2-ϕr > ϕr-FastMultipole.π2+ε_Δϕ_χ && (return ε_Nϕ_χ-1)
    return ε_Nϕ_χ
end

# IΧ2 = preintegrate_multipole_χ(30, ε_Nθ_χ, ε_Nϕ_χ)
# FastMultipole.write_file("src/multipole_integrals_chi.csv", IΧ2)
const IΧ = FastMultipole.read_file("src/multipole_integrals_chi.csv")

function multipole_check(x_target, n_sources; seed=123, expansion_order=5, method="final χ", χ_method = "loop", v_method="loop")
    # source system
    # source_system = generate_vortex(seed, n_sources)

    # Random.seed!(123)
    # position = zeros(3,n_sources)
    # position .= 0.7,0.4,0.5
    # position .= 0.5,0.5,0.5
    position = rand(3,n_sources) .* 2
    strength = rand(3,n_sources)
    # strength[1,1] = 1.0 # a single vortex of unit strength aligned with the x axis
                        # means the velocity should be in the -y direction
    strength = rand(3,n_sources)
    source_system = VortexParticles(position, strength)

    # source branch
    bodies_index, n_branches, branch_index, i_parent, i_leaf_index = 1:n_sources, 0, 1:0, 0, 1
    source_center = target_center = SVector{3}(0.5,0.5,0.5) * 2
    source_radius = target_radius = 0.5 * sqrt(3) * 2
    source_box = target_box = SVector{3}(0.5,0.5,0.5) * 2
    # source_center = target_center = SVector{3}(0.5,0.5,0.5)
    # source_radius = target_radius = 0.5 * sqrt(3)
    # source_box = target_box = SVector{3}(0.5,0.5,0.5)

    branch = Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, expansion_order+20)

    # get multipole coefficients
    FastMultipole.body_to_multipole!(branch, source_system, branch.harmonics, Val(expansion_order+20))
    # branch.multipole_expansion[:,1,:] .= 0.0 # we are testing just χ
    # branch.multipole_expansion[:,2,:] .= 0.0 # we are testing just ϕ

    # evaluate expansion
    ϕ, v, _ = evaluate_multipole(x_target, source_center, branch.multipole_expansion, DerivativesSwitch(true,true,false), Val(expansion_order), Val(true))
    ϕtrue, vtrue, _ = evaluate_multipole(x_target, source_center, branch.multipole_expansion, DerivativesSwitch(true,true,false), Val(expansion_order+20), Val(true))

    # check velocity
    v_analytic = SVector{3}(0.0,0,0)
    for i_body in branch.bodies_index
        dx = x_target-source_system[i_body,FastMultipole.Position()]
        q = source_system[i_body,FastMultipole.Strength()]
        r = norm(dx)
        gamma_over_R = q / r / (4*pi)
        jacobian = zeros(3,3)
        gamma_over_R /= r^2
        for j_potential in 1:3
            for i_r in 1:3
                jacobian[i_r,j_potential] -= gamma_over_R[j_potential] * dx[i_r]
            end
        end
        v_analytic += flatten_derivatives!(jacobian)
    end

    # rotate frame
    r⃗ = x_target - source_center
    r = norm(r⃗)
    r̂ = r⃗ / r
    x̂, ŷ, ẑ, R = rotate(r̂)
    @show R

    @assert isapprox(dot(x̂,ŷ),0.0; atol=1e-12)
    @assert isapprox(dot(x̂,ẑ),0.0; atol=1e-12)
    @assert isapprox(dot(ŷ,ẑ),0.0; atol=1e-12)

    ωx, ωy, ωz = FastMultipole.vortex_from_multipole(branch.multipole_expansion)
    ω⃗ = SVector{3}(ωx, ωy, ωz)
    ωxr = dot(ω⃗, x̂)
    ωyr = dot(ω⃗, ŷ)
    ωzr = dot(ω⃗, ẑ)
    ω⃗r = SVector{3}(ωxr, ωyr, ωzr)

    v̂x = dot(vtrue, x̂) / norm(vtrue)
    v̂y = dot(vtrue, ŷ) / norm(vtrue)

    # rotated system
    rotated_position = similar(position)
    rotated_strength = similar(strength)
    for i in 1:size(rotated_position,2)
        rotated_position[:,i] = source_center + R * (position[:,i] - source_center)
        rotated_strength[:,i] = R * strength[:,i]
    end
    rotated_system = VortexParticles(rotated_position, rotated_strength)

    rotated_branch = Branch(bodies_index, n_branches, branch_index, i_parent, i_leaf_index, source_center, target_center, source_radius, target_radius, source_box, target_box, expansion_order+20)

    # get multipole coefficients
    FastMultipole.body_to_multipole!(rotated_branch, rotated_system, branch.harmonics, Val(expansion_order+20))

    # evaluate expansion
    x_target_rotated = source_center + R * (x_target - source_center)
    ϕtruer, vtruer, _ = evaluate_multipole(x_target_rotated, source_center, rotated_branch.multipole_expansion, DerivativesSwitch(true,true,false), Val(expansion_order+20), Val(true))

    println("\nVerifying rotated χ")
    @show R * vtrue - vtruer

    # calculate v using single sum over n
    if method == "simplest ϕ"
        v̂ = v / norm(v)
        vmag = 0.0
        inm = 3
        r = norm(x_target - source_center)

        for n in 1:expansion_order
            ϕn1 = branch.multipole_expansion[1,1,inm] + branch.multipole_expansion[2,1,inm]*im
            vmag += ϕn1 * (n+1) / r^(n+2)
            inm += n+1
        end
        vmag *= (-v̂[2] + v̂[1]*im) / 4 / pi

    elseif method == "2 terms ϕ"

        v̂ = v / norm(v)
        vmag = 0.0
        inm = 3
        r = norm(x_target - source_center)

        for n in 1:expansion_order
            ϕn0 = branch.multipole_expansion[1,1,inm-1] + branch.multipole_expansion[2,1,inm-1]*im
            ϕn1 = branch.multipole_expansion[1,1,inm] + branch.multipole_expansion[2,1,inm]*im
            vmag += ϕn0 * (n+1) / r^(n+2) * -v̂[3]
            vmag += ϕn1 * (n+1) / r^(n+2) * (-v̂[2] + v̂[1]*im)
            inm += n+1
        end
        vmag /= 4 * pi

    elseif method == "3 terms ϕ"

        v̂ = v / norm(v)
        vmag = 0.0
        inm = 3
        r = norm(x_target - source_center)

        for n in 1:expansion_order
            ϕn0 = branch.multipole_expansion[1,1,inm-1] + branch.multipole_expansion[2,1,inm-1]*im
            ϕn1 = branch.multipole_expansion[1,1,inm] + branch.multipole_expansion[2,1,inm]*im
            ϕn_1 = -conj(ϕn1)
            vmag += ϕn0 * (n+1) / r^(n+2) * -v̂[3]
            vmag += ϕn1 * (n+1) / r^(n+2) * 0.5 * (-v̂[2] + v̂[1]*im)
            vmag += ϕn_1 * (n+1) / r^(n+2) * 0.5 * (v̂[2] + v̂[1]*im)
            inm += n+1
        end
        vmag /= 4 * pi

    elseif method == "final ϕ"

        vmag = 0.0
        v̂ = v / norm(v)
        r⃗ = x_target - source_center
        r, θ, ϕ = FastMultipole.cartesian_to_spherical(r⃗)
        inm = 1
        in0 = 1
        for n in 0:expansion_order
            # m=0
            ϕn0 = branch.multipole_expansion[1,1,in0] + branch.multipole_expansion[2,1,in0]*im
            Sn1m = factorial(n+1) / r^(n+2)
            vmag += ϕn0 * -v̂[3]*Sn1m

            # m=1
            ϕn1 = (branch.multipole_expansion[1,1,in0+1] + branch.multipole_expansion[2,1,in0+1]*im) * (n>0)
            vmag += ϕn1 * -((v̂[2] - im*v̂[1]) * factorial(n+1) / r^(n+2))

            # recurse indices
            in0 += n+1

        end

        vmag /= 4*pi

    elseif method == "final χ"

        vmag = 0.0
        v̂ = vtruer / norm(vtruer) # all quantities in rotated frame
        this_r⃗ = R * (x_target - source_center)
        r, θ, ϕ = FastMultipole.cartesian_to_spherical(this_r⃗)
        @show θ, ϕ
        in1 = 3
        for n in 1:expansion_order+20
            # m=1
            χn1 = (rotated_branch.multipole_expansion[1,2,in1] + rotated_branch.multipole_expansion[2,2,in1]*im)
            Sn1 = Float64(factorial(big(n))) / r^(n+1)
            Sn1_check = Snm(r,θ,ϕ,n,0)
            @assert isapprox(Sn1,Sn1_check;rtol=1e-10) "Sn1 = $(Sn1), Sn1_check=$Sn1_check"
            vmag += -χn1 * (v̂[1] + im*v̂[2]) * (n+1) * Float64(factorial(big(n))) / r^(n+1)

            # recurse indices
            in1 += n+1
        end

        vmag /= 4*pi

    end

    @show vmag norm(v) norm(vtrue) norm(vtruer)


    # calculate error
    εϕ, εχ = 0.0, 0.0

    # check error using known χ
    in1 = (((expansion_order+1) * (expansion_order + 2)) >> 1) + 2
    for n in expansion_order+1:expansion_order+2
        εχ += Complex{Float64}((branch.multipole_expansion[1,2,in1] + branch.multipole_expansion[2,2,in1]*im) * factorial(big(n+1)) / r^(n+1))
        @assert in1 == ((n*(n+1))>>1) + 2
        in1 += n+1
    end
    εχ *= -(v̂x + im*v̂y) / 4 / pi

    println("\nChecking error estimate with known χ")
    @show εχ
    @show norm(εχ) / norm(v-vtrue)

    ω⃗_check = SVector{3}(0.0,0,0)
    for i_source in branch.bodies_index
        ω⃗_check += source_system[i_source, FastMultipole.Strength()]
    end
    println("\nEstimating ω")
    @show sum(abs.(ω⃗ - ω⃗_check))
    @show sum(abs.(ω⃗r - R * ω⃗_check))
    @show ω⃗r

    εχ2 = 0.0
    εχ_loop = 0.0

    println("\nEstimating χ")
    # if χ_method == "preintegration"
    #     ITns = preintegrate_multipole_χ(r̂, n_max; s=1.0)
    # end

    for n in expansion_order+1:expansion_order+2 # +2
        in1 = (((n + 1) * n) >> 1) + 2

        if χ_method == "integration"
            nx = 100
            dx = source_box[1] * 2 / nx
            dy = source_box[2] * 2 / nx
            dz = source_box[3] * 2 / nx
            x0 = -source_box[1] + dx * 0.5
            y0 = -source_box[2] + dy * 0.5
            z0 = -source_box[3] + dz * 0.5
            ITnx = ITny = ITnz = 0.0 + 0im
            for x in range(x0, step=dx, length=nx)
                for y in range(y0, step=dy, length=nx)
                    for z in range(z0, step=dz, length=nx)
                        ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(R * SVector{3}(x,y,z))
                        cθ = cos(θ)
                        # ITnx += ρ^(n-1) * im * 0.5 * (-(n+1) * Plm(cos(θ),n-1,0) + 1/n * Plm(cos(θ),n-1,2) * exp(-2*im*ϕ))
                        # ITny += ρ^(n-1) * 0.5 * (-(n+1) * Plm(cos(θ),n-1,0) - 1/n * Plm(cos(θ),n-1,2) * exp(-2*im*ϕ))
                        # ITnz += ρ^(n-1) * (1 + 1/n) * im * Plm(cos(θ),n-1,1) * exp(-im*ϕ)
                        ITnx += ρ^(n-1) * im * 0.5 * (-(n+1) * abs(Plm(cos(θ),n-1,0)) + 1/n * abs(Plm(cos(θ),n-1,2)) * exp(-2*im*ϕ))
                        ITny += ρ^(n-1) * 0.5 * (-(n+1) * abs(Plm(cos(θ),n-1,0)) - 1/n * abs(Plm(cos(θ),n-1,2)) * exp(-2*im*ϕ))
                        ITnz += ρ^(n-1) * (1 + 1/n) * im * abs(Plm(cos(θ),n-1,1)) * exp(-im*ϕ)
                    end
                end
            end
            dV = dx * dy * dz
            V = dV * nx * nx * nx
            ITnx *= dV / V # / Float64(factorial(big(n+1)))
            ITny *= dV / V # / Float64(factorial(big(n+1)))
            ITnz *= dV / V # / Float64(factorial(big(n+1)))

            @show n, ITnx, ITny, ITnz

            # χn1 = Complex{Float64}(ωxr * ITnx + ωyr * ITny + ωzr * ITnz)
            χn1 = Complex{Float64}(ωxr * ITnx + ωyr * ITny) # + ωzr * ITnz)

        elseif χ_method == "preintegration"

            # s = mean(source_box)*2
            # χn1 = Complex{Float64}(ωxr * ITns[1,n] + ωyr * ITns[2,n] + ωzr * ITns[3,n])

            # get iθ, iϕ
            _, θr, ϕr = FastMultipole.cartesian_to_spherical(r⃗)
            iθ = get_iθ_χ(θr)
            iϕ = get_iϕ_χ(ϕr)
            println("\nVerifying lookup table indices:")
            @show (iθ-1) * ε_Δθ_χ, θr, iθ
            @show (iϕ-1) * ε_Δϕ_χ, ϕr, iϕ

            # estimate coefficient
            χn1_real = (ωxr * IΧ[1,1,n,iθ,iϕ] + ωyr * IΧ[1,2,n,iθ,iϕ]) # + ωzr * IΧ[1,3,n,iθ,iϕ])
            χn1_imag = (ωxr * IΧ[2,1,n,iθ,iϕ] + ωyr * IΧ[2,2,n,iθ,iϕ]) # + ωzr * IΧ[2,3,n,iθ,iϕ])
            χn1 = χn1_real + im * χn1_imag

            @show n IΧ[1,1,n,iθ,iϕ] + im*IΧ[2,1,n,iθ,iϕ] IΧ[1,2,n,iθ,iϕ] + im*IΧ[2,2,n,iθ,iϕ]# , IΧ[1,3,n,iθ,iϕ] + im*IΧ[2,3,n,iθ,iϕ]

            println()
        end
        # elseif χ_method == "loop"

            χn1_loop = 0.0
            for i_body in branch.bodies_index
                ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(R * (source_system[i_body,Position()] - source_center))
                this_ω⃗ = R * source_system[i_body, Strength()]
                χn1_loop += this_ω⃗[1] * ρ^(n-1) * im * 0.5 * (-(n+1) * Plm(cos(θ),n-1,0) + 1/n * Plm(cos(θ),n-1,2) * exp(-2*im*ϕ))
                χn1_loop += this_ω⃗[2] * ρ^(n-1) * 0.5 * (-(n+1) * Plm(cos(θ),n-1,0) - 1/n * Plm(cos(θ),n-1,2) * exp(-2*im*ϕ))
                χn1_loop += this_ω⃗[3] * ρ^(n-1) * (1 + 1/n) * im * Plm(cos(θ),n-1,1) * exp(-im*ϕ)
            end
            # χn1 /= factorial(big(n+1))
            # χn1 = Complex{Float64}(χn1)

        # end

        # check value
        χn1_check = rotated_branch.multipole_expansion[1,2,in1] + im * rotated_branch.multipole_expansion[2,2,in1]

        # @show χn1 / Float64(factorial(big(n+1))) χn1_check χn1 / Float64(factorial(big(n+1))) / χn1_check
        s = χ_method == "preintegration" ? mean(source_box) * 2 : 1.0
        @show χn1 / Float64(factorial(big(n+1))) / χn1_check * s^(n-1)
        @show χn1_loop / Float64(factorial(big(n+1))) / χn1_check

        @show s^(n-1) / r^(n+1) χn1
        εχ2 += χn1 / r^(n+1) * s^(n-1)
        εχ_loop += χn1_loop / r^(n+1)
        # εχ2 += χn1 * Float64(factorial(big(n+1))) / r^(n+1)
        in1 += n+1
    end
    εχ2 *= -(v̂x + im*v̂y) / 4 / pi
    εχ_loop *= -(v̂x + im*v̂y) / 4 / pi

    println("\nChecking error estimate with estimated χ")
    @show norm(εχ2) / norm(v-vtrue)
    @show norm(εχ_loop) / norm(v-vtrue)

    # test FastMultipole
    εχ_fmm = FastMultipole.multipole_error(r⃗, branch, expansion_order, FastMultipole.LambHelmholtzΧVelocity())

    @show εχ_fmm / norm(v-vtrue)

    return ϕ, v, v_analytic, εϕ, εχ, branch, source_system
end

x_target = SVector{3}(0.0,1.0,1.0) + SVector{3}(0.5,0.5,0.5)
x_target = -10*SVector{3}(1.0,0.8,1.2) + SVector{3}(0.5,0.5,0.5)
x_target *= 2
expansion_order = 5
n_sources = 30
seed = 123

ϕ, v, v_analytic, εϕ, εχ, branch, source_system = multipole_check(x_target, n_sources; seed, expansion_order, χ_method="integration")
println()

# errs = Float64[]
# for θ in range(0, stop=pi/2, length=10)
#     for ϕ in range(0, stop=pi, length=10)
#         x_target = SVector{3}(sin(θ)*cos(ϕ),sin(θ)*sin(ϕ),cos(θ))
#         x_target = x_target / norm(x_target) * 4.0
#         for seed in 1:1000
#             ϕ, v, v_analytic, εϕ, εχ = multipole_check(x_target, n_sources; seed, expansion_order)
#             push!(errs,norm(v-v_analytic)/(εϕ+εχ))
#         end
#     end
# end
# @show err10
