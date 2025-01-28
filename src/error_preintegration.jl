using Statistics
using LegendrePolynomials
#------- write/read files -------#

function write_file(fname, arr::AbstractArray)
	arr_size = size(arr)
	open(fname, "w") do io
        for s in arr_size
            print(io,s,",")
        end
        print(io,'\n')
		for i in eachindex(arr)
			println(io,arr[i])
		end
	end
end

function read_file(fname)
	arr = open(fname, "r") do io
        arr_size = Tuple(parse.(Int64,split(readline(io),','; keepempty=false)))
		arr = zeros(Float64, arr_size)
		i = 1
		for line in eachline(io)
            arr[i] = parse(Float64,strip(line))
			i += 1
		end
		return arr
	end
	return arr
end

function read_write_multipole(fpath)
    if isfile(fpath)
        multipole_integrals = read_file(fpath)
        if size(multipole_integrals) != (ε_MAX_N,ε_Nθ,ε_Nϕ)
            println("----- writing multipole error integrals -----")
            multipole_integrals = integrate_multipole(ε_MAX_N, 1.0, ε_NX, ε_Nθ, ε_Nϕ)
            write_file(fpath,multipole_integrals)
        end
    else
        println("----- writing multipole error integrals -----")
        multipole_integrals = integrate_multipole(ε_MAX_N, 1.0, ε_NX, ε_Nθ, ε_Nϕ)
        write_file(fpath,multipole_integrals)
    end
    return multipole_integrals
end

function read_write_local(fpath)
    if isfile(fpath)
        local_integrals = read_file(fpath)
        if size(local_integrals) != (ε_MAX_N,ε_Nω,ε_Nγ)
            println("----- writing local error integrals -----")
            local_integrals = integrate_local(ε_MAX_N, ε_Nω, ε_Nγ)
            write_file(fpath,local_integrals)
        end
    else
        println("----- writing local error integrals -----")
        local_integrals = integrate_local(ε_MAX_N, ε_Nω, ε_Nγ)
        write_file(fpath,local_integrals)
    end
    return local_integrals
end

function read_write_multipole_χ(fpath)
    if isfile(fpath)
        multipole_integrals_χ = read_file(fpath)
        if size(multipole_integrals_χ) != (3,ε_MAX_N,ε_Nθ_χ,ε_Nϕ_χ)
            println("----- writing multipole Lamb-Helmholtz error integrals -----")
            multipole_integrals_χ = preintegrate_multipole_χ(ε_MAX_N, ε_Nθ_χ, ε_Nϕ_χ; s=1.0, θ_max=π, ϕ_max=2*π)
            write_file(fpath, multipole_integrals_χ)
        end
    else
        println("----- writing multipole Lamb-Helmholtz error integrals -----")
        multipole_integrals_χ = preintegrate_multipole_χ(ε_MAX_N, ε_Nθ_χ, ε_Nϕ_χ; s=1.0, θ_max=π, ϕ_max=2*π)
        write_file(fpath, multipole_integrals_χ)
    end
    return multipole_integrals_χ
end

function read_write_local_ϕχ(fpath)
    if isfile(fpath)
        local_integrals = read_file(fpath)
        if size(local_integrals) != (5, ε_MAX_N, ε_Nθ_ϕχ, ε_Nϕ_ϕχ, ε_Nω_ϕχ)
            println("----- writing local Lamb-Helmholtz error integrals -----")
            local_integrals = integrate_local_ϕχ(ε_MAX_N, ε_Nθ_ϕχ, ε_Nϕ_ϕχ, ε_Nω_ϕχ)
            write_file(fpath,local_integrals)
        end
    else
        println("----- writing local Lamb-Helmholtz error integrals -----")
        local_integrals = integrate_local_ϕχ(ε_MAX_N, ε_Nθ_ϕχ, ε_Nϕ_ϕχ, ε_Nω_ϕχ)
        write_file(fpath,local_integrals)
    end
    return local_integrals
end

#------- multipole error -------#

function integrate_multipole(rvec, n_max, s, nx)
    dx = s/nx
    x = y = z = -0.5*s + dx*0.5
    R2 = s*s*0.25 # squared radius of largest enclosed sphere
    r = norm(rvec)
    r̂x, r̂y, r̂z = rvec / r # unit vector pointing to target
    rinv = 1/r

    res = zeros(n_max) # preallocate result
    for ix in 1:nx
        y = -0.5*s + dx*0.5
        for iy in 1:nx
            z = -0.5*s + dx*0.5
            for iz in 1:nx
                ρ2 = x*x + y*y + z*z
                if ρ2 > R2 # only integrate outside the sphere
                    ρ = sqrt(ρ2)
                    cθ = (r̂x * x + r̂y * y + r̂z * z) / ρ # cosine theta
                    ρn = ρ # ρ^n
                    Pnm1 = 1.0 # Legendre polynomial of degree n-1
                    Pn = cθ # Legendre polynomial of degree n
                    ρr = ρ * rinv
                    ρr_n = ρr

                    # update integral for each n
                    for n in 1:n_max
                        res[n] += ρr_n * abs(Pn) * rinv # actual integral for ρ/r
                        # res[n] += ρn * Pn

                        # next Legendre polynomial
                        Pnp1 = ((2n+1) * cθ * Pn - n * Pnm1) / (n+1)

                        # recurse
                        Pnm1 = Pn
                        Pn = Pnp1
                        ρn *= ρ
                        ρr_n *= ρr
                    end
                end
                z += dx # increment z
            end
            y += dx # increment y
        end
        x += dx # increment x
    end

    # finish integration
    res .*= dx * dx * dx / (s * s * s) # note that s=1 makes the division unnecessary

    return res
end

function integrate_multipole(n_max, s, nx, nθ, nϕ)
    res = zeros(n_max, nθ, nϕ)
    for (iθ,θ) in enumerate(range(0, stop=pi, length=nθ))
    # for (icθ,cθ) in enumerate(range(1.0, stop=-1.0, length=nθ))
        for (iϕ,ϕ) in enumerate(range(0, stop=pi, length=nϕ))
            # sθ = sqrt(1-cθ*cθ)
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            r = SVector{3}(sθ*cϕ, sθ*sϕ, cθ)
            # res[:,icθ,iϕ] .= multipole_preintegration(r, n_max, s, nx)
            res[:,iθ,iϕ] .= integrate_multipole(r, n_max, s, nx)
        end
    end
    return res
end

function get_iθ(θr)
    θ = 0.0
    if θr > π_over_2
        θr = π - θr
    end
    for i in 1:ε_Nθ
        if θ >= θr
            θ-θr > θr-θ+ε_Δθ && (return i-1)
            return i
        end
        θ += ε_Δθ
    end
    π_over_2-θr > θr-π_over_2+ε_Δθ && (return ε_Nθ-1)
    return ε_Nθ
end

function get_iϕ(ϕr)
    # ϕr < 0.0 && (ϕr = -ϕr)
    if ϕr > π
        ϕr = π2 - ϕr
    end
    ϕr = abs(ϕr)
    ϕ = 0.0
    for i in 1:ε_Nϕ
        if ϕ >= ϕr
            ϕ-ϕr > ϕr-ϕ+ε_Δϕ && (return i-1)
            return i
        end
        ϕ += ε_Δϕ
    end
    π-ϕr > ϕr-π+ε_Δϕ && (return ε_Nϕ-1)
    return ε_Nϕ
end

#------- local error -------#

# function integrate_local(dtheta::Float64, n::Int; R=1.0)
"""
Calculates the integral of Pn(γ)/ρ^(n+1) over a sphere at a distance of 1 and radius forming a cone of angle ω
"""
function integrate_local(ω, γ, n_max; R=1.0, nx=100)
    ω < 1e-6 && (ω += 1e-6)
    s_over_2 = R*sin(ω)
    dx = s_over_2*2/nx
    x0 = -s_over_2 + dx*0.5
    y0 = -s_over_2 + dx*0.5
    z0 = R - s_over_2 + dx*0.5
    cx, cy, cz = 0.0, 0.0, R
    sγ, cγ = sincos(γ)
    r_vec = SVector{3}(sγ,0.0,cγ) # assume π >= γ >= 0
    val = zeros(n_max)
    for x in x0:dx:x0+2*s_over_2
        for y in y0:dx:y0+2*s_over_2
            for z in z0:dx:z0+2*s_over_2
                if (x^2+y^2+(z-cz)^2) <= s_over_2 * s_over_2 # inside the sphere
                    # radial part
                    ρ_inv = 1/sqrt(x*x+y*y+z*z)
                    ρ_nm1 = ρ_inv*ρ_inv

                    # Legendre polynomial
                    cθ = (x*sγ + z*cγ) * ρ_inv
                    Pnm1 = 1.0 # Legendre polynomial of degree n-1
                    Pn = cθ # Legendre polynomial of degree n

                    for n in 1:n_max
                        # val[n] += rho^(-n-1)
                        val[n] += ρ_nm1 * abs(Pn)

                        # next Legendre polynomial
                        Pnp1 = ((2n+1) * cθ * Pn - n * Pnm1) / (n+1)

                        # recurse
                        Pnm1 = Pn
                        Pn = Pnp1

                        ρ_nm1 *= ρ_inv
                    end
                end
            end
        end
    end
    val .*= dx*dx*dx / (8*s_over_2*s_over_2*s_over_2)
    return val
end

function integrate_local(n_max::Int, n_ω::Int, n_γ)
    res = zeros(n_max, n_ω, n_γ)
    for (iγ,γ) in enumerate(range(0, stop=π, length=n_γ))
        for (iω,ω) in enumerate(range(0, stop=π_over_2, length=n_ω))
            res[:,iω,iγ] .= integrate_local(ω,γ,n_max)
        end
    end
    return res
end

function get_iω(ωr)
    ω = 0.0
    for i in 1:ε_Nω
        if ω >= ωr
            ω-ωr > ωr-ω+ε_Δω && (return i-1)
            return i
        end
        ω += ε_Δω
    end
    π_over_2-ωr > ωr-π_over_2+ε_Δω && (return ε_Nω-1)
    return ε_Nω
end

function get_iγ(γr)
    γ = 0.0
    for i in 1:ε_Nγ
        if γ >= γr
            γ-γr > γr-γ+ε_Δγ && (return i-1)
            return i
        end
        γ += ε_Δγ
    end
    π-γr > γr-π+ε_Δγ && (return ε_Nγ-1)
    return ε_Nγ
end

#------- Lamb-Helmholtz multipole error -------#

function rotate(r̂)
    ẑ = r̂
    x̂ = SVector{3,eltype(r̂)}(1.0,0,0) - dot(SVector{3,eltype(r̂)}(1.0,0,0), r̂) * r̂
    xnorm = zero(eltype(r̂))
    xnorm += abs(x̂[1])
    xnorm += abs(x̂[2])
    xnorm += abs(x̂[3])
    if xnorm < 1e-1
        x̂ = SVector{3,eltype(r̂)}(0.0,1.0,0.0) - dot(SVector{3,eltype(r̂)}(0.0,1,0), r̂) * r̂
    end
    x̂ /= norm(x̂)
    ŷ = cross(ẑ, x̂)
    # R = vcat(x̂', ŷ', ẑ')
    R = SMatrix{3,3,eltype(r̂)}(x̂[1], ŷ[1], ẑ[1], x̂[2], ŷ[2], ẑ[2], x̂[3], ŷ[3], ẑ[3])
    return R
end

@inline function cartesian_to_spherical(x; EPSILON=1e-10)
    return cartesian_to_spherical(x[1], x[2], x[3]; EPSILON)
end

@inline function cartesian_to_spherical(x, y, z; EPSILON=1e-10)
    x2y2 = x*x + y*y
    r2 = x2y2 + z*z
    r = iszero(r2) ? r2 : sqrt(r2)
    z_r = z/r
    if r > 0
        theta = x2y2 > 0 ? acos(z_r) : π * (z < 0)
    else
        theta = zero(r)
    end
    phi = iszero(x2y2) ? zero(x2y2) : atan(y, x)
    return r, theta, phi
end

function preintegrate_multipole_χ(r̂, n_max; s=1.0, abs_val=true)
    # rotate frame
    R = rotate(r̂)

    # preallocate results
    ITns = zeros(Complex{Float64},3,n_max)
    # ITns = zeros(Complex{Float64},7,n_max)

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
                eimϕ = exp(-im*ϕ)
                my_eimϕ = abs(real(eimϕ)) + im * abs(imag(eimϕ))
                e2imϕ = eimϕ * eimϕ
                my_e2imϕ = abs(real(e2imϕ)) + im * abs(imag(e2imϕ))

                # Legendre polynomial recursions
                P_nm1_0 = one(ρ)
                P_n_0 = cθ
                P_n_1 = -sθ

                # n=1
                if abs_val
                    n = 1
                    ITns[1,1] += im * abs(P_nm1_0)
                    ITns[2,1] += abs(P_nm1_0)
                    # ITns[4,1] += ρ * abs(P_n_1) * my_eimϕ * 0.5 # Φₙ0
                    # tmp = n * 0.5 * ρ * abs(P_n_0) * 0.5
                    # ITns[5,1] += tmp # Φₙ1x
                    # ITns[6,1] += im * tmp # Φₙ1y
                    # ITns[7,1] += ρ / (n+1) * P_n_1 * my_eimϕ # Φₙ1z
                else
                    ITns[1,1] += -im * P_nm1_0
                    ITns[2,1] += -P_nm1_0
                    # ITns[4,1] += ρ * P_n_1 * conj(eimϕ)
                    # ITns[5,1] += n * 0.5 * ρ * P_n_0
                    # ITns[6,1] += -im * n * 0.5 * ρ * P_n_0
                    # ITns[7,1] += ρ / (n+1) * P_n_1 * eimϕ
                end

                # recurse
                P_nm1_0 = cθ
                P_nm2_0 = one(cθ)
                P_n_0 = (3 * cθ * P_nm1_0 - P_nm2_0) * 0.5
                # @assert isapprox(P_n_0, Plm(cθ,2,0); atol=1e-12)

                P_np1_1 = 3 * cθ * P_n_1 # n=1 just for this one
                P_nm1_1 = P_n_1 # n=2 as usual
                P_n_1 = P_np1_1 # n=2 as usual
                # @assert isapprox(P_n_1, Plm(cθ,2,1); atol=1e-12)

                P_n_2 = 3 * sθ * sθ

                # n=2

                n=2
                if abs_val
                    ITns[1,2] += 1.5 * ρ * im * abs(P_nm1_0) * 0.5 # 0.5 factor because P_1^0 is odd
                    ITns[2,2] += 1.5 * ρ * abs(P_nm1_0) * 0.5 # 0.5 factor because P_1^0 is odd
                    tmp = 1.5 * ρ * abs(P_nm1_1) * my_eimϕ * 0.5 # 0.5 factor because eimϕ is odd
                    ITns[3,2] += imag(tmp) + im * real(tmp)
                    # ITns[4,1] += ρ * ρ * abs(P_n_1) * my_eimϕ * 0.5
                    # tmp = n * 0.5 * ρ * ρ * abs(P_n_0) + 0.5 / (n+1) * ρ * ρ * abs(P_n_2) * my_e2imϕ * 0.5
                    # ITns[5,1] += tmp
                    # ITns[6,1] += im * real(tmp) + imag(tmp)
                    # ITns[7,1] += ρ * ρ / (n+1) * abs(P_n_1) * my_eimϕ * 0.5
                else
                    ITns[1,2] += -1.5 * ρ * im * P_nm1_0
                    ITns[2,2] += -1.5 * ρ * P_nm1_0
                    ITns[3,2] += 1.5 * ρ * im * P_nm1_1 * eimϕ
                    ITns[4,1] += ρ * ρ * P_n_1 * conj(eimϕ)
                    ITns[5,1] += n * 0.5 * ρ * ρ * P_n_0 + 0.5 / (n+1) * ρ * ρ * P_n_2 * e2imϕ
                    ITns[6,1] += -im * n * 0.5 * ρ * ρ * P_n_0 + im * 0.5 / (n+1) * ρ * ρ * P_n_2 * e2imϕ
                    ITns[7,1] += ρ * ρ / (n+1) * P_n_1 * eimϕ
                end

                # next polynomial of order 0
                n = 2
                P_n_0 = ( (2*n-1) * cθ * P_nm1_0 - (n-1) * P_nm2_0 ) / n

                # recurse (n=3)
                P_nm2_0 = P_nm1_0
                P_nm1_0 = P_n_0
                P_n_0 = (5*cθ*P_nm1_0 - 2 * P_nm2_0) / 3
                # @assert isapprox(P_n_0, Plm(cθ,3,0); atol=1e-12)

                # next polynomials of order 1 (n=2)
                P_np1_1 = 0.5 * (5*cθ*P_n_1 - 3 * P_nm1_1)

                # recurse (n=3)
                P_nm1_1 = P_n_1
                P_n_1 = P_np1_1

                # @assert isapprox(P_nm1_1, Plm(cθ,2,1); atol=1e-12)
                # @assert isapprox(P_n_1, Plm(cθ,3,1); atol=1e-12)

                # next polynomial of order 2 (n=3)
                P_nm2_2 = 0.0
                P_nm1_2 = 3 * sθ * sθ
                P_n_2 = 5*cθ*P_nm1_2 - 4*P_nm2_2

                # n>2
                ρnm1 = ρ * ρ
                for n in 3:n_max
                    if abs_val
                        fac = iseven(n) ? 0.5 : 1.0 # whenever the function is odd
                        fac2 = isodd(n) ? 0.5 : 1.0
                        tmp = ρnm1 * 0.5 * ((n+1) * abs(P_nm1_0) * fac + 1/n * abs(P_nm1_2) * my_e2imϕ * 0.5) # 0.5 because of e2imϕ
                        ITns[1,n] += imag(tmp) + im * real(tmp)
                        ITns[2,n] += tmp
                        tmp = ρnm1 * (1+1/n) * abs(P_nm1_1) * my_eimϕ * 0.5
                        ITns[3,n] += imag(tmp) + im * real(tmp)
                        # ITns[4,n] += ρnm1 * ρ * abs(P_n_1) * my_eimϕ * 0.5
                        # tmp = n * 0.5 * ρnm1 * ρ * abs(P_n_0) * fac2 + 0.5 / (n+1) * ρnm1 * ρ * abs(P_n_2) * my_e2imϕ * 0.5
                        # ITns[5,n] += tmp
                        # ITns[6,n] += im * real(tmp) + imag(tmp)
                        # ITns[7,n] += ρnm1 * ρ / (n+1) * abs(P_n_1) * my_eimϕ * 0.5
                    else
                        ITns[1,n] += ρnm1 * im * 0.5 * (-(n+1) * P_nm1_0 + 1/n * P_nm1_2 * e2imϕ)
                        ITns[2,n] += ρnm1 * 0.5 * (-(n+1) * P_nm1_0 - 1/n * P_nm1_2 * e2imϕ)
                        ITns[3,n] += ρnm1 * im * (1+1/n) * P_nm1_1 * eimϕ
                        # ITns[4,n] += ρnm1 * ρ * P_n_1 * conj(eimϕ)
                        # ITns[5,n] += n * 0.5 * ρnm1 * ρ * P_n_0 + 0.5 / (n+1) * ρnm1 * ρ * P_n_2 * e2imϕ
                        # ITns[6,n] += -im * n * 0.5 * ρnm1 * ρ * P_n_0 + im * 0.5 / (n+1) * ρnm1 * ρ * P_n_2 * e2imϕ
                        # ITns[7,n] += ρnm1 * ρ / (n+1) * P_n_1 * eimϕ
                    end

                    #--- recurse ---#

                    # ρ^(n-1)
                    ρnm1 *= ρ

                    # order 0
                    P_np1_0 = ((2*n+1) * cθ * P_n_0 - n * P_nm1_0) / (n+1)
                    P_nm2_0 = P_nm1_0
                    P_nm1_0 = P_n_0
                    P_n_0 = P_np1_0
                    # @assert isapprox(P_np1_0, Plm(cθ,n+1,0); atol=1e-12)

                    # order 1
                    P_np1_1 = ((2*n+1) * cθ * P_n_1 - (n+1) * P_nm1_1) / n
                    P_nm1_1 = P_n_1
                    P_n_1 = P_np1_1
                    # @assert isapprox(P_np1_1, Plm(cθ,n+1,1); atol=1e-12)

                    # order 2
                    P_np1_2 = ((2*n+1) * cθ * P_n_2 - (n+2) * P_nm1_2) / (n-1)
                    P_nm2_2 = P_nm1_2
                    P_nm1_2 = P_n_2
                    P_n_2 = P_np1_2
                    # @assert isapprox(P_np1_2, Plm(cθ,n+1,2); atol=1e-12)
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
    ITns = zeros(Float64,3,n_max,nθ,nϕ)
    container = zeros(Complex{Float64},3,n_max)
    for (iϕ,ϕ) in enumerate(range(0.0, ϕ_max, nϕ))
        println("iϕ = $iϕ / $nϕ")
        for (iθ,θ) in enumerate(range(0.0, θ_max, nθ))
            println("iθ = $iθ / $nθ")
            sθ, cθ = sincos(θ)
            sϕ, cϕ = sincos(ϕ)
            r̂ = SVector{3}(sθ * cϕ, sθ * sϕ, cθ)
            container .= preintegrate_multipole_χ(r̂, n_max; s)
            for n in 1:n_max
                for i in 1:3 # χ1x,χ1y,χ1z,ϕ0,ϕ0x,ϕ0y,ϕ0z
                    # ITns[1,i,n,iθ,iϕ] = real(container[i,n])
                    # ITns[2,i,n,iθ,iϕ] = imag(container[i,n])
                    ITns[i,n,iθ,iϕ] = abs(container[i,n])
                end
            end
        end
    end
    return ITns
end

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

#------- Lamb-Helmholtz Local Error -------#

function integrate_local_ϕχ(θ, ϕ, ω, n_max; R=1.0, nx=100)
    # preallocate integrals
    res = zeros(Complex{Float64},7,n_max)

    # perform integration
    integrate_local_ϕχ!(res, θ, ϕ, ω, n_max; R=1.0, nx=100)
    return res
end

"""
Calculates the integrals over a source sphere at a distance of R, polar and azimutha angles θ and ϕ, and radius forming a cone of angle ω. Assume the evaluation point is aligned with the z axis.
"""
function integrate_local_ϕχ!(res, θ, ϕ, ω, n_max; R=1.0, nx=100)
    ω < 1e-6 && (ω += 1e-6)

    # get target radius
    target_radius = R*sin(ω)
    dx = 2 * target_radius / nx

    # get unit vectors
    sθ, cθ = sincos(θ)
    sϕ, cϕ = sincos(ϕ)
    ΔC = R * SVector{3}(sθ * cϕ, sθ * sϕ, cθ)
    Ĉz = ΔC / norm(ΔC)
    Ĉx = SVector{3}(1.0,0,0) - Ĉz * Ĉz[1]
    if abs(Ĉx[1]) + abs(Ĉx[2]) + abs(Ĉx[3]) < 1e-1
        Ĉx = SVector{3}(0,1.0,0) - Ĉz * Ĉz[2]
    end
    Ĉx /= norm(Ĉx)
    Ĉy = cross(Ĉz, Ĉx)

    # integration bounds
    x⃗ = -target_radius * Ĉx + dx * 0.5 * Ĉx
    δx⃗ = Ĉx * dx
    y⃗0 = -target_radius * Ĉy + dx * 0.5 * Ĉy
    δy⃗ = Ĉy * dx
    z⃗0 = -target_radius * Ĉz + dx * 0.5 * Ĉz
    δz⃗ = Ĉz * dx

    # perform integration
    for ix in 1:nx
        y⃗ = y⃗0
        for iy in 1:nx
            z⃗ = z⃗0
            for iz in 1:nx
                ρ⃗ = ΔC + x⃗ + y⃗ + z⃗
                if norm(ρ⃗ - ΔC) <= target_radius
                    ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(ρ⃗)
                    ρnp1 = ρ * ρ
                    ρnp2 = ρ * ρnp1
                    eimϕ = exp(-im*ϕ)
                    my_eimϕ = abs(real(eimϕ)) + im * abs(imag(eimϕ))
                    e2imϕ = eimϕ * eimϕ
                    my_e2imϕ = abs(real(e2imϕ)) + im * abs(imag(e2imϕ))
                    sθ, cθ = sincos(θ)

                    # legendre polynomials, n=1
                    P_nm1_0 = 1.0
                    P_n_0 = cθ

                    P_nm1_1 = 0.0
                    P_n_1 = -sθ

                    P_nm1_2 = 0.0
                    P_n_2 = 0.0
                    P_np1_2 = 3 * sθ * sθ

                    # update integrals
                    for n in 1:n_max

                        # next legendre polynomials
                        P_np1_0 = ( (2*n+1) * cθ * P_n_0 - n * P_nm1_0 ) / (n+1)
                        P_np1_1 = ( (2*n+1) * cθ * P_n_1 - (n+1) * P_nm1_1 ) / n
                        P_np2_2 = ( (2*n+3) * cθ * P_np1_2 - (n+3) * P_n_2 ) / n

                        # update integrals
                        res[1,n] += im / ρnp1 * abs(P_n_1) * (real(my_eimϕ) - im*imag(my_eimϕ)) # Φₙ
                        res[2,n] += (n+1) * 0.5 / ρnp1 * abs(P_n_0) # Φₙ0
                        res[3,n] += 1.0 / (n*ρnp1) * abs(P_n_1) * my_eimϕ # Φₙ1
                        res[4,n] += 1.0 / (2*n*ρnp1) * abs(P_n_2) * my_e2imϕ # Φₙ2
                        res[5,n] += n * 0.5 / ρnp2 * abs(P_np1_0) # Χₙ0
                        res[6,n] += n / ((n+1)*ρnp2) * abs(P_np1_1) * my_eimϕ # Χₙ1
                        res[7,n] += 0.5 / ((n+1)*ρnp2) * abs(P_np1_2) * my_e2imϕ # Χₙ2

                        # recurse
                        ρnp1 *= ρ
                        ρnp2 *= ρ

                        # Pn0
                        P_nm1_0 = P_n_0
                        P_n_0 = P_np1_0

                        # Pn1
                        P_nm1_1 = P_n_1
                        P_n_1 = P_np1_1

                        # Pn2
                        P_nm1_2 = P_n_2
                        P_n_2 = P_np1_2
                        P_np1_2 = P_np2_2

                    end
                end
                z⃗ += δz⃗
            end
            y⃗ += δy⃗
        end
        x⃗ += δx⃗
    end

    # volume
    V = 4/3*pi*target_radius*target_radius*target_radius

    res .*= dx * dx * dx / V

    return res
end

function integrate_local_ϕχ(n_max::Int, n_θ::Int, n_ϕ::Int, n_ω::Int)
    res = zeros(5, n_max, n_θ, n_ϕ, n_ω)
    integrals = zeros(Complex{Float64},7,n_max)
    for (iω,ω) in enumerate(range(0, stop=π_over_2, length=n_ω))
        println("iω = $iω / $n_ω")
        for (iϕ,ϕ) in enumerate(range(0,stop=π2, length=n_ϕ))
            println("iϕ = $iϕ / $n_ϕ")
            for (iθ,θ) in enumerate(range(0,stop=π,length=n_θ))
                # integrate
                integrate_local_ϕχ!(integrals, θ,ϕ,ω,n_max)

                # save Φₙx, Φₙy, Φₙz
                for n in 1:n_max
                    Φₙ = integrals[1,n]
                    Φₙ0 = integrals[2,n]
                    Φₙ1 = integrals[3,n]
                    Φₙ2 = integrals[4,n]
                    Χₙ0 = integrals[5,n]
                    Χₙ1 = integrals[6,n]
                    Χₙ2 = integrals[7,n]
                    res[1,n,iθ,iϕ,iω] = real(Φₙ0) + imag(Φₙ0) + real(Φₙ2) + imag(Φₙ2) + real(Φₙ)
                    res[2,n,iθ,iϕ,iω] = real(Φₙ0) + imag(Φₙ0) + real(Φₙ2) + imag(Φₙ2) + imag(Φₙ)
                    res[3,n,iθ,iϕ,iω] = real(Φₙ1) + imag(Φₙ1)
                    res[4,n,iθ,iϕ,iω] = abs(Χₙ0+Χₙ2)
                    res[5,n,iθ,iϕ,iω] = abs(Χₙ1)
                end
            end
        end
    end
    return res
end

function get_iθ_ϕχ(θr)
    θ = 0.0
    for i in 1:ε_Nθ_ϕχ
        if θ >= θr
            θ-θr > θr-θ+ε_Δθ_ϕχ && (return i-1)
            return i
        end
        θ += ε_Δθ_ϕχ
    end
    π_over_2-θr > θr-π_over_2+ε_Δθ_ϕχ && (return ε_Nθ_ϕχ-1)
    return ε_Nθ_ϕχ
end

function get_iϕ_ϕχ(ϕr)
    ϕ = 0.0
    ϕr = abs(ϕr)
    for i in 1:ε_Nϕ_ϕχ
        if ϕ >= ϕr
            ϕ-ϕr > ϕr-ϕ+ε_Δϕ_ϕχ && (return i-1)
            return i
        end
        ϕ += ε_Δϕ_ϕχ
    end
    π_over_2-ϕr > ϕr-π_over_2+ε_Δϕ_ϕχ && (return ε_Nϕ_ϕχ-1)
    return ε_Nϕ_ϕχ
end

function get_iω_ϕχ(ωr)
    ω = 0.0
    for i in 1:ε_Nω_ϕχ
        if ω >= ωr
            ω-ωr > ωr-ω+ε_Δω_ϕχ && (return i-1)
            return i
        end
        ω += ε_Δω_ϕχ
    end
    π_over_2-ωr > ωr-π_over_2+ε_Δω_ϕχ && (return ε_Nω_ϕχ-1)
    return ε_Nω_ϕχ
end

