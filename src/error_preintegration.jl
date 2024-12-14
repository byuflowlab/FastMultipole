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
            multipole_integrals = integrate_multipole(ε_MAX_N, 1.0, ε_NX, ε_Nθ, ε_Nϕ)
            write_file(fpath,multipole_integrals)
        end
    else
        multipole_integrals = integrate_multipole(ε_MAX_N, 1.0, ε_NX, ε_Nθ, ε_Nϕ)
        write_file(fpath,multipole_integrals)
    end
    return multipole_integrals
end

function read_write_local(fpath)
    if isfile(fpath)
        local_integrals = read_file(fpath)
        if size(local_integrals) != (ε_MAX_N,ε_Nω,ε_Nγ)
            local_integrals = integrate_local(ε_MAX_N, ε_Nω, ε_Nγ)
            write_file(fpath,local_integrals)
        end
    else
        local_integrals = integrate_local(ε_MAX_N, ε_Nω, ε_Nγ)
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
    θr > π_over_2 && (θr = π - θr)
    for i in 1:ε_Nθ
        if θ >= θr
            θ-θr > θr-θ+ε_Δθ && (return i-1)
            return i
        end
        θ += ε_Δθ
    end
    π-θr > θr-π+ε_Δθ && (return ε_Nθ-1)
    return ε_Nθ
end

function get_iϕ(ϕr)
    # ϕr < 0.0 && (ϕr = -ϕr)
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
