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
        @show fpath
        local_integrals = read_file(fpath)
        if size(local_integrals) != (ε_MAX_N,ε_NΔθ)
            local_integrals = integrate_local(ε_MAX_N, ε_NΔθ)
            write_file(fpath,local_integrals)
        end
    else
        local_integrals = integrate_local(ε_MAX_N, ε_NΔθ)
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
    dV = dx * dx * dx # differential volume
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
                        res[n] += ρr_n * Pn * rinv # actual integral for ρ/r
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
    res .*= dV

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

function integrate_local(dtheta::Float64, n::Int; R=1.0)
    dtheta < 1e-6 && (dtheta += 1e-6)
    drho = R*sin(dtheta)
    nx = 30
    dx = drho*2/nx
    x0 = -drho + dx*0.5
    y0 = -drho + dx*0.5
    z0 = R - drho + dx*0.5
    cx, cy, cz = 0.0, 0.0, R
    val = 0.0
    for x in x0:dx:x0+2*drho
        for y in y0:dx:y0+2*drho
            for z in z0:dx:z0+2*drho
                if (x^2+y^2+(z-cz)^2) <= drho*drho
                    rho = sqrt(x*x+y*y+z*z)
                    val += rho^(-n-1)
                end
            end
        end
    end
    val *= dx*dx*dx
    return val
end

function integrate_local(n_max::Int, nθ::Int)
    res = zeros(n_max, nθ)
    for (iθ,Δθ) in enumerate(range(0, stop=pi/2, length=nθ))
        for n in 1:n_max
            res[n,iθ] = integrate_local(Δθ,n)
        end
    end
    return res
end

function get_iΔθ(Δθr)
    Δθ = 0.0
    for i in 1:ε_NΔθ
        if Δθ >= Δθr
            Δθ-Δθr > Δθr-Δθ+ε_dΔθ && (return i-1)
            return i
        end
        Δθ += ε_dΔθ
    end
    pi05-Δθr > Δθr-pi05+ε_dΔθ && (return ε_NΔθ-1)
    return ε_NΔθ
end
