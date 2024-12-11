#--- choose expansion order based on error tolerance ---#

function get_P(Δx, Δy, Δz, target_branch, source_branch, expansion_order::Int, error_method::ErrorMethod)
    return expansion_order
end

function get_P(Δx, Δy, Δz, target_branch, source_branch, ::Dynamic{PMAX,ε}, error_method::Union{EqualSpheres, UnequalSpheres, UnequalBoxes, UniformUnequalSpheres, UniformUnequalBoxes}) where {PMAX,ε}
    ΔC2 = Δx*Δx+Δy*Δy+Δz*Δz
    r_min, r_max, ρ_min, ρ_max = get_r_ρ(target_branch, source_branch, ΔC2)
    A = get_A(source_branch, SVector{3}(Δx,Δy,Δz))
    return get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, PMAX, A, ε, error_method)
end

function get_P(Δx, Δy, Δz, target_branch, source_branch, ::Dynamic{PMAX,ε}, error_method::UniformCubes) where {PMAX,ε}
return nothing
end

@inline function warn_Pmax()
    if WARNING_FLAG_PMAX[]
        @warn "error tolerance not met with dynamic expansion order; try increasing `Pmax`"
        WARNING_FLAG_PMAX[] = false
    end
end


"""
    get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, Pmax, ε_rel, error_method)

Returns the smallest expansion order not greater than `Pmax` and satisfying the specified relative error tolerance.

# Inputs

* `r_min::Flaot64`: distance from the multipole expansion to the closest target
* `r_max::Float64`: distance from the local expansion to the farthest target
* `ρ_min::Float64`: distance from the local expansion to the closest source
* `ρ_max::Float64`: distance from the multipole expansion to the farthest source
* `ΔC2::Float64`: distance squared between multipole and local centers
* `Pmax::Int64`: maximum allowable expansion order
* `ε_rel::Float64`: relative error tolerance
* `error_method::ErrorMethod`: type used to dispatch on the desired error method

# Ouputs

* `P::Int`: the smallest expansion order to satisfy the error tolerance

"""
function get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, Pmax, A, ε, error_method::Union{EqualSpheres, UnequalSpheres, UnequalBoxes})
    ρ_max_over_r_min = ρ_max / r_min
    r_max_over_ρ_min = r_max / ρ_min
    t1 = A / (r_min - ρ_max) * ρ_max_over_r_min # absolute error
    # t1 = ρ_max / (r_min - ρ_max) # relative error
    t2 = A / (ρ_min - r_max) * r_max_over_ρ_min # absolute error
    # t2 = ρ_min / (ρ_min - r_max) # relative error
    for P in 0:Pmax-1
        t1 + t2 < ε && (return P)
        t1 *= ρ_max_over_r_min
        t2 *= r_max_over_ρ_min
    end

    warn_Pmax()

    return Pmax
end

function get_P(r_min, r_max, ρ_min, ρ_max, ΔC2, Pmax, A, ε, ::Union{UniformUnequalSpheres, UniformUnequalBoxes})

    # multipole error
    ρ_max_over_r_min = ρ_max / r_min
    # ε_multipole = 1.5 * ρ_max / (r_min - ρ_max) # relative error
    ε_multipole = 1.5 * A / (r_min - ρ_max) * ρ_max_over_r_min # absolute error

    # local error
    ρ_max2 = ρ_max * ρ_max
    # Γ = 3 * r_max * r_min / (2 * ρ_max2 * ρ_max) # relative error
    Γ = 3 * A * r_max / (2 * ρ_max2 * ρ_max) # absolute error

    # distance from multipole center to point closest to the local expansion
    ΔC = sqrt(ΔC2)
    ρ_max_line = ΔC - ρ_min

    # distance from local center to farthest multipole point
    ΔC_plus = ρ_min + ρ_max_line + ρ_max_line

    # distance from local center to nearest multipole point
    ΔC_minus = ρ_min
    t_plus = r_max / ΔC_plus
    one_over_ΔC_minus = 1/ΔC_minus
    t_minus = r_max * one_over_ΔC_minus
    L = log((1-t_plus)/(1-t_minus))
    ΔC2_inv = 1/(2*ΔC)
    ρ2_ΔC2_2ΔC = (ρ_max2 - ΔC2) * ΔC2_inv
    r2_ΔC2 = r_max * r_max * ΔC2_inv

    # recursively tracked quantities
    Lζ_sum_pm1 = L
    Lζ_sum_pm2 = L
    Lζ_sum_pm3 = L
    t_plus_n = t_plus
    t_minus_n = t_minus
    η_p = log(ΔC_plus * one_over_ΔC_minus)
    r_inv = 1/r_max
    η_pm1 = η_p + (ΔC_plus - ΔC_minus) * r_inv
    η_pm2 = η_pm1 + (ΔC_plus*ΔC_plus - ΔC_minus*ΔC_minus) * r_inv * r_inv * 0.5

    #--- test p=0 ---#

    # test error
    ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))
    ε_multipole * 0.25 + ε_local < ε_rel && (return 0)

    # recurse multipole
    ε_multipole *= ρ_max_over_r_min

    # recurse local
    η_pm2 = η_pm1
    η_pm1 = η_p
    η_p = zero(r_min)

    #--- test p>0 ---#

    for P in 1:Pmax-1
        # get local error
        ε_local = Γ * (r_max * (η_pm1 + Lζ_sum_pm2) + ρ2_ΔC2_2ΔC * (η_p + Lζ_sum_pm1) - r2_ΔC2 * (η_pm2 + Lζ_sum_pm3))

        # check error
        ε_multipole / (P+4) + ε_local < ε_rel && (return P)

        # recurse multipole
        ε_multipole *= ρ_max_over_r_min

        # recurse local
        Lζ_sum_pm3 = Lζ_sum_pm2
        Lζ_sum_pm2 = Lζ_sum_pm1
        Lζ_sum_pm1 += (t_plus_n - t_minus_n) / P
        t_plus_n *= t_plus
        t_minus_n *= t_minus
        η_pm2 = η_pm1
        η_pm1 = η_p
        η_p = zero(r_min)
    end

    warn_Pmax()

    return Pmax
end

#--- dispatching dynamic expansion order ---#

@inline function get_Pmax(expansion_order::Int)
    return expansion_order
end

@inline function get_Pmax(expansion_order::Dynamic{PMAX,<:Any}) where PMAX
    return PMAX
end
