# Defines rrules for ChainRules and registers them with ReverseDiff.
tracked_type = Union{ReverseDiff.TrackedReal,ReverseDiff.TrackedArray} # tracked type used for registering ChainRulesCore rrules with ReverseDiff.

"""
function regular_harmonic!(harmonics, rho, theta, phi, P)
"""

function ChainRulesCore.rrule(::typeof(regular_harmonic!), _harmonics, rho, theta, phi, P)

    harmonics = copy(_harmonics)
    y,x = sincos(theta)
    fact = 1.0
    pl = 1.0
    p = 0.0
    rhom = 1.0 # rho^l / (l+m)! * (-1)^l
    ei = exp(im*phi)
    eim = 1.0

    pl = 1.0
    lpl = 1
    lpm = 1
    p1 = 1.0
    #_p = zeros(eltype(theta),P+1) # used in recurrance relation for derivatives
    #_pl = zeros(eltype(theta),P+1) # used in recurrance relation for derivatives
    for m=0:P # l=m up here
        p = pl
        lpl = m * m + 2 * m + 1
        lml = m * m + 1

        harmonics[lpl] = rhom*p*eim
        harmonics[lml] = conj(harmonics[lpl])

        p1 = p
        p = x * (2 * m + 1) * p1
        #_pl[m+1] = pl # m starts at zero so I have to adjust these indices
        #_p[m+1] = p
        rhom *= rho
        rhol = rhom
        for l=m+1:P # l>m in here
            lpm = l * l + l + m + 1
            lmm = l * l + l - m + 1
            rhol /= -(l + m)
            
            harmonics[lpm] = rhol*p*eim
            harmonics[lmm] = conj(harmonics[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        eim *= ei
    end

    function harmonics_pullback(h̄)

        #@show h̄
        #y,x = sincos(theta)
        ydx = tan(theta)
        xdy = cot(theta)
        T = eltype(harmonics)
        r̄ho = zero(T)
        t̄heta = zero(T)
        p̄hi = zero(T)
        H̄conjH = zero(T)
        _pm2 = zero(eltype(theta))
        _pm1 = zero(eltype(theta))
        _pm0 = zero(eltype(theta)) # used in updating _pm1 and _pm2

        for m=0:P
            lpl = m * m + 2 * m + 1
            # Precomputing conj(H[lpl])*H̄[lpl] since it shows up everywhere.
            H̄conjH = conj(h̄[lpl])*harmonics[lpl]
            # r̄ho[lpl] = (m/rho)*conj(H[lpl])*H̄[lpl]
            r̄ho += 2*m/rho*H̄conjH

            # t̄heta[lpl] = (-tan(theta) + m*cot(theta))*conj(H[lpl])*H̄[lpl]
            t̄heta += 2*(ydx + m*xdy)*H̄conjH

            # p̄hi[lpl] = -i*m*conj(H[lpl])*H̄[lpl]
            p̄hi += -im*m*H̄conjH

            # _pl is pre-recorded.
            _pm2 = x*(2*m+1)*pl#_p[m+1]
            _pm1 = (_pm2*x*(4*m+5) - (2*m+2))*pl/3 # offset indexing by 1 because m starts at zero

            for l=m+1:P
                lpm = l * l + l + m + 1
                # again, precompute this because it gets used several times.
                H̄conjH = conj(h̄[lpm])*harmonics[lpm]

                # r̄ho[lpm] = ((P-2)/rho)*conj(H[lpm])*H̄[lpm]
                r̄ho += 2*(P-2)/rho*H̄conjH

                # t̄heta[lpm] = ∂p/∂θ * 1/p * conj(H[lpm])*H̄[lpm]      
                # ∂p/∂θ * 1/p = (l - m + 1)/sin(θ) - (l+1)*(_pm1/_pm0)*cot(θ)
                _pm0 = (x * (2 * l + 1) * _pm1 - (l + m) * _pm2) / (l - m + 1)
                _pm2 = _pm1
                _pm1 = _pm0
                t̄heta += 2*((l - m + 1) / y - (l + 1) * (_pm1 / _pm0) * xdy) * H̄conjH

                #p̄hi[lpm] = -i*m*conj(H[lpm])*H̄[lpm]
                p̄hi += -2*m*H̄conjH

            end
            pl = -1*(2*m+1)*y*pl
        end
        h̄armonics = zero(T) # harmonics is completely overwritten.
        s̄elf = NoTangent() # not a closure
        P̄ = NoTangent() # P is the multipole expansion order
        # The sums result from the fact that the repeated conj(H[lpm])*H̄[lpm] term is actually a dot product... I think. I'll have to do more work on this.
        return s̄elf, h̄armonics, r̄ho, t̄heta, p̄hi, P̄
    end
    return harmonics, harmonics_pullback

end
#ReverseDiff.@grad_from_chainrules regular_harmonic!(harmonics, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)
ReverseDiff.@grad_from_chainrules regular_harmonic!(harmonics::Vector{Complex{ReverseDiff.TrackedReal}}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)

"""
function spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
"""

function get_drjdxi(theta,phi,rho)
    s_theta,c_theta = sincos(theta)
    s_phi,c_phi = sincos(phi)
    return @SMatrix [s_theta*c_phi    c_theta*c_phi/rho      -s_phi/rho/s_theta
            s_theta * s_phi  c_theta * s_phi / rho  c_phi / rho / s_theta
            c_theta         -s_theta / rho          0                    ]
end

function get_drkdxidxj(theta,phi,rho,k_coord)
    s_theta,c_theta = sincos(theta)
    s_phi,c_phi = sincos(phi)

    if k_coord == 1 # r coordinate
        return @SMatrix [
            (1-c_phi^2 * s_theta^2)/rho -s_theta^2*c_phi*s_phi/rho -s_theta*c_phi*c_theta/rho;
            (-s_theta^2*c_phi*s_phi)/rho (1-s_theta^2*s_phi^2)/rho -s_theta*s_phi*c_theta/rho;
            -s_theta*c_phi*c_theta/rho -s_theta*s_phi*c_theta/rho s_theta^2/rho
        ]
    elseif k_coord == 2 # theta coordinate
        return @SMatrix [
            c_theta/s_theta*(1-c_phi^2*(1+2*s_theta^2))/rho^2 -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_phi*(1-2*c_theta^2)/rho^2;
            -c_theta/s_theta*s_phi*c_phi*(1+2*s_theta^2)/rho^2 c_theta/s_theta*(1-s_phi^2*(1+2*s_theta^2))/rho^2 (2*s_theta^2-1)/rho^2*s_phi;
            c_phi*(1-2*c_theta^2)/rho^2 (2*s_theta^2-1)/rho^2*s_phi 2*s_theta*c_theta/rho^2
        ]
    else # phi coordinate
        return @SMatrix [
            2*c_phi*s_phi/rho^2/s_theta^2 (2*s_phi^2-1)/rho^2/s_theta^2 0;
            (2*s_phi^2-1)/rho^2/s_theta^2 -2*s_phi*c_phi/rho^2/s_theta^2 0;
            0 0 0
        ]
    end

end

function ChainRulesCore.rrule(::typeof(s2c_hess!),potential_jacobian, potential_hessian, workspace, rho, theta, phi)

    potential_hessian_out = s2c_hess!(potential_jacobian,potential_hessian,workspace,rho,theta,phi)
    T(A) = permutedims(A,(3,1,2))

    function H_pullback(H̄2)
        
        R = get_drjdxi(theta,phi,rho)
        # partial derivatives of drjdxi ≡ R:
        dRdr = ForwardDiff.derivative((_rho)->get_drjdxi(theta,phi,_rho),rho)
        dRdt = ForwardDiff.derivative((_theta)->get_drjdxi(_theta,phi,rho),theta)
        dRdp = ForwardDiff.derivative((_phi)->get_drjdxi(theta,_phi,rho),phi)

        # partial derivatives of drkdxidxj ≡ [Rr; Rj; Rp]
        dRrdr = ForwardDiff.derivative((_rho)->get_drkdxidxj(theta,phi,_rho,1),rho)
        dRtdt = ForwardDiff.derivative((_theta)->get_drkdxidxj(_theta,phi,rho,2),theta)
        dRpdp = ForwardDiff.derivative((_phi)->get_drkdxidxj(theta,_phi,rho,3),phi)
        r̄ho = zero(eltype(rho))
        t̄heta = zero(eltype(theta))
        p̄hi = zero(eltype(phi))
        H̄ = zeros(size(potential_hessian))
        J̄ = zeros(size(potential_jacobian))
        #@show size(dRdr) size(potential_hessian) size(R) size(dRrdr) size(potential_jacobian) size(H̄2)
        for ind = 1:4
            r̄ho += sum(dRdr'*H̄2[:,:,ind]*R*potential_hessian[:,:,ind]' + dRdr*potential_hessian[:,:,ind]'*R'*H̄2[:,:,ind] + dRrdr'*H̄2[:,:,ind]*potential_jacobian[1,ind]')
            t̄heta += sum(dRdt'*H̄2[:,:,ind]*R*potential_hessian[:,:,ind]' + dRdt*potential_hessian[:,:,ind]'*R'*H̄2[:,:,ind] + dRtdt'*H̄2[:,:,ind]*potential_jacobian[2,ind]')
            p̄hi += sum(dRdp'*H̄2[:,:,ind]*R*potential_hessian[:,:,ind]' + dRdp*potential_hessian[:,:,ind]'*R'*H̄2[:,:,ind] + dRpdp'*H̄2[:,:,ind]*potential_jacobian[3,ind]')
            #@show J̄ get_drkdxidxj(theta,phi,rho,ind)'*H̄2[:,:,ind]
            J̄[:,ind] .+= sum(get_drkdxidxj(theta,phi,rho,ind)'*H̄2[:,:,ind]) # hopefully this line is correct -_-
            H̄[:,:,ind] .+= R'*H̄2[:,:,ind]*R
        end
        #H̄ = R'*H̄2*R
        W̄ = zeros(size(workspace))'*H̄2[:,:,1] # just memory for storing calculations in. May need to be ZeroTangent() instead.
        #@show size(J̄) size(H̄) size(W̄) size(r̄ho) size(t̄heta) size(p̄hi)
        s̄elf = NoTangent() # not a closure
        
        return s̄elf, J̄, H̄, W̄, r̄ho, t̄heta, p̄hi
    end
    return potential_hessian_out, H_pullback

end
ReverseDiff.@grad_from_chainrules s2c_hess!(potential_jacobian, potential_hessian, workspace, rho::tracked_type, theta::tracked_type, phi::tracked_type)

function ChainRulesCore.rrule(::typeof(s2c_jac!),potential_jacobian, workspace, rho, theta, phi)

    potential_jacobian_out = s2c_jac!(copy(potential_jacobian),workspace,rho,theta,phi)

    function J_pullback(J̄2)

        R = get_drjdxi(theta,phi,rho)
        # partial derivatives of drjdxi ≡ R:
        dRdr = ForwardDiff.derivative((_rho)->get_drjdxi(theta,phi,_rho),rho)
        dRdt = ForwardDiff.derivative((_theta)->get_drjdxi(_theta,phi,rho),theta)
        dRdp = ForwardDiff.derivative((_phi)->get_drjdxi(theta,_phi,rho),phi)

        r̄ho = dRdr'*J̄2
        t̄heta = dRdt'*J̄2
        p̄hi = dRdp'*J̄2
        J̄ = R'*J̄2
        W̄ = zeros(size(workspace))'*J̄2 # just memory for storing calculations in.
        s̄elf = NoTangent() # not a closure
        
        return s̄elf, J̄, W̄, sum(r̄ho), sum(t̄heta), sum(p̄hi)

    end
    return potential_jacobian_out, J_pullback

end
ReverseDiff.@grad_from_chainrules s2c_jac!(potential_jacobian, workspace, rho::tracked_type, theta::tracked_type, phi::tracked_type)

function ChainRulesCore.rrule(::typeof(flatten_jacobian!),jacobian)

    j = flatten_jacobian!(copy(jacobian))
    
    function j_pullback(J̄)

        # I need to write real pullbacks for this function and the next one
        C = zeros(size(jacobian))
        C[1:3,1] .= -1.0
        C[2,4] = 1.0
        C[3,2] = 1.0
        C[1,3] = 1.0
        C[3,3] = -1.0
        C[1,4] = -1.0
        C[2,2] = -1.0
        J̄_out = J̄.*C
        #@show size(J̄_out)
        return NoTangent(),J̄_out

    end
    return j,j_pullback

end
ReverseDiff.@grad_from_chainrules flatten_jacobian!(jacobian::tracked_type)

# specialized for the 3x3x4 hessian.
function ChainRulesCore.rrule(::typeof(flatten_hessian!),hessian)

    h = flatten_hessian!(copy(hessian))
    
    function h_pullback(H̄)

        C = zeros(size(hessian))
        C[1:3,1:3,1] .= -1.0
        C[2,1:3,4] .= 1.0
        C[3,1:3,2] .= 1.0
        C[1,1:3,3] .= 1.0
        C[3,1:3,3] .= -1.0
        C[1,1:3,4] .= -1.0
        C[2,1:3,2] .= -1.0

        H̄_out = H̄.*C
        #@show size(H̄_out)
        return NoTangent(), H̄_out#permutedims(H̄_out,(2,1,3))

    end
    return h,h_pullback

end
ReverseDiff.@grad_from_chainrules flatten_hessian!(hessian::tracked_type)

function ChainRulesCore.rrule(::typeof(update_scalar_potential!),scalar_potential,LE,h,P)
    error("the pullback probably needs to be recalculated")
    scalar_potential2 = update_scalar_potential!(copy(scalar_potential),LE,h,P)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = zeros(size(scalar_potential))
        L̄E = zeros(size(LE))
        h̄ = zeros(size(h))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            L̄E[nms] += real(p̄[]*h[nms])
            h̄[nms] += real(p̄[]*LE[nms])
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                L̄E[nms] += 2*real(p̄[]*h[nms])
                h̄[nms] += 2*real(p̄[]*LE[nms])
                
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄

    end
    return scalar_potential2,potential_pullback

end
# Both of these are needed for the same function call due to the SubArray input.
# why does this cause a segfault unless loaded last
ReverseDiff.@grad_from_chainrules update_scalar_potential!(scalar_potential::Array{<:ReverseDiff.TrackedReal, 0},
                                                           LE::Matrix{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           h::Vector{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           P)
ReverseDiff.@grad_from_chainrules update_scalar_potential!(scalar_potential::SubArray{<:ReverseDiff.TrackedReal, 0, Vector{<:ReverseDiff.TrackedReal}, Tuple{Int64}, true},
                                                           LE::Matrix{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           h::Vector{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           P)

function ChainRulesCore.rrule(::typeof(update_vector_potential!),vector_potential,LE,h,P)

    vector_potential2 = update_vector_potential!(copy(vector_potential),LE,h,P)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = zeros(size(vector_potential))
        L̄E = zeros(size(LE))
        h̄ = zeros(size(h))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            L̄E[nms,1] += real(p̄[1]*h[nms])
            L̄E[nms,2] += real(p̄[2]*h[nms])
            L̄E[nms,3] += real(p̄[3]*h[nms])
            h̄[nms] += real(p̄[1]*LE[nms,1])
            h̄[nms] += real(p̄[2]*LE[nms,2])
            h̄[nms] += real(p̄[3]*LE[nms,3])
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                L̄E[nms,1] += 2*real(p̄[1]*h[nms])
                L̄E[nms,2] += 2*real(p̄[2]*h[nms])
                L̄E[nms,3] += 2*real(p̄[3]*h[nms])
                h̄[nms] += 2*real(p̄[1]*LE[nms,1])
                h̄[nms] += 2*real(p̄[2]*LE[nms,2])
                h̄[nms] += 2*real(p̄[3]*LE[nms,3])
                
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄

    end
    return vector_potential2,potential_pullback

end
# Both of these are needed for the same function call due to the SubArray input.
#=
ReverseDiff.@grad_from_chainrules update_vector_potential!(vector_potential::Array{ReverseDiff.TrackedReal{Float64, Float64, Nothing}, 1},
                                                           LE::Matrix{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           h::Vector{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           P)
ReverseDiff.@grad_from_chainrules update_vector_potential!(vector_potential::SubArray{ReverseDiff.TrackedReal{Float64, Float64, Nothing}, 1, Vector{ReverseDiff.TrackedReal{Float64, Float64, Nothing}}, Tuple{UnitRange{Int64}}, true},
                                                           LE::Matrix{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           h::Vector{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           P)
                                                            =#

function ChainRulesCore.rrule(::typeof(update_potential_jacobian!),vector_potential,LE,h,ht,P,r)

    potential_jacobian2 = update_potential_jacobian!(copy(potential_jacobian),LE,h,P,r)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = zeros(size(potential_jacobian))
        L̄E = zeros(size(LE))
        h̄ = zeros(size(h))
        h̄t = zeros(size(ht))
        r̄ = zero(eltype(r))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            for ind in 1:4
                # contribution from J[1,ind]
                r̄ += -potential_jacobian[nms,ind]/r*p̄[nms,ind]
                L̄E[nms,ind] += n/r*(real(h[nms]))*p̄[nms,ind]
                h̄[nms] += n/r*(real(LE[nms,ind]))*p̄[nms,ind]
                # contribution from J[2,ind]
                L̄E[nms,ind] += real(ht[nms])*p̄[nms,ind]
                h̄t[nms,ind] += real(LE[nms,ind])*p̄[nms,ind]
            end
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                for ind in 1:4
                    # contribution from J[1,ind]
                    L̄E[nms,ind] += 2*n/r*(real(h[nms]))*p̄[nms,ind]
                    h̄[nms] += 2*n/r*(real(LE[nms,ind]))*p̄[nms,ind]
                    # contribution from J[2,ind]
                    L̄E[nms,ind] += 2*real(ht[nms])*p̄[nms,ind]
                    h̄t[nms,ind] += 2*real(LE[nms,ind])*p̄[nms,ind]
                    # contribution from J[3,ind]
                    L̄E[nms,ind] += 2*m*imag(h[nms])*p̄[nms,ind] # sign flip due to complex conjugate in inner product
                    h̄[nms] += 2*m*imag(LE[nms,ind])*p̄[nms,ind] # sign flip due to complex conjugate in inner product

                end
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, h̄t, P̄, r̄

    end
    return potential_jacobian2,potential_pullback

end
#=
# Both of these are needed for the same function call due to the SubArray input.
ReverseDiff.@grad_from_chainrules update_potential_jacobian!(vector_potential::Array{ReverseDiff.TrackedReal{Float64, Float64, Nothing}, 1},
                                                           LE::Matrix{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           h::Vector{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           P,
                                                           r<:ReverseDiff.TrackedReal)
ReverseDiff.@grad_from_chainrules update_potential_jacobian!(vector_potential::SubArray{ReverseDiff.TrackedReal{Float64, Float64, Nothing}, 1, Vector{ReverseDiff.TrackedReal{Float64, Float64, Nothing}}, Tuple{UnitRange{Int64}}, true},
                                                           LE::Matrix{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           h::Vector{Complex{ReverseDiff.TrackedReal{Float64, Float64, ReverseDiff.TrackedArray{Float64, Float64, 1, Vector{Float64}, Vector{Float64}}}}},
                                                           P,
                                                           r<:ReverseDiff.TrackedReal)=#