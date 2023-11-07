# Defines rrules for ChainRules and registers them with ReverseDiff.
tracked_type = Union{ReverseDiff.TrackedReal,ReverseDiff.TrackedArray} # tracked type used for registering ChainRulesCore rrules with ReverseDiff.

"""
function regular_harmonic!(harmonics, rho, theta, phi, P)
"""
# wrong results
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
        r̄ho = zero(eltype(rho))
        t̄heta = zero(eltype(theta))
        p̄hi = zero(eltype(phi))
        H̄conjH = zero(eltype(_harmonics))
        _pm2 = zero(eltype(theta))
        _pm1 = zero(eltype(theta))
        _pm0 = zero(eltype(theta)) # used in updating _pm1 and _pm2

        for m=0:P
            lpl = m * m + 2 * m + 1
            # Precomputing conj(H[lpl])*H̄[lpl] since it shows up everywhere.
            H̄conjH = conj(h̄[lpl])*_harmonics[lpl]
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
                H̄conjH = conj(h̄[lpm])*_harmonics[lpm]

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
        h̄armonics = h̄ # harmonics is completely overwritten.
        s̄elf = NoTangent() # not a closure
        P̄ = NoTangent() # P is the multipole expansion order
        return s̄elf, h̄armonics, r̄ho, t̄heta, p̄hi, P̄ # complex values are returned where they should be real...
    end
    return harmonics, harmonics_pullback

end
#ReverseDiff.@grad_from_chainrules regular_harmonic!(harmonics, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)
#ReverseDiff.@grad_from_chainrules regular_harmonic!(harmonics::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)

"""
function spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, rho, theta, phi)
"""

function get_drjdxi(rho,theta,phi)
    s_theta,c_theta = sincos(theta)
    s_phi,c_phi = sincos(phi)
    return @SMatrix [s_theta*c_phi    c_theta*c_phi/rho      -s_phi/rho/s_theta
            s_theta * s_phi  c_theta * s_phi / rho  c_phi / rho / s_theta
            c_theta         -s_theta / rho          0                    ]
end

function get_drkdxidxj(rho,theta,phi,k_coord)
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

function get_RR(rho,theta,phi)

    s_theta,c_theta = sincos(theta)
    s_phi,c_phi = sincos(phi)
    return cat(get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,1),
               get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,2),
               get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,3),
               dims=3)

end

# finally passes tests.
function ChainRulesCore.rrule(::typeof(s2c_hess!),potential_jacobian, potential_hessian, workspace, rho, theta, phi)

    s_theta, c_theta = sincos(theta)
    s_phi, c_phi = sincos(phi)

    R = get_drjdxi(rho,s_theta,s_phi,c_theta,c_phi)
    # we only need to evaluate these functions once instead of 4 times
    RR = get_RR(rho,theta,phi)
    #Rr = get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,1)
    #Rt = get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,2)
    #Rp = get_drkdxidxj(rho,s_theta,s_phi,c_theta,c_phi,3)
    #hessian_out = copy(potential_hessian)
    hessian_out = zeros(eltype(potential_hessian), size(potential_hessian))

    # convert Hessian to cartesian coordinates
    workspace3x3 = view(workspace,:,1:3)
    for ind in 1:4
        workspace3x3 .= potential_hessian[:,:,ind]
        hessian_out[:,:,ind] .= R * workspace3x3 * transpose(R)
        for ν in 1:3
            hessian_out[:,:,ind] .+= RR[:,:,ν] * potential_jacobian[ν,ind]
        end
    end
    #temp = 0.0
    function H_pullback(H̄2)
        
        #R = get_drjdxi(rho,theta,phi)
        # partial derivatives of drjdxi ≡ R:
        dRdxν = cat(ForwardDiff.derivative((_rho)->get_drjdxi(_rho,theta,phi),rho),
                    ForwardDiff.derivative((_theta)->get_drjdxi(rho,_theta,phi),theta),
                    ForwardDiff.derivative((_phi)->get_drjdxi(rho,theta,_phi),phi),
                    dims = 3)
        
        # partial derivatives of drkdxidxj
        dRRσdxν = cat(ForwardDiff.derivative((_rho)->get_RR(_rho,theta,phi),rho),
                     ForwardDiff.derivative((_theta)->get_RR(rho,_theta,phi),theta),
                     ForwardDiff.derivative((_phi)->get_RR(rho,theta,_phi),phi),
                     dims=4)
        x̄ = zeros(eltype(rho),3) # rho, theta, and phi cotangents
        H̄ = zeros(eltype(potential_hessian),size(potential_hessian))
        J̄ = zeros(eltype(potential_jacobian),size(potential_jacobian))
        for μ = 1:4 # μ ≡ ind
            for i=1:3
                for j = 1:3
                    for k=1:3
                        for l=1:3
                            #H̄[j,l,μ] += H̄2[i,k,μ]*R[i,j]*R[k,l]
                            #H̄[j,k,μ] += R[j,i]*H̄2[i,l,μ]*R[k,l]
                            H̄[j,k,μ] += R[i,j]*H̄2[i,l,μ]*R[l,k]
                        end
                    end
                end
            end
        end
        for μ = 1:4
            for σ = 1:3
                #J̄[σ,μ] += H̄2[:,:,μ]*RR[:,:,σ]
                for i = 1:3
                    for k = 1:3
                        J̄[σ,μ] += H̄2[i,k,μ]*RR[k,i,σ]
                    end
                end
            end
        end
        for ν = 1:3
            for i=1:3
                for l=1:3
                    for μ=1:4
                        for j=1:3
                            for k = 1:3
                                #x̄[ν] += R[k,l]*potential_hessian[k,j,μ]*dRdxν[j,i,ν]*H̄2[i,l,μ]
                                #x̄[ν] += dRdxν[k,l,ν]*potential_hessian[k,j,μ]*R[j,i]*H̄2[i,l,μ]
                                x̄[ν] += R[l,k]*potential_hessian[k,j,μ]*dRdxν[i,j,ν]*H̄2[l,i,μ]
                                x̄[ν] += dRdxν[l,k,ν]*potential_hessian[k,j,μ]*R[i,j]*H̄2[l,i,μ]
                            end
                        end
                    end
                end
            end
        end
        for ν = 1:3
            for σ = 1:3
                for μ = 1:4
                    for i=1:3
                        for k=1:3
                            x̄[ν] += H̄2[i,k,μ]*dRRσdxν[k,i,σ,ν]*potential_jacobian[σ,μ]
                            #x̄[ν] += H̄2[i,l,μ]*dRRσdxν[l,i,σ,ν]*potential_jacobian[σ,μ]
                        end
                    end
                end
            end
        end
        W̄ = zeros(size(workspace)) # just memory for storing calculations in.
        s̄elf = NoTangent() # not a closure
        
        return s̄elf, J̄, H̄, W̄, x̄[1], x̄[2], x̄[3]
    end
    return hessian_out, H_pullback

end
ReverseDiff.@grad_from_chainrules s2c_hess!(potential_jacobian, potential_hessian, workspace, rho::tracked_type, theta::tracked_type, phi::tracked_type)

# passes tests
function ChainRulesCore.rrule(::typeof(s2c_jac!),potential_jacobian, workspace, rho, theta, phi)
    #println("initial call:")
    #@show potential_jacobian
    function J_pullback(J̄2)
        #println("pullback call:")
        #@show potential_jacobian

        R = get_drjdxi(rho,theta,phi)
        # partial derivatives of drjdxi ≡ R:
        dRdr = ForwardDiff.derivative((_rho)->get_drjdxi(_rho,theta,phi),rho)
        dRdt = ForwardDiff.derivative((_theta)->get_drjdxi(rho,_theta,phi),theta)
        dRdp = ForwardDiff.derivative((_phi)->get_drjdxi(rho,theta,_phi),phi)
        r̄ho = 0.0
        t̄heta = 0.0
        p̄hi = 0.0
        C = 0.0
        # r̄ho' = J2'[j,i]*dRdr[i,k]*J[k,j] (although in the end it's just a scalar so the transpose doesn't matter for r̄ho')
        # and similarly for t̄heta and p̄hi
        for i=1:size(potential_jacobian)[1]
            for j=1:size(potential_jacobian)[2]
                for k=1:size(potential_jacobian)[1]
                    C = J̄2[i,j]*potential_jacobian[k,j]
                    r̄ho += dRdr[i,k]*C
                    t̄heta += dRdt[i,k]*C
                    p̄hi += dRdp[i,k]*C
                end
            end
        end
        J̄ = R'*J̄2
        W̄ = zeros(size(workspace)) # just memory for storing calculations in.
        s̄elf = NoTangent() # not a closure
        #@show J̄ sum(r̄ho) sum(t̄heta) sum(p̄hi)
        
        return s̄elf, J̄, W̄, sum(r̄ho), sum(t̄heta), sum(p̄hi)

    end
    return s2c_jac!(copy(potential_jacobian),workspace,rho,theta,phi), J_pullback
    #potential_jacobian .= s2c_jac!(potential_jacobian,workspace,rho,theta,phi)
    #return potential_jacobian, J_pullback

end
ReverseDiff.@grad_from_chainrules s2c_jac!(potential_jacobian, workspace, rho::tracked_type, theta::tracked_type, phi::tracked_type)

function ChainRulesCore.rrule(::typeof(flatten_jacobian!),jacobian)

    jacobian = flatten_jacobian!(jacobian)
    
    function j_pullback(J̄)

        J̄_out = similar(jacobian)
        J̄_out .= J̄
        # column 1:
        J̄_out[1:3,1] .*= -1
        # column 2:
        #J̄_out[1,2] += 0
        J̄_out[2,2] -= J̄[3,1]
        J̄_out[3,2] += J̄[2,1]
        # column 3:
        J̄_out[1,3] += J̄[3,1]
        #J̄_out[2,3] += 0
        J̄_out[3,3] -= J̄[1,1]
        # column 4
        J̄_out[1,4] -= J̄[2,1]
        J̄_out[2,4] += J̄[1,1]
        #J̄_out[3,4] += 0
        #@show size(J̄_out)
        return NoTangent(),J̄_out

    end
    return jacobian,j_pullback

end
ReverseDiff.@grad_from_chainrules flatten_jacobian!(jacobian::AbstractArray{<:ReverseDiff.TrackedReal})

# not fully checked, but should work as long as there are no typos
function ChainRulesCore.rrule(::typeof(flatten_hessian!),hessian)

    hessian = flatten_hessian!(hessian)
    
    function h_pullback(H̄)

        H̄_out = copy(H̄) # 3x3x4
        H̄_out[1:3,1:3,1] .*= -1

        H̄_out[2,1,4] += H̄[1,1,1]
        H̄_out[3,1,3] -= H̄[1,1,1]
        
        H̄_out[3,1,2] += H̄[2,1,1]
        H̄_out[1,1,4] -= H̄[2,1,1]
        
        H̄_out[1,1,3] += H̄[3,1,1]
        H̄_out[2,1,2] -= H̄[3,1,1]
        
        H̄_out[2,2,4] += H̄[1,2,1]
        H̄_out[3,2,3] -= H̄[1,2,1]
        
        H̄_out[3,2,2] += H̄[2,2,1]
        H̄_out[1,2,4] -= H̄[2,2,1]
        
        H̄_out[1,2,3] += H̄[3,2,1]
        H̄_out[2,2,2] -= H̄[3,2,1]
        
        H̄_out[2,3,4] += H̄[1,3,1]
        H̄_out[3,3,3] -= H̄[1,3,1]
        
        H̄_out[3,3,2] += H̄[2,3,1]
        H̄_out[1,3,4] -= H̄[2,3,1]
        
        H̄_out[1,3,3] += H̄[3,3,1]
        H̄_out[2,3,2] -= H̄[3,3,1]
        
        return NoTangent(), H̄_out

    end
    return hessian,h_pullback

end
ReverseDiff.@grad_from_chainrules flatten_hessian!(hessian::AbstractArray{<:ReverseDiff.TrackedReal})

# works
function ChainRulesCore.rrule(::typeof(update_scalar_potential),scalar_potential,LE,h,P)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            L̄E[1,nms] += p̄*conj(h[nms])
            h̄[nms] += p̄*conj(LE[1,nms])
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                L̄E[1,nms] += 2*p̄*conj(h[nms])
                h̄[nms] += 2*p̄*conj(LE[1,nms])
                
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄

    end
    return update_scalar_potential(copy(scalar_potential),LE,h,P),potential_pullback

end

ReverseDiff.@grad_from_chainrules update_scalar_potential!(scalar_potential::ReverseDiff.TrackedReal,
                                                           LE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           P)
                                                           

function ChainRulesCore.rrule(::typeof(update_vector_potential!),vector_potential,LE,h,P)

    vector_potential2 = update_vector_potential!(copy(vector_potential),LE,h,P)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            L̄E[2,nms] += p̄[1]*conj(h[nms])
            L̄E[3,nms] += p̄[2]*conj(h[nms])
            L̄E[4,nms] += p̄[3]*conj(h[nms])
            h̄[nms] += p̄[1]*conj(LE[2,nms])
            h̄[nms] += p̄[2]*conj(LE[3,nms])
            h̄[nms] += p̄[3]*conj(LE[4,nms])
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                L̄E[2,nms] += 2*p̄[1]*conj(h[nms])
                L̄E[3,nms] += 2*p̄[2]*conj(h[nms])
                L̄E[4,nms] += 2*p̄[3]*conj(h[nms])
                h̄[nms] += 2*p̄[1]*conj(LE[2,nms])
                h̄[nms] += 2*p̄[2]*conj(LE[3,nms])
                h̄[nms] += 2*p̄[3]*conj(LE[4,nms])
                
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄

    end
    return vector_potential2,potential_pullback

end

ReverseDiff.@grad_from_chainrules update_vector_potential!(vector_potential::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           LE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           P)
                                                           
                                                          
# works
function ChainRulesCore.rrule(::typeof(update_potential_jacobian!),potential_jacobian,LE,h,ht,P,r)
    potential_jacobian_out = update_potential_jacobian!(copy(potential_jacobian),LE,h,ht,P,r)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        h̄t = zeros(eltype(ht),size(ht))
        r̄ = zero(eltype(r))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        # no need to loop for r̄, since it factors out of the loops.
        #for ind in 1:4
        #    r̄ += -potential_jacobian_out[1,ind]/r*p̄[1,ind]
        #    r̄ += -2*potential_jacobian_out[1,ind]/r*p̄[1,ind]
        #end

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            for ind in 1:4

                # contribution from J[1,ind]
                r̄ += -n/r^2*real(conj(h[nms]*LE[ind,nms]*p̄[1,ind]))
                L̄E[ind,nms] += n/r*conj(h[nms])*p̄[1,ind]
                h̄[nms] += n/r*conj(LE[ind,nms])*p̄[1,ind]
                # contribution from J[2,ind]
                L̄E[ind,nms] += conj(ht[nms])*p̄[2,ind]
                h̄t[nms] += conj(LE[ind,nms])*p̄[2,ind]
            end
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                for ind in 1:4
                    # contribution from J[1,ind]
                    r̄ += -2*n/r^2*real(conj(h[nms]*LE[ind,nms]*p̄[1,ind]))
                    L̄E[ind,nms] += 2*n/r*conj(h[nms])*p̄[1,ind]
                    h̄[nms] += 2*n/r*conj(LE[ind,nms])*p̄[1,ind]
                    # contribution from J[2,ind]
                    L̄E[ind,nms] += 2*conj(ht[nms])*p̄[2,ind]
                    h̄t[nms] += 2*conj(LE[ind,nms])*p̄[2,ind]
                    # contribution from J[3,ind]
                    L̄E[ind,nms] += 2*m*conj(h[nms]*im)*p̄[3,ind]
                    h̄[nms] += 2*m*conj(LE[ind,nms]*im)*p̄[3,ind]

                end
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, h̄t, P̄, r̄

    end
    #@show size(potential_jacobian) size(LE) size(h) size(ht) size(P) size(r)
    return potential_jacobian_out,potential_pullback

end

ReverseDiff.@grad_from_chainrules update_potential_jacobian!(potential_jacobian::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           LE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           ht::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           P,
                                                           r::ReverseDiff.TrackedReal)
                                                           
# passes tests
function ChainRulesCore.rrule(::typeof(update_potential_hessian!),potential_hessian,LE,h,ht,ht2,P,r)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        h̄t = zeros(eltype(ht),size(ht))
        h̄t2 = zeros(eltype(ht2),size(ht2))
        r̄ = zero(eltype(r))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        # this appears to break, but I'm not sure why.
        #=for ind in 1:4
            r̄ += 2*potential_hessian[1,1,ind]*p̄[1,1,ind]
            r̄ += potential_hessian[1,2,ind]*p̄[2,1,ind]
            r̄ += potential_hessian[2,1,ind]*p̄[1,2,ind]
            r̄ += potential_hessian[1,3,ind]*p̄[3,1,ind]
            r̄ += potential_hessian[3,1,ind]*p̄[1,3,ind]
        end
        r̄ *= -1/r=#

        # everything else still needs to be looped over.
        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            #c1 = n/r
            #c2 = c1*(n-1)/r
            for ind in 1:4
                # contribution from H[1,1,ind]
                r̄ += -2*n*(n-1)/r^3*real(h[nms]*LE[ind,nms])*p̄[1,1,ind]
                L̄E[ind,nms] += n*(n-1)/r^2*(real(h[nms]) - im*imag(h[nms]))*p̄[1,1,ind]
                h̄[nms] += n*(n-1)/r^2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[1,1,ind]
                # contribution from H[1,2,ind] and H[2,1,ind]
                r̄ += -n/r^2*real(ht[nms]*LE[ind,nms])*(p̄[1,2,ind] + p̄[2,1,ind])
                L̄E[ind,nms] += n/r*(real(ht[nms]) - im*imag(ht[nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                h̄t[nms] += n/r*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                # contriubiton from H[2,2,ind]
                L̄E[ind,nms] += (real(ht2[nms]) - im*imag(ht2[nms]))*p̄[2,2,ind]
                h̄t2[nms] += (real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[2,2,ind]
                #h̄t2[nms] += real(LE[ind,nms])*p̄[2,2,ind]
            end
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                for ind in 1:4
                    # H[1,1,ind]
                    r̄ += -2*2*n*(n-1)/r^3*real(h[nms]*LE[ind,nms])*p̄[1,1,ind]
                    L̄E[ind,nms] += 2*n*(n-1)/r^2*(real(h[nms]) - im*imag(h[nms]))*p̄[1,1,ind]
                    h̄[nms] += 2*n*(n-1)/r^2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[1,1,ind]
                    # H[1,2,ind] and H[2,1,ind]
                    r̄ += -2*n/r^2*real(ht[nms]*LE[ind,nms])*(p̄[1,2,ind] + p̄[2,1,ind])
                    L̄E[ind,nms] += 2*n/r*(real(ht[nms]) - im*imag(ht[nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                    h̄t[nms] += 2*n/r*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                    # H[2,2,ind]
                    L̄E[ind,nms] += 2*(real(ht2[nms]) - im*imag(ht2[nms]))*p̄[2,2,ind]
                    h̄t2[nms] += 2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[2,2,ind]
                    # H[1,3,ind] and H[3,1,ind]
                    r̄ += 2*m*n/r^2*imag(h[nms]*LE[ind,nms])*(p̄[1,3,ind] + p̄[3,1,ind])
                    L̄E[ind,nms] += -2*m*n/r*(imag(h[nms]) + im*real(h[nms]))*(p̄[1,3,ind] + p̄[3,1,ind])
                    h̄[nms] += -2*m*n/r*(imag(LE[ind,nms]) + im*real(LE[ind,nms]))*(p̄[1,3,ind] + p̄[3,1,ind])
                    # H[2,3,ind] and H[3,2,ind]
                    L̄E[ind,nms] += -2*m*(imag(ht[nms]) + im*real(ht[nms]))*(p̄[2,3,ind] + p̄[3,2,ind])
                    h̄t[nms] += -2*m*(imag(LE[ind,nms]) + im*real(LE[ind,nms]))*(p̄[2,3,ind] + p̄[3,2,ind])
                    # H[3,3,ind]
                    L̄E[ind,nms] += -2*m^2*(real(h[nms]) - im*imag(h[nms]))*p̄[3,3,ind]
                    h̄[nms] += -2*m^2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[3,3,ind]
                end
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, h̄t, h̄t2, P̄, r̄

    end
    return update_potential_hessian!(copy(potential_hessian),LE,h,ht,ht2,P,r),potential_pullback

end

ReverseDiff.@grad_from_chainrules update_potential_hessian!(potential_hessian::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           LE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           ht::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           ht2::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                                           P,
                                                           r::ReverseDiff.TrackedReal)

function ChainRulesCore.rrule(::typeof(M2L_loop!),LE,L,ME,h,expansion_order)

    function LE_pullback(L̄E2)
        s̄elf = NoTangent()
        L̄E = L̄E2
        L̄ = zeros(eltype(L),size(L))
        M̄E = zeros(eltype(ME),size(ME))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent()
        for j in 0:expansion_order
            Cnm = odd_or_even(j)
            for k in 0:j
                jks = (j * (j + 1)) >> 1 + k + 1
                
                for n in 0:expansion_order
                    for m in -n:-1
                        nms = (n * (n+1)) >> 1 - m + 1
                        jnkm = (j + n)^2 + j + n + m - k + 1
                        # jnkm_max = (P + P)^2 + P + P + -1 - 0 + 1 = (2P)^2 + 2P = 2P(2P+1)
                        for dim in 1:4
                            M̄E[dim,nms] += conj(L̄E2[dim,jks])*Cnm*h[jnkm]
                            h̄[jnkm] += L̄E2[dim,jks]*Cnm*ME[dim,nms]
                        end
                    end
                    for m in 0:n
                        nms = (n * (n+1)) >> 1 + m + 1
                        jnkm = (j + n) * (j + n) + j + n + m - k + 1
                        # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                        Cnm2 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m)
                        for dim in 1:4
                            M̄E[dim,nms] += L̄E2[dim,jks]*Cnm2*conj(h[jnkm])
                            h̄[jnkm] += L̄E2[dim,jks]*Cnm2*conj(ME[dim,nms])
                        end
                    end
                end
            end
        end
        return s̄elf,L̄E,L̄,M̄E,h̄,P̄
    end
    return M2L_loop!(copy(LE),L,ME,h,expansion_order), LE_pullback
end

ReverseDiff.@grad_from_chainrules M2L_loop!(LE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            L::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            ME::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            expansion_order)

function ChainRulesCore.rrule(::typeof(M2M_loop!),BM,CM,h,P)

    function BM_pullback(B̄M2)
        s̄elf = NoTangent()
        B̄M = B̄M2
        C̄M = zeros(eltype(CM),size(CM))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent()

        #M = zeros(eltype(BM), 4)
        for j in 0:P # iterate over new multipole coefficients B_j^k
            for k in 0:j
                i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
                for l in 0:j
                    for m in max(-l,-j+k+l):min(k-1,l)
                        jlkms = (((j-l) * (j-l+1)) >> 1) + k - m + 1
                        lm = l * l + l - m + 1
                        ipow = ipow2l(m)
                        oddeven = odd_or_even(l)
                        C = ipow * oddeven
                        for dim in 1:4
                            C̄M[dim,jlkms] += B̄M2[dim,i_jk]*C*conj(h[lm])
                            h̄[lm] += B̄M2[dim,i_jk]*C*conj(CM[dim,jlkms])
                        end
                    end
                    for m in k:min(l,j+k-l)
                        jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                        lm = l * l + l - m + 1
                        oddeven = odd_or_even(k + l + m)
                        for dim in 1:4
                            C̄M[dim,jlkms] += conj(B̄M2[dim,i_jk])*h[lm]*oddeven
                            h̄[lm] += B̄M2[dim,i_jk]*CM[dim,jlkms]*oddeven
                        end
                    end
                end
            end
        end
        return s̄elf, B̄M, C̄M, h̄, P̄
        
    end
    return M2M_loop!(copy(BM),CM,h,P),BM_pullback

end
ReverseDiff.@grad_from_chainrules M2M_loop!(BM::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            CM::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            P)

# incorrect results
function ChainRulesCore.rrule(::typeof(L2L_loop!),CLE,BLE,h,L,P)

    function CLE_pullback(C̄LE2)
        s̄elf = NoTangent()
        C̄LE = C̄LE2
        B̄LE = zeros(eltype(BLE),size(BLE))
        h̄ = zeros(eltype(h),size(h))
        L̄ = zeros(eltype(L),size(L))
        P̄ = NoTangent()
        for j in 0:P
            for k in 0:j
                jks = (j * (j + 1)) >> 1 + k + 1
                #L = zeros(eltype(CLE[1]), 4)
                for n in j:P
                    for m in j+k-n:-1
                        jnkm = (n-j) * (n-j) + n - j + m - k + 1
                        nms = (n * (n + 1)) >> 1 - m + 1
                        oddeven = odd_or_even(k)
                        for dim in 1:4
                            #L[dim] += conj(BLE[dim,nms]) * h[jnkm] * oddeven
                            B̄LE[dim,nms] += conj(C̄LE2[dim,jks])*h[jnkm]*oddeven
                            h̄[jnkm] += C̄LE2[dim,jks]*BLE[dim,nms]*oddeven
                        end
                    end
                    for m in 0:n
                        if n-j >= abs(m-k)
                            jnkm = (n - j) * (n - j) + n - j + m - k + 1
                            nms = (n * (n + 1)) >> 1 + m + 1
                            oddeven = odd_or_even((m-k) * (1 >> (m >= k)))
                            for dim in 1:4
                                #L[dim] += BLE[dim,nms] * h[jnkm] * oddeven
                                B̄LE[dim,nms] += C̄LE2[dim,jks]*conj(h[jnkm])*oddeven
                                h̄[jnkm] += C̄LE2[dim,jks]*conj(BLE[dim,nms])*oddeven
                                #M̄E[dim,nms] += L̄E2[dim,jks]*Cnm2*conj(h[jnkm])
                                #h̄[jnkm] += L̄E2[dim,jks]*Cnm2*conj(ME[dim,nms])
                            end
                        end
                    end
                end
                #CLE[:,jks] .+= L
            end
        end
        return s̄elf, C̄LE, B̄LE, h̄, L̄, P̄

    end
    return L2L_loop!(copy(CLE),BLE,h,L,P),CLE_pullback
    
end
ReverseDiff.@grad_from_chainrules L2L_loop!(CLE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            BLE::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            L::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            P)