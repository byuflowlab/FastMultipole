# Defines rrules for ChainRules and registers them with ReverseDiff.
tracked_type = Union{ReverseDiff.TrackedReal,ReverseDiff.TrackedArray} # tracked type used for registering ChainRulesCore rrules with ReverseDiff.

"""
function regular_harmonic!(harmonics, rho, theta, phi, P)
"""

function ChainRulesCore.rrule(::typeof(regular_harmonic!), harmonics, rho, theta, phi, P)
    regular_harmonic!(harmonics,rho,theta,phi,P)
    function harmonics_pullback(h̄)

        y,x = sincos(theta)
        xdy = cot(theta)
        r̄ho = zero(eltype(rho))
        t̄heta = zero(eltype(theta))
        p̄hi = zero(eltype(phi))
        plp1 = 1.0
        pl = 1.0
        plm1 = 1.0
        pm = 1.0
        h̄conjh_lpl = zeros(eltype(harmonics),2)
        h̄conjh_lml = zeros(eltype(harmonics),2)
        h̄conjh_lpm = zeros(eltype(harmonics),2)
        h̄conjh_lmm = zeros(eltype(harmonics),2)

        for m=0:P
            lpl = m * m + 2 * m + 1
            lml = m * m + 1
            # save a few multiplications by only computing these products once
            #h̄conjh_lpl = h̄[lpl]*conj(harmonics[lpl])
            #h̄conjh_lml = h̄[lml]*conj(harmonics[lml])
            h̄conjh_lpl = [h̄[1,lpl]*harmonics[1,lpl] + h̄[2,lpl]*harmonics[2,lpl], -h̄[1,lpl]*harmonics[2,lpl] + h̄[2,lpl]*harmonics[1,lpl]]
            h̄conjh_lml = [h̄[1,lml]*harmonics[1,lml] + h̄[2,lml]*harmonics[2,lml], -h̄[1,lml]*harmonics[2,lml] + h̄[2,lml]*harmonics[1,lml]]
            r̄ho += h̄conjh_lpl[1]*(m/rho)
            r̄ho += h̄conjh_lml[1]*(m/rho)
            t̄heta += h̄conjh_lpl[1]*(m*xdy)
            t̄heta += h̄conjh_lml[1]*(m*xdy)
            #p̄hi += -h̄conjh_lpl*im*m
            #p̄hi += h̄conjh_lml*im*m
            p̄hi += h̄conjh_lpl[2]*m
            p̄hi += -h̄conjh_lml[2]*m

            # initialize these values
            plm1 = pm
            pl = x*(2*m+1)*pm
            for l=m+1:P
                lpm = l * l + l + m + 1
                lmm = l * l + l - m + 1
                
                #h̄conjh_lpm = h̄[lpm]*conj(harmonics[lpm])
                #h̄conjh_lmm = h̄[lmm]*conj(harmonics[lmm])
                h̄conjh_lpm = [h̄[1,lpm]*harmonics[1,lpm] + h̄[2,lpm]*harmonics[2,lpm], -h̄[1,lpm]*harmonics[2,lpm] + h̄[2,lpm]*harmonics[1,lpm]]
                h̄conjh_lmm = [h̄[1,lmm]*harmonics[1,lmm] + h̄[2,lmm]*harmonics[2,lmm], -h̄[1,lmm]*harmonics[2,lmm] + h̄[2,lmm]*harmonics[1,lmm]]

                r̄ho += h̄conjh_lpm[1]*(l/rho)
                t̄heta += h̄conjh_lpm[1]*(l*xdy - (l+m)*plm1/(y*pl))
                #p̄hi += -h̄conjh_lpm*im*m
                p̄hi += m*h̄conjh_lpm[2]
                if lpm !== lmm # this avoids double-counting rho, phi, and theta contributions when lpm == lmm.
                    r̄ho += h̄conjh_lmm[1]*(l/rho)
                    t̄heta += h̄conjh_lmm[1]*(l*xdy - (l+m)*plm1/(y*pl))
                    #p̄hi += h̄conjh_lmm*im*m
                    p̄hi += -m*h̄conjh_lmm[2]

                end
                plp1 = (x*(2*l+1)*pl - (l+m)*plm1) / (l-m+1)
                plm1 = pl
                pl = plp1

            end
            pm *= -1*y*(2*m+1)
        end
        h̄armonics = zeros(eltype(harmonics), size(harmonics)) # harmonics is completely overwritten.
        s̄elf = NoTangent() # not a closure
        P̄ = NoTangent() # P is the multipole expansion order
        return s̄elf, h̄armonics, r̄ho, t̄heta, p̄hi, P̄
    end
    return harmonics, harmonics_pullback

end
ReverseDiff.@grad_from_chainrules regular_harmonic!(harmonics::AbstractArray{<:ReverseDiff.TrackedReal}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)
#=
# all correct except for an incredibly painful third derivative with respect to θ. also accounts for ~15% of tape entries, so it will need to be done eventually.
# update: I think this works? It produces the right answer in practice but fails tests.

function ChainRulesCore.rrule(::typeof(regular_harmonic_derivs!),harmonics,rho,theta,phi,P)
    #@show typeof(harmonics)
    #harmonics = regular_harmonic_derivs!(harmonics,rho,theta,phi,P)
    #regular_harmonic_derivs!(harmonics,rho,theta,phi,P)
    
    function h_pullback(h̄)

        # h̄[1,:] -> harmonics
        # h̄[2,:] -> dhdθ
        # h̄[3,:] -> d2hdθ2
        h̄armonics = zeros(eltype(harmonics), size(harmonics)) # completely overwritten -> just set to zeros.
        h̄armonics .= h̄ # for testing purposes
        r̄ho = zero(eltype(rho))
        t̄heta = zero(eltype(theta))
        p̄hi = zero(eltype(phi))
        
        y,x = sincos(theta)
        x2 = x^2
        x3 = x^3
        xdy = cot(theta)
        #xdy2 = xdy^2
        invY = y == 0 ? 0 : 1 / y
        #invY2 = invY^2
        #invY3 = invY^3
        #y2 = y^2
        pl = 1.0
        a = 0.0
        #rhom = 1.0
        #ei = exp(im * phi)
        #eim = 1.0
        for m=0:P
            p = pl
            lpl = (m * (m + 1)) >> 1 + m + 1
            #h̄armonics[1:3,lpl] .= zero(eltype(h̄armonics))
            h̄armonics[1:2,1:3,lpl] .= zero(eltype(h̄armonics)) # explicitly setting this to zero avoids uninitialized memory showing up.# also this line should be redundant since I already initialize the harmonic cotangents to zero.
            # derivatives of harmonics with respect to rho and phi do not depend on the first index
            for σ=1:3
                #r̄ho += real(h̄[σ,lpl]*conj(m/rho*harmonics[σ,lpl]))
                #p̄hi += real(h̄[σ,lpl]*conj(im*m*harmonics[σ,lpl]))
                r̄ho += m/rho*(h̄[1,σ,lpl]*harmonics[1,σ,lpl] + h̄[2,σ,lpl]*harmonics[2,σ,lpl])
                p̄hi += m*(-h̄[1,σ,lpl]*harmonics[2,σ,lpl] + h̄[2,σ,lpl]*harmonics[1,σ,lpl])#real(h̄[σ,lpl]*conj(im*harmonics[σ,lpl]))
            end
            p1 = p
            p = x * (2 * m + 1) * p1
            #t̄heta += real(h̄[1,lpl]*conj(m*xdy*harmonics[1,lpl]))
            #t̄heta += (-m*invY^2 + m^2*xdy^2)*real(h̄[2,lpl]*conj(harmonics[1,lpl]))
            #t̄heta += (-m^2*xdy + m*(m-1)*(-2*xdy*invY^2 + m*xdy^3))*real(h̄[3,lpl]*conj(harmonics[1,lpl]))
            t̄heta += m*xdy*(h̄[1,1,lpl]*harmonics[1,1,lpl] + h̄[2,1,lpl]*harmonics[2,1,lpl])
            t̄heta += (-m*invY^2 + m^2*xdy^2)*(h̄[1,2,lpl]*harmonics[1,1,lpl] + h̄[2,2,lpl]*harmonics[2,1,lpl])
            t̄heta += (-m^2*xdy + m*(m-1)*(-2*xdy*invY^2 + m*xdy^3))*(h̄[1,3,lpl]*harmonics[1,1,lpl] + h̄[2,3,lpl]*harmonics[2,1,lpl])
            # lpl = m * m + 2 * m + 1
            # lml = m * m + 1

            for l=m+1:P
                lpm = (l * (l + 1)) >> 1 + m + 1
                h̄armonics[1:2,1:3,lpm] .= zero(eltype(h̄armonics))
                # lpm = l * l + l + m + 1
                # lmm = l * l + l - m + 1
                for σ=1:3
                    #r̄ho += real(h̄[σ,lpm]*conj(l/rho*harmonics[σ,lpm]))
                    #p̄hi += real(h̄[σ,lpm]*conj(im*m*harmonics[σ,lpm]))
                    r̄ho += l/rho*(h̄[1,σ,lpm]*harmonics[1,σ,lpm] + h̄[2,σ,lpm]*harmonics[2,σ,lpm])
                    p̄hi += m*(-h̄[1,σ,lpm]*harmonics[2,σ,lpm] + h̄[2,σ,lpm]*harmonics[1,σ,lpm])
                end
                #harmonics[1,lpm] = rhol * p * eim
                # harmonics[lmm] = conj(harmonics[lpm])
                #t̄heta += real(h̄[1,lpm]*conj(harmonics[1,lpm]))*(l*xdy*p - (l+m)*invY*p1)/p
                t̄heta += (h̄[1,1,lpm]*harmonics[1,1,lpm] + h̄[2,1,lpm]*harmonics[2,1,lpm])*(l*xdy*p - (l+m)*invY*p1)/p
                p2 = p1
                p1 = p
                p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
                #harmonics[2,lpm] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim
                #harmonics[3,lpm] = rhol * ((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2 * eim
                #t̄heta += real(h̄[2,lpm]*conj(harmonics[2,lpm]))*(((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2)/(((l - m + 1) * p - (l + 1) * x * p1) * invY) # not the most elegant solution but it works.
                t̄heta += (h̄[1,2,lpm]*harmonics[1,2,lpm] + h̄[2,2,lpm]*harmonics[2,2,lpm])*(((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2)/(((l - m + 1) * p - (l + 1) * x * p1) * invY)
                a = ((8*x + 9*l*x - 5*m^2*x - l^3*x + l*m^2*x - l*x3 + l^3*x3)*p + (3*l + 4*l^2 + 3*l*m - m^2 - 8*x2 + l^3 + l^2*m - l*m^2 - 11*l*x2 - m^3 - 8*m*x2 - 4*l^2*x2 - 3*l*m*x2 - l^3*x2 - l^2*m*x2)*p1)*(invY/(2*x2+1))
                t̄heta += (h̄[1,3,lpm]*harmonics[1,3,lpm] + h̄[2,3,lpm]*harmonics[2,3,lpm])*a/(((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2)

                #t̄heta += real(h̄[3,lpm]*conj(harmonics[3,lpm]))*a/(((m-l-1) * x * p + (m^2 - l*(l+1) + (l+1)^2 * x^2) * p1) * invY^2)
            end
            pl = -pl * (2m+1) * y
    
        end

        return NoTangent(),h̄armonics,r̄ho,t̄heta,p̄hi,NoTangent()

    end
    #@show sum(harmonics)
    return regular_harmonic_derivs!(harmonics,rho,theta,phi,P), h_pullback
    #return harmonics, h_pullback
    #return nothing, h_pullback

end
ReverseDiff.@grad_from_chainrules regular_harmonic_derivs!(harmonics::AbstractArray{<:ReverseDiff.TrackedReal}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)
#@grad_from_chainrules_extended (1,) regular_harmonic_derivs!(harmonics::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)
=#
function ChainRulesCore.rrule(::typeof(irregular_harmonic!),harmonics, rho, theta, phi, P)
    function h_pullback(h̄)
        h̄armonics = zeros(eltype(harmonics),size(harmonics))
        r̄ho = zero(eltype(rho))
        t̄heta = zero(eltype(rho))
        p̄hi = zero(eltype(phi))
        plp1 = 1.0
        pl = 1.0
        plm1 = 1.0
        pm = 1.0
        y,x = sincos(theta)
        xdy = cot(theta)
        h̄conjh_npl = zeros(eltype(harmonics),2)
        h̄conjh_nml = zeros(eltype(harmonics),2)

        for m=0:P
            npl = m * m + 2 * m + 1
            nml = m * m + 1

            #h̄conjh_npl = h̄[npl]*conj(harmonics[npl])
            h̄conjh_npl = [h̄[1,npl]*harmonics[1,npl] + h̄[2,npl]*harmonics[2,npl],-h̄[1,npl]*harmonics[2,npl] + h̄[2,npl]*harmonics[1,npl]]
            #r̄ho += h̄conjh_npl*(-m-1)/rho
            r̄ho += h̄conjh_npl[1]*(-m-1)/rho # we only need the real parts of these in the end.
            #t̄heta += h̄conjh_npl*m*xdy
            t̄heta += h̄conjh_npl[1]*m*xdy
            #p̄hi += -h̄conjh_npl*im*m
            p̄hi += h̄conjh_npl[2]*m
            if npl !== nml # this avoids double-counting rho, phi, and theta contributions when npl == nml.
                #h̄conjh_nml = h̄[nml]*conj(harmonics[nml])
                h̄conjh_nml = [h̄[1,nml]*harmonics[1,nml] + h̄[2,nml]*harmonics[2,nml],-h̄[1,nml]*harmonics[2,nml] + h̄[2,nml]*harmonics[1,nml]]
                #r̄ho += h̄conjh_nml*(-m-1)/rho
                #t̄heta += h̄conjh_nml*m*xdy
                #p̄hi += h̄conjh_nml*im*m
                r̄ho += h̄conjh_nml[1]*(-m-1)/rho
                t̄heta += h̄conjh_nml[1]*m*xdy
                p̄hi += -h̄conjh_nml[2]*m
            end

            plm1 = pm
            pl = x*(2*m+1)*pm

            for l=m+1:P
                npm = l * l + l + m + 1
                nmm = l * l + l - m + 1

                #h̄conjh_npm = h̄[npm]*conj(harmonics[npm])
                #h̄conjh_nmm = h̄[nmm]*conj(harmonics[nmm])
                h̄conjh_npm = [h̄[1,npm]*harmonics[1,npm] + h̄[2,npm]*harmonics[2,npm],-h̄[1,npm]*harmonics[2,npm] + h̄[2,npm]*harmonics[1,npm]]

                #r̄ho += h̄conjh_npm*(-l-1)/rho
                #t̄heta += h̄conjh_npm*(l*xdy - (l+m)*plm1/(y*pl))
                #p̄hi += -h̄conjh_npm*im*m
                r̄ho += h̄conjh_npm[1]*(-l-1)/rho
                t̄heta += h̄conjh_npm[1]*(l*xdy - (l+m)*plm1/(y*pl))
                p̄hi += h̄conjh_npm[2]*m
                
                if npm !== nmm # this avoids double-counting rho, phi, and theta contributions when npm == nmm.
                    h̄conjh_nmm = [h̄[1,nmm]*harmonics[1,nmm] + h̄[2,nmm]*harmonics[2,nmm],-h̄[1,nmm]*harmonics[2,nmm] + h̄[2,nmm]*harmonics[1,nmm]]
                    #r̄ho += h̄conjh_nmm*(-l-1)/rho
                    #t̄heta += h̄conjh_nmm*(l*xdy - (l+m)*plm1/(y*pl))
                    #p̄hi += h̄conjh_nmm*im*m
                    r̄ho += h̄conjh_nmm[1]*(-l-1)/rho
                    t̄heta += h̄conjh_nmm[1]*(l*xdy - (l+m)*plm1/(y*pl))
                    p̄hi += -h̄conjh_nmm[2]*m

                end
                plp1 = (x*(2*l+1)*pl - (l+m)*plm1) / (l-m+1)
                plm1 = pl
                pl = plp1
            end
            pm *= -1*y*(2*m+1)
        end
        
        s̄elf = NoTangent()
        P̄ = NoTangent()
        return s̄elf, h̄armonics, r̄ho, t̄heta, p̄hi, P̄
    end
    return irregular_harmonic!(harmonics, rho, theta, phi, P), h_pullback

end
ReverseDiff.@grad_from_chainrules irregular_harmonic!(harmonics::AbstractArray{<:ReverseDiff.TrackedReal}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)

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
#=
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
=#

#=
# passes tests
function ChainRulesCore.rrule(::typeof(s2c_jac!),potential_jacobian, workspace, rho, theta, phi)\
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
=#

#=
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
    flatten_hessian!(hessian)
    
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
=#

#=
# works
function ChainRulesCore.rrule(::typeof(update_scalar_potential),scalar_potential,LE,h,expansion_order::Val{P}) where P
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            #L̄E[1,nms] += p̄*conj(h[nms])
            #h̄[nms] += p̄*conj(LE[1,nms])
            L̄E[1,1,nms] += p̄*h[1,1,nms]
            L̄E[2,1,nms] -= p̄*h[2,1,nms]
            h̄[1,1,nms] += p̄*LE[1,1,nms]
            h̄[2,1,nms] -= p̄*LE[2,1,nms]
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                #L̄E[1,nms] += 2*p̄*conj(h[nms])
                #h̄[nms] += 2*p̄*conj(LE[1,nms])
                L̄E[1,1,nms] += 2*p̄*h[1,1,nms]
                L̄E[2,1,nms] -= 2*p̄*h[2,1,nms]
                h̄[1,1,nms] += 2*p̄*LE[1,1,nms]
                h̄[2,1,nms] -= 2*p̄*LE[2,1,nms]
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄

    end
    return update_scalar_potential(scalar_potential,LE,h,expansion_order),potential_pullback

end
=#

#=
ReverseDiff.@grad_from_chainrules update_scalar_potential(scalar_potential,
                                                           LE::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           h::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           P::Val)
                                                           

function ChainRulesCore.rrule(::typeof(update_vector_potential!),vector_potential,LE,h,expansion_order::Val{P}) where P

    vector_potential2 = update_vector_potential!(copy(vector_potential),LE,h,expansion_order)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent() # Order of the multipole expansion; non-differentiable integer.

        for n in 0:P
            nms = (n * (n+1)) >> 1 + 1
            #=L̄E[2,nms] += p̄[1]*conj(h[nms])
            L̄E[3,nms] += p̄[2]*conj(h[nms])
            L̄E[4,nms] += p̄[3]*conj(h[nms])
            h̄[nms] += p̄[1]*conj(LE[2,nms])
            h̄[nms] += p̄[2]*conj(LE[3,nms])
            h̄[nms] += p̄[3]*conj(LE[4,nms])=#
            for ν=1:3 # this would be 12 lines fully expanded, but the structure is simple so I stuck it in some loops. ν loops over the entry of p̄ and σ is for real/imaginary parts.
                for σ = 1:2
                    s = (σ == 1 ? 1 : -1)
                    L̄E[σ,ν+1,nms] += p̄[ν]*h[σ,1,nms]*s
                    h̄[σ,1,nms] += p̄[ν]*LE[σ,ν+1,nms]*s
                end
            end

            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                #=L̄E[2,nms] += 2*p̄[1]*conj(h[nms])
                L̄E[3,nms] += 2*p̄[2]*conj(h[nms])
                L̄E[4,nms] += 2*p̄[3]*conj(h[nms])
                h̄[nms] += 2*p̄[1]*conj(LE[2,nms])
                h̄[nms] += 2*p̄[2]*conj(LE[3,nms])
                h̄[nms] += 2*p̄[3]*conj(LE[4,nms])=#
                for ν=1:3 # this would be 12 lines fully expanded, but the structure is simple so I stuck it in some loops. ν loops over the entry of p̄ and σ is for real/imaginary parts.
                    for σ = 1:2
                        s = (σ == 1 ? 2 : -2)
                        L̄E[σ,ν+1,nms] += p̄[ν]*h[σ,1,nms]*s
                        h̄[σ,1,nms] += p̄[ν]*LE[σ,ν+1,nms]*s
                    end
                end
                
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄

    end
    return vector_potential2,potential_pullback

end
=#

#=
ReverseDiff.@grad_from_chainrules update_vector_potential!(vector_potential::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           LE::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           h::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           expansion_order::Val)
                                                           
                                
function ChainRulesCore.rrule(::typeof(update_potential_jacobian!),potential_jacobian,LE,h,expansion_order::Val{P},r) where P
    potential_jacobian_out = update_potential_jacobian!(copy(potential_jacobian),LE,h,expansion_order,r)
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        #h̄t = zeros(eltype(ht),size(ht))
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
                #r̄ += -n/r^2*real(conj(h[nms]*LE[ind,nms]*p̄[1,ind]))
                r̄ += -n/r^2*p̄[1,ind]*(h[1,1,nms]*LE[1,ind,nms] - h[2,1,nms]*LE[2,ind,nms])
                #L̄E[ind,nms] += n/r*conj(h[nms])*p̄[1,ind]
                #h̄[nms] += n/r*conj(LE[ind,nms])*p̄[1,ind]
                L̄E[1,ind,nms] += n/r*h[1,1,nms]*p̄[1,ind]
                L̄E[2,ind,nms] -= n/r*h[2,1,nms]*p̄[1,ind]
                h̄[1,1,nms] += n/r*LE[1,ind,nms]*p̄[1,ind]
                h̄[2,1,nms] -= n/r*LE[2,ind,nms]*p̄[1,ind]

                # contribution from J[2,ind]
                #L̄E[ind,nms] += conj(ht[nms])*p̄[2,ind]
                #h̄t[nms] += conj(LE[ind,nms])*p̄[2,ind]
                L̄E[1,ind,nms] += h[1,2,nms]*p̄[2,ind]
                L̄E[2,ind,nms] -= h[2,2,nms]*p̄[2,ind]
                h̄[1,2,nms] += LE[1,ind,nms]*p̄[2,ind]
                h̄[2,2,nms] -= LE[2,ind,nms]*p̄[2,ind]
            end
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                for ind in 1:4
                    # contribution from J[1,ind]
                    #r̄ += -2*n/r^2*real(conj(h[nms]*LE[ind,nms]*p̄[1,ind]))
                    r̄ += -2*n/r^2*p̄[1,ind]*(h[1,1,nms]*LE[1,ind,nms] - h[2,1,nms]*LE[2,ind,nms])
                    #L̄E[ind,nms] += 2*n/r*conj(h[nms])*p̄[1,ind]
                    #h̄[nms] += 2*n/r*conj(LE[ind,nms])*p̄[1,ind]
                    L̄E[1,ind,nms] += 2*n/r*h[1,1,nms]*p̄[1,ind]
                    L̄E[2,ind,nms] -= 2*n/r*h[2,1,nms]*p̄[1,ind]
                    h̄[1,1,nms] += 2*n/r*LE[1,ind,nms]*p̄[1,ind]
                    h̄[2,1,nms] -= 2*n/r*LE[2,ind,nms]*p̄[1,ind]

                    # contribution from J[2,ind]
                    #L̄E[ind,nms] += 2*conj(ht[nms])*p̄[2,ind]
                    #h̄t[nms] += 2*conj(LE[ind,nms])*p̄[2,ind]
                    L̄E[1,ind,nms] += 2*h[1,2,nms]*p̄[2,ind]
                    L̄E[2,ind,nms] -= 2*h[2,2,nms]*p̄[2,ind]
                    h̄[1,2,nms] += 2*LE[1,ind,nms]*p̄[2,ind]
                    h̄[2,2,nms] -= 2*LE[2,ind,nms]*p̄[2,ind]

                    # contribution from J[3,ind]
                    #L̄E[ind,nms] += 2*m*conj(h[nms]*im)*p̄[3,ind]
                    #h̄[nms] += 2*m*conj(LE[ind,nms]*im)*p̄[3,ind]
                    L̄E[1,ind,nms] -= 2*m*h[2,1,nms]*p̄[3,ind]
                    L̄E[2,ind,nms] -= 2*m*h[1,1,nms]*p̄[3,ind]
                    h̄[1,1,nms] -= 2*m*LE[2,ind,nms]*p̄[3,ind]
                    h̄[2,1,nms] -= 2*m*LE[1,ind,nms]*p̄[3,ind]

                end
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄, r̄

    end
    return potential_jacobian_out,potential_pullback

end
=#

#=
ReverseDiff.@grad_from_chainrules update_potential_jacobian!(potential_jacobian::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           LE::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           h::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           expansion_order::Val,
                                                           r::ReverseDiff.TrackedReal)
                                    
function ChainRulesCore.rrule(::typeof(update_potential_hessian!),potential_hessian,LE,h,expansion_order::Val{P},r) where P
    function potential_pullback(p̄)
        s̄elf = NoTangent() # not a closure
        p̄otential = p̄
        L̄E = zeros(eltype(LE),size(LE))
        h̄ = zeros(eltype(h),size(h))
        #h̄t = zeros(eltype(ht),size(ht))
        #h̄t2 = zeros(eltype(ht2),size(ht2))
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
                #r̄ += -2*n*(n-1)/r^3*real(h[nms]*LE[ind,nms])*p̄[1,1,ind]
                r̄ += -2*n*(n-1)/r^3*(h[1,1,nms]*LE[1,ind,nms] - h[2,1,nms]*LE[2,ind,nms])*p̄[1,1,ind]
                #L̄E[ind,nms] += n*(n-1)/r^2*(real(h[nms]) - im*imag(h[nms]))*p̄[1,1,ind]
                #h̄[nms] += n*(n-1)/r^2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[1,1,ind]
                L̄E[1,ind,nms] += n*(n-1)/r^2*h[1,1,nms]*p̄[1,1,ind]
                L̄E[2,ind,nms] -= n*(n-1)/r^2*h[2,1,nms]*p̄[1,1,ind]
                h̄[1,1,nms] += n*(n-1)/r^2*LE[1,ind,nms]*p̄[1,1,ind]
                h̄[2,1,nms] -= n*(n-1)/r^2*LE[2,ind,nms]*p̄[1,1,ind]

                # contribution from H[1,2,ind] and H[2,1,ind]
                #r̄ += -n/r^2*real(ht[nms]*LE[ind,nms])*(p̄[1,2,ind] + p̄[2,1,ind])
                r̄ += -n/r^2*(h[1,2,nms]*LE[1,ind,nms] - h[2,2,nms]*LE[2,ind,nms])*(p̄[1,2,ind] + p̄[2,1,ind])
                #L̄E[ind,nms] += n/r*(real(ht[nms]) - im*imag(ht[nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                #h̄t[nms] += n/r*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                L̄E[1,ind,nms] += n/r*h[1,2,nms]*(p̄[1,2,ind] + p̄[2,1,ind])
                L̄E[2,ind,nms] -= n/r*h[2,2,nms]*(p̄[1,2,ind] + p̄[2,1,ind])
                h̄[1,2,nms] += n/r*LE[1,ind,nms]*(p̄[1,2,ind] + p̄[2,1,ind])
                h̄[2,2,nms] -= n/r*LE[2,ind,nms]*(p̄[1,2,ind] + p̄[2,1,ind])

                # contribution from H[2,2,ind]
                #L̄E[ind,nms] += (real(ht2[nms]) - im*imag(ht2[nms]))*p̄[2,2,ind]
                #h̄t2[nms] += (real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[2,2,ind]
                L̄E[1,ind,nms] += h[1,3,nms]*p̄[2,2,ind]
                L̄E[2,ind,nms] -= h[2,3,nms]*p̄[2,2,ind]
                h̄[1,3,nms] += LE[1,ind,nms]*p̄[2,2,ind]
                h̄[2,3,nms] -= LE[2,ind,nms]*p̄[2,2,ind]
                #h̄t2[nms] += real(LE[ind,nms])*p̄[2,2,ind]
            end
            for m in 1:n
                nms = (n * (n + 1)) >> 1 + m + 1
                for ind in 1:4
                    # H[1,1,ind]
                    #r̄ += -2*2*n*(n-1)/r^3*real(h[nms]*LE[ind,nms])*p̄[1,1,ind]
                    r̄ += -2*2*n*(n-1)/r^3*(h[1,1,nms]*LE[1,ind,nms] - h[2,1,nms]*LE[2,ind,nms])*p̄[1,1,ind]
                    #L̄E[ind,nms] += 2*n*(n-1)/r^2*(real(h[nms]) - im*imag(h[nms]))*p̄[1,1,ind]
                    #h̄[nms] += 2*n*(n-1)/r^2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[1,1,ind]
                    L̄E[1,ind,nms] += 2*n*(n-1)/r^2*h[1,1,nms]*p̄[1,1,ind]
                    L̄E[2,ind,nms] -= 2*n*(n-1)/r^2*h[2,1,nms]*p̄[1,1,ind]
                    h̄[1,1,nms] += 2*n*(n-1)/r^2*LE[1,ind,nms]*p̄[1,1,ind]
                    h̄[2,1,nms] -= 2*n*(n-1)/r^2*LE[2,ind,nms]*p̄[1,1,ind]

                    # H[1,2,ind] and H[2,1,ind]
                    #r̄ += -2*n/r^2*real(ht[nms]*LE[ind,nms])*(p̄[1,2,ind] + p̄[2,1,ind])
                    r̄ += -2*n/r^2*(h[1,2,nms]*LE[1,ind,nms] - h[2,2,nms]*LE[2,ind,nms])*(p̄[1,2,ind] + p̄[2,1,ind])
                    #L̄E[ind,nms] += 2*n/r*(real(ht[nms]) - im*imag(ht[nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                    #h̄t[nms] += 2*n/r*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*(p̄[1,2,ind] + p̄[2,1,ind])
                    L̄E[1,ind,nms] += 2*n/r*h[1,2,nms]*(p̄[1,2,ind] + p̄[2,1,ind])
                    L̄E[2,ind,nms] -= 2*n/r*h[2,2,nms]*(p̄[1,2,ind] + p̄[2,1,ind])
                    h̄[1,2,nms] += 2*n/r*LE[1,ind,nms]*(p̄[1,2,ind] + p̄[2,1,ind])
                    h̄[2,2,nms] -= 2*n/r*LE[2,ind,nms]*(p̄[1,2,ind] + p̄[2,1,ind])

                    # H[2,2,ind]
                    #L̄E[ind,nms] += 2*(real(ht2[nms]) - im*imag(ht2[nms]))*p̄[2,2,ind]
                    #h̄t2[nms] += 2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[2,2,ind]
                    L̄E[1,ind,nms] += 2*h[1,3,nms]*p̄[2,2,ind]
                    L̄E[2,ind,nms] -= 2*h[2,3,nms]*p̄[2,2,ind]
                    h̄[1,3,nms] += 2*LE[1,ind,nms]*p̄[2,2,ind]
                    h̄[2,3,nms] -= 2*LE[2,ind,nms]*p̄[2,2,ind]

                    # H[1,3,ind] and H[3,1,ind]
                    #r̄ += 2*m*n/r^2*imag(h[nms]*LE[ind,nms])*(p̄[1,3,ind] + p̄[3,1,ind])
                    r̄ += 2*m*n/r^2*(h[1,1,nms]*LE[2,ind,nms] + h[2,1,nms]*LE[1,ind,nms])*(p̄[1,3,ind] + p̄[3,1,ind])
                    #L̄E[ind,nms] += -2*m*n/r*(imag(h[nms]) + im*real(h[nms]))*(p̄[1,3,ind] + p̄[3,1,ind])
                    #h̄[nms] += -2*m*n/r*(imag(LE[ind,nms]) + im*real(LE[ind,nms]))*(p̄[1,3,ind] + p̄[3,1,ind])
                    L̄E[1,ind,nms] += -2*m*n/r*h[2,1,nms]*(p̄[1,3,ind] + p̄[3,1,ind])
                    L̄E[2,ind,nms] += -2*m*n/r*h[1,1,nms]*(p̄[1,3,ind] + p̄[3,1,ind])
                    h̄[1,1,nms] += -2*m*n/r*LE[2,ind,nms]*(p̄[1,3,ind] + p̄[3,1,ind])
                    h̄[2,1,nms] += -2*m*n/r*LE[1,ind,nms]*(p̄[1,3,ind] + p̄[3,1,ind])

                    # H[2,3,ind] and H[3,2,ind]
                    #L̄E[ind,nms] += -2*m*(imag(ht[nms]) + im*real(ht[nms]))*(p̄[2,3,ind] + p̄[3,2,ind])
                    #h̄t[nms] += -2*m*(imag(LE[ind,nms]) + im*real(LE[ind,nms]))*(p̄[2,3,ind] + p̄[3,2,ind])
                    L̄E[1,ind,nms] += -2*m*h[2,2,nms]*(p̄[2,3,ind] + p̄[3,2,ind])
                    L̄E[2,ind,nms] += -2*m*h[1,2,nms]*(p̄[2,3,ind] + p̄[3,2,ind])
                    h̄[1,2,nms] += -2*m*LE[2,ind,nms]*(p̄[2,3,ind] + p̄[3,2,ind])
                    h̄[2,2,nms] += -2*m*LE[1,ind,nms]*(p̄[2,3,ind] + p̄[3,2,ind])

                    # H[3,3,ind]
                    #L̄E[ind,nms] += -2*m^2*(real(h[nms]) - im*imag(h[nms]))*p̄[3,3,ind]
                    #h̄[nms] += -2*m^2*(real(LE[ind,nms]) - im*imag(LE[ind,nms]))*p̄[3,3,ind]
                    L̄E[1,ind,nms] += -2*m^2*h[1,1,nms]*p̄[3,3,ind]
                    L̄E[2,ind,nms] -= -2*m^2*h[2,1,nms]*p̄[3,3,ind]
                    h̄[1,1,nms] += -2*m^2*LE[1,ind,nms]*p̄[3,3,ind]
                    h̄[2,1,nms] -= -2*m^2*LE[2,ind,nms]*p̄[3,3,ind]
                end
            end
        end
        return s̄elf, p̄otential, L̄E, h̄, P̄, r̄

    end
    return update_potential_hessian!(copy(potential_hessian),LE,h,expansion_order,r),potential_pullback

end

ReverseDiff.@grad_from_chainrules update_potential_hessian!(potential_hessian::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           LE::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           h::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           #ht::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           #ht2::AbstractArray{<:ReverseDiff.TrackedReal},
                                                           expansion_order::Val,
                                                           r::ReverseDiff.TrackedReal)
                                                           =#

function ChainRulesCore.rrule(::typeof(M2L_loop!),LE,L,ME,h,expansion_order::Val{P}) where P

    function LE_pullback(L̄E2)
        s̄elf = NoTangent()
        L̄E = L̄E2
        L̄ = zeros(eltype(L),size(L))
        M̄E = zeros(eltype(ME),size(ME))
        h̄ = zeros(eltype(h),size(h))
        P̄ = NoTangent()
        for j in 0:P
            Cnm = odd_or_even(j)
            for k in 0:j
                jks = (j * (j + 1)) >> 1 + k + 1
                
                for n in 0:P
                    for m in -n:-1
                        nms = (n * (n+1)) >> 1 - m + 1
                        jnkm = (j + n)^2 + j + n + m - k + 1
                        # jnkm_max = (P + P)^2 + P + P + -1 - 0 + 1 = (2P)^2 + 2P = 2P(2P+1)
                        for dim in 1:4
                            #M̄E[dim,nms] += conj(L̄E2[dim,jks])*Cnm*h[jnkm]
                            M̄E[1,dim,nms] += Cnm*(L̄E2[1,dim,jks]*h[1,jnkm] + L̄E2[2,dim,jks]*h[2,jnkm])
                            M̄E[2,dim,nms] += Cnm*(L̄E2[1,dim,jks]*h[2,jnkm] - L̄E2[2,dim,jks]*h[1,jnkm])
                            #h̄[jnkm] += L̄E2[dim,jks]*Cnm*ME[dim,nms]
                            h̄[1,jnkm] += Cnm*(L̄E2[1,dim,jks]*ME[1,dim,nms] - L̄E2[2,dim,jks]*ME[2,dim,nms])
                            h̄[2,jnkm] += Cnm*(L̄E2[1,dim,jks]*ME[2,dim,nms] + L̄E2[2,dim,jks]*ME[1,dim,nms])
                        end
                    end
                    for m in 0:n
                        nms = (n * (n+1)) >> 1 + m + 1
                        jnkm = (j + n) * (j + n) + j + n + m - k + 1
                        # jnkm_max = 2P * 2P + 2P + P + P - 0 + 1 = (2P)^2 + 2P + 2P + 1 = 4P^2 + 4P + 1 = (2P + 1)^2
                        Cnm2 = Cnm * odd_or_even((k-m) * (1 >> (k>=m)) + m)
                        for dim in 1:4
                            #M̄E[dim,nms] += L̄E2[dim,jks]*Cnm2*conj(h[jnkm])
                            M̄E[1,dim,nms] += Cnm2*(L̄E2[1,dim,jks]*h[1,jnkm] + L̄E2[2,dim,jks]*h[2,jnkm])
                            M̄E[2,dim,nms] += Cnm2*(-L̄E2[1,dim,jks]*h[2,jnkm] + L̄E2[2,dim,jks]*h[1,jnkm])
                            #h̄[jnkm] += L̄E2[dim,jks]*Cnm2*conj(ME[dim,nms])
                            h̄[1,jnkm] += Cnm2*(L̄E2[1,dim,jks]*ME[1,dim,nms] + L̄E2[2,dim,jks]*ME[2,dim,nms])
                            h̄[2,jnkm] += Cnm2*(-L̄E2[1,dim,jks]*ME[2,dim,nms] + L̄E2[2,dim,jks]*ME[1,dim,nms])
                        end
                    end
                end
            end
        end
        return s̄elf,L̄E,L̄,M̄E,h̄,P̄
    end
    
    return M2L_loop!(copy(LE),L,ME,h,expansion_order), LE_pullback
end

ReverseDiff.@grad_from_chainrules M2L_loop!(LE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            L::AbstractArray{<:ReverseDiff.TrackedReal},
                                            ME::AbstractArray{<:ReverseDiff.TrackedReal},
                                            h::AbstractArray{<:ReverseDiff.TrackedReal},
                                            expansion_order::Val)

#=@grad_from_chainrules_extended (1,) M2L_loop!(LE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            L::AbstractArray{<:ReverseDiff.TrackedReal},
                                            ME::AbstractArray{<:ReverseDiff.TrackedReal},
                                            h::AbstractArray{<:ReverseDiff.TrackedReal},
                                            expansion_order)=#

function ChainRulesCore.rrule(::typeof(M2M_loop!),BM,CM,h,expansion_order::Val{P}) where P

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
                            #C̄M[dim,jlkms] += B̄M2[dim,i_jk]*C*conj(h[lm])
                            C̄M[1,dim,jlkms] += C*(B̄M2[1,dim,i_jk]*h[1,lm] + B̄M2[2,dim,i_jk]*h[2,lm]) # B̄M is complex iirc
                            C̄M[2,dim,jlkms] += C*(-B̄M2[1,dim,i_jk]*h[2,lm] + B̄M2[2,dim,i_jk]*h[1,lm])
                            #h̄[lm] += B̄M2[dim,i_jk]*C*conj(CM[dim,jlkms])
                            h̄[1,lm] += C*(B̄M2[1,dim,i_jk]*CM[1,dim,jlkms] + B̄M2[2,dim,i_jk]*CM[2,dim,jlkms])
                            h̄[2,lm] += C*(-B̄M2[1,dim,i_jk]*CM[2,dim,jlkms] + B̄M2[2,dim,i_jk]*CM[1,dim,jlkms])
                        end
                    end
                    for m in k:min(l,j+k-l)
                        jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                        lm = l * l + l - m + 1
                        oddeven = odd_or_even(k + l + m)
                        for dim in 1:4
                            #C̄M[dim,jlkms] += conj(B̄M2[dim,i_jk])*h[lm]*oddeven
                            C̄M[1,dim,jlkms] += oddeven*(B̄M2[1,dim,i_jk]*h[1,lm] + B̄M2[2,dim,i_jk]*h[2,lm])
                            C̄M[2,dim,jlkms] += oddeven*(B̄M2[1,dim,i_jk]*h[2,lm] - B̄M2[2,dim,i_jk]*h[1,lm])
                            #h̄[lm] += B̄M2[dim,i_jk]*CM[dim,jlkms]*oddeven
                            h̄[1,lm] += oddeven*(B̄M2[1,dim,i_jk]*CM[1,dim,jlkms] - B̄M2[2,dim,i_jk]*CM[2,dim,jlkms])
                            h̄[2,lm] += oddeven*(B̄M2[1,dim,i_jk]*CM[2,dim,jlkms] + B̄M2[2,dim,i_jk]*CM[1,dim,jlkms])
                        end
                    end
                end
            end
        end
        return s̄elf, B̄M, C̄M, h̄, P̄
        
    end
    
    return M2M_loop!(copy(BM),CM,h,expansion_order),BM_pullback

end
ReverseDiff.@grad_from_chainrules M2M_loop!(BM::AbstractArray{<:ReverseDiff.TrackedReal},
                                            CM::AbstractArray{<:ReverseDiff.TrackedReal},
                                            h::AbstractArray{<:ReverseDiff.TrackedReal},
                                            expansion_order::Val)
#=@grad_from_chainrules_extended (1,) M2M_loop!(BM::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            CM::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            h::AbstractArray{<:Complex{<:ReverseDiff.TrackedReal}},
                                            P)=#

function ChainRulesCore.rrule(::typeof(L2L_loop!),CLE,BLE,h,L,expansion_order::Val{P}) where P

    function CLE_pullback(C̄LE2)
        s̄elf = NoTangent()
        C̄LE = zeros(eltype(CLE),size(CLE))
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
                            #B̄LE[dim,nms] += conj(C̄LE2[dim,jks])*h[jnkm]*oddeven
                            B̄LE[1,dim,nms] += oddeven*(C̄LE2[1,dim,jks]*h[1,jnkm] + C̄LE2[2,dim,jks]*h[2,jnkm])
                            B̄LE[2,dim,nms] += oddeven*(C̄LE2[1,dim,jks]*h[2,jnkm] - C̄LE2[2,dim,jks]*h[1,jnkm])
                            #h̄[jnkm] += C̄LE2[dim,jks]*BLE[dim,nms]*oddeven
                            h̄[1,jnkm] += oddeven*(C̄LE2[1,dim,jks]*BLE[1,dim,nms] - C̄LE2[2,dim,jks]*BLE[2,dim,nms])
                            h̄[2,jnkm] += oddeven*(C̄LE2[1,dim,jks]*BLE[2,dim,nms] + C̄LE2[2,dim,jks]*BLE[1,dim,nms])
                        end
                    end
                    for m in 0:n
                        if n-j >= abs(m-k)
                            jnkm = (n - j) * (n - j) + n - j + m - k + 1
                            nms = (n * (n + 1)) >> 1 + m + 1
                            oddeven = odd_or_even((m-k) * (1 >> (m >= k)))
                            for dim in 1:4
                                #L[dim] += BLE[dim,nms] * h[jnkm] * oddeven
                                #B̄LE[dim,nms] += C̄LE2[dim,jks]*conj(h[jnkm])*oddeven
                                B̄LE[1,dim,nms] += oddeven*(C̄LE2[1,dim,jks]*h[1,jnkm] + C̄LE2[2,dim,jks]*h[2,jnkm])
                                B̄LE[2,dim,nms] += oddeven*(-C̄LE2[1,dim,jks]*h[2,jnkm] + C̄LE2[2,dim,jks]*h[1,jnkm])
                                #h̄[jnkm] += C̄LE2[dim,jks]*conj(BLE[dim,nms])*oddeven
                                h̄[1,jnkm] += oddeven*(C̄LE2[1,dim,jks]*BLE[1,dim,nms] + C̄LE2[2,dim,jks]*BLE[2,dim,nms])
                                h̄[2,jnkm] += oddeven*(-C̄LE2[1,dim,jks]*BLE[2,dim,nms] + C̄LE2[2,dim,jks]*BLE[1,dim,nms])
                            end
                        end
                    end
                end
                #CLE[:,jks] .+= L
            end
        end
        return s̄elf, C̄LE, B̄LE, h̄, L̄, P̄

    end
    return L2L_loop!(copy(CLE),BLE,h,L,expansion_order),CLE_pullback
    #L2L_loop!(copyCLE,BLE,h,L,P)
    #@show sum(CLE)
    #return CLE, CLE_pullback
    
end
ReverseDiff.@grad_from_chainrules L2L_loop!(CLE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            BLE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            h::AbstractArray{<:ReverseDiff.TrackedReal},
                                            L::AbstractArray{<:ReverseDiff.TrackedReal},
                                            expansion_order::Val)

#=@grad_from_chainrules_extended (1,) L2L_loop!(CLE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            BLE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            h::AbstractArray{<:ReverseDiff.TrackedReal},
                                            L::AbstractArray{<:ReverseDiff.TrackedReal},
                                            P)=#

# TODO: Remove depreciated rrules.

function L2B(body_position, expansion_center::SVector{3,TF}, local_expansion, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF<:ReverseDiff.TrackedReal,PS,VPS,VS,GS,P}

    return _L2B(body_position, expansion_center, local_expansion; derivatives_switch=derivatives_switch, expansion_order=expansion_order)

end

function _L2B(body_position, expansion_center::SVector{3,TF}, local_expansion; derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF<:ReverseDiff.TrackedReal,PS,VPS,VS,GS,P}
    error("should not run?")
end

# forwarddiff implementation. Less efficient, but really easy to put together and test.
function ChainRulesCore.rrule(::typeof(_L2B),body_position, expansion_center::SVector{3,TF}, local_expansion; derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF,PS,VPS,VS,GS,P}
    scalar_potential, vector_potential, velocity, gradient = L2B(body_position, expansion_center, local_expansion, derivatives_switch, expansion_order)
    function pullbacks(ȳ)
        println("pullback ran")
        s̄p = view(ȳ[1],:)
        v̄p = view(ȳ[2],:)
        v̄ = view(ȳ[3],:)
        ḡ = view(ȳ[4],:)
        #sp_deriv = zero(TF)
        #vp_derivs = zeros(TF,3)
        #v_derivs = zeros(TF,3)
        #g_derivs = zeros(TF,(3,3))

        body_position_bar = zeros(length(body_position))
        expansion_center_bar = zeros(length(expansion_center))
        local_expansion_bar = zeros(length(local_expansion))
        args = (body_position, expansion_center, local_expansion, derivatives_switch, expansion_order)

        sp_deriv,vp_derivs,v_derivs,g_derivs = forwarddiff_jacobian(L2B,args,1)
        for i=1:length(body_position_bar)
            body_position_bar[i] = s̄p*sp_deriv[i] + sum(v̄p.*vp_derivs[3*(i-1)+1:3*(i-1)+3]) + sum(v̄.*v_derivs[3*(i-1)+1:3*(i-1)+3]) + sum(ḡ.*g_derivs[9*(i-1)+1:9*(i-1)*9])
            expansion_center_bar[i] = -body_position_bar[i] # shortcut the calculation
        end
        sp_deriv,vp_derivs,v_derivs,g_derivs = forwarddiff_jacobian(L2B,args,3)
        for i=1:length(local_expansion_bar)
            local_expansion_bar[i] = s̄p*sp_deriv[i] + sum(v̄p.*vp_derivs[3*(i-1)+1:3*(i-1)+3]) + sum(v̄.*v_derivs[3*(i-1)+1:3*(i-1)+3]) + sum(ḡ.*g_derivs[9*(i-1)+1:9*(i-1)*9])
        end

        return NoTangent(),body_position_bar,expansion_center_bar,local_expansion_bar,NoTangent(),NoTangent()
    end
    return (scalar_potential, vector_potential, velocity, gradient), pullbacks
end

@grad_from_chainrules_multiple_returns L2B(::AbstractArray{<:ReverseDiff.TrackedReal},
                                   ::AbstractArray{<:ReverseDiff.TrackedReal},
                                   ::AbstractArray{<:ReverseDiff.TrackedReal})


function rot(theta,phi)
    st, ct = sincos(theta)
    sp, cp = sincos(phi)
    return SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
end

# take derivative with respect to an input. Not necessary to define, but it makes the code a lot more readable.
@inline forwarddiff_deriv(f,args,idx) = ForwardDiff.derivative(_arg->f(args[1:idx-1]...,_arg,args[idx+1:length(args)]...),args[idx])
# same for the gradient
@inline forwarddiff_grad(f,args,idx) = ForwardDiff.gradient(_arg->f(args[1:idx-1]...,_arg,[args,idx+1:length(args)]...),args[idx])
# and for the jacobian
@inline forwarddiff_jacobian(f,args,idx) = ForwardDiff.jacobian(_arg->f(args[1:idx-1]...,_arg,args[idx+1:length(args)]...),args[idx])

#=
function ChainRulesCore.rrule(::typeof(L2B),body_position, expansion_center::SVector{3,TF}, local_expansion, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF,PS,VPS,VS,GS,P}

    function pullbacks(ȳ)
        #s̄p,v̄p,v̄,ḡ = ȳ # I don't actually want to allocate here... would a view be better?
        s̄p = view(ȳ[1],:)
        v̄p = view(ȳ[2],:)
        v̄ = view(ȳ[3],:)
        ḡ = view(ȳ[4],:)

        # initialize cotangent vectors
        body_position_bar = zeros(length(body_position))
        expansion_center_bar = zeros(length(expansion_center))
        local_expansion_bar = zeros(length(local_expansion))

        # polar coordinate coordinates
        r_bar = zero(TF)
        theta_bar = zero(TF)
        phi_bar = zero(TF)

        # precompute some of the same things as the original function.

        dx, dy, dz = body_position - expansion_center
        r, theta, phi = cartesian_2_spherical(dx, dy, dz)

        scalar_potential = zero(TF)
        vector_potential = zero(SVector{3,TF})
        velocity = zero(SVector{3,TF})
        gradient = zero(SMatrix{3,3,TF,9})
        dudr = zero(SVector{3,TF})
        dudt_r = zero(SVector{3,TF})
        dudp_r_st = zero(SVector{3,TF})
        dudr_bar = zero(SVector{3,TF})
        dudt_r_bar = zero(SVector{3,TF})
        dudp_r_st_bar = zero(SVector{3,TF})

        st, ct = sincos(theta)
        sp, cp = sincos(phi)

        if VS || GS
            R = rot(theta,phi) # need this as a function for later.
            # these next two derivatives are also needed
            dRdtheta = forwarddiff_deriv(rot,(theta,phi),1) # 3x3
            dRdphi = forwarddiff_deriv(rot,(theta,phi),2) # 3x3
        end

        # m = 0, n = 0
        n = m = 0
        rn = 1.0            # r^n
        rnm1 = 0.0          # r^(n-1)
        rnm2 = 0.0          # r^(n-2)
        C_n_m = 1.0           # (-1)^n / (n + |m|)!
        eimp_real = TF(1.0) # real(e^(i m phi))
        eimp_imag = TF(0.0) # imag(e^(i m phi))
        Pnm = 1.0           # associated legendre polynomial of degree n, order m
        alpha_n_m = 0.0     # alpha = Pnm / sin^2(theta)
        beta_n_m = 0.0      # beta = Pnm / sin(theta)
        dbetadt_n_m = 0.0
        dPdt_n_m = 0.0      # derivative w.r.t. theta of P_n_m
        d2Pdt2_n_m = 0.0    # second derivative w.r.t. theta of P_n_m
        index = 1           # index of the local expansion corresponding to n and m

        if PS
            #scalar_potential += L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            # the sign change for the second term is because the adjoint of a complex number is the complex conjugate
            args = (local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            local_expansion_bar[1,1,index] += s̄p*forwarddiff_deriv(L2B_scalar_potential,args,1)
            local_expansion_bar[1,2,index] -= s̄p*forwarddiff_deriv(L2B_scalar_potential,args,2)
            #local_expansion_bar[1,1,index] += ForwardDiff.derivative(LE_real->L2B_scalar_potential(LE_real, local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[1,1,index])
            #local_expansion_bar[1,2,index] -= ForwardDiff.derivative(LE_imag->L2B_scalar_potential(local_expansion[1,1,index], LE_imag, eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[1,2,index])
            # for this case r, θ, ϕ are fixed so derivatives that depend on them are zero.
            # We only need the real part here. The m^2*cos(mϕ) term comes from two applications of the chain rule. TODO: double check this is correct.
            #phi_bar -= m^2*cos(m*phi)*ForwardDiff.derivative(_eimp_imag->L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, _eimp_imag, rn, Pnm, C_n_m),eimp_imag)
            #theta_bar += 0.0 # dPdθ is zero for this case, so there is no theta bar contribution.
            # r_bar += 0.0 # dfdr is zero for this case since r = 1.0 for n=m=0
        end

        if VPS
            #vector_potential += L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            #local_expansion_bar[1,2,index] += ForwardDiff.derivative(LE2_real->L2B_vector_potential(LE2_real, local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[1,2,index])
            #local_expansion_bar[2,2,index] -= ForwardDiff.derivative(LE2_imag->L2B_vector_potential(local_expansion[1,2,index], LE2_imag, local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[2,2,index])
            #local_expansion_bar[1,3,index] += ForwardDiff.derivative(LE3_real->L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], LE3_real, local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[1,3,index])
            #local_expansion_bar[2,3,index] -= ForwardDiff.derivative(LE3_imag->L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], LE3_imag, local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[2,3,index])
            #local_expansion_bar[1,4,index] += ForwardDiff.derivative(LE4_real->L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], LE4_real, local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[1,4,index])
            #local_expansion_bar[2,4,index] -= ForwardDiff.derivative(LE4_imag->L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], LE4_imag, eimp_real, eimp_imag, rn, Pnm, C_n_m),local_expansion[2,4,index])
            args = (local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            for i=2:4
                # Using sum is easier to write than looping over indices of v̄p. A lower-allocation implementation would allocate space to evaluate the gradients in-place. 
                # Also, a lower-allocation version would evaluate analytic derivatives but that's a lot of work.
                local_expansion_bar[1,i,index] += sum(v̄p .* forwarddiff_grad(L2B_vector_potential,args,2*i-1))
                local_expansion_bar[2,i,index] -= sum(v̄p .* forwarddiff_grad(L2B_vector_potential,args,2*i))
            end

            #phi_bar -= m^2*cos(m*phi)*ForwardDiff.derivative(_eimp_imag->L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, _eimp_imag, rn, Pnm, C_n_m), eimp_imag)
        end

        #--- m = 0, n = 1 ---#

        if P > 0
            n = 1
            m = 0
            
            #####
            ##### recurse
            #####

            # r
            rnm2 = rnm1
            rnm1 = rn
            rn *= r
            
            # Pnm
            Pnm = ct
            dPdt_n_m = -st
            d2Pdt2_n_m = -ct

            # C_n_m
            C_n_m *= -1

            # eimp is unchanged
            # alpha = 0
            beta_n_m = zero(TF)
            
            # index
            index += 1

            #####
            ##### evaluate
            #####

            if PS
                # scalar_potential += L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                args = (local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                local_expansion_bar[1,1,index] += s̄p*forwarddiff_deriv(L2B_scalar_potential,args,1)
                local_expansion_bar[2,1,index] -= s̄p*forwarddiff_deriv(L2B_scalar_potential,args,2)
                # e^imϕ is still 1, so no derivative for ϕ
                # however, r and θ are nonzero now.

                r_bar += s̄p*forwarddiff_deriv(L2B_scalar_potential,args,5)#*n*r^(n-1) n=1 -> n*r^(n-1) = 1
                theta_bar += s̄p*dPdt_n_m*forwarddiff_deriv(L2B_scalar_potential,args,6)

            end

            if VPS
                #vector_potential += L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                args = (local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                for i=2:4
                    local_expansion_bar[1,i,index] += sum(v̄p .* forwarddiff_deriv(L2B_vector_potential,args,2*i-1))
                    local_expansion_bar[2,i,index] -= sum(v̄p .* forwarddiff_deriv(L2B_vector_potential,args,2*i))
                end
                #phi, r, theta:
                #phi_bar += 0.0
                r_bar += sum(v̄p .* forwarddiff_deriv(L2B_vector_potential,args,9))#*n*r^(n-1) n=1 -> n*r^(n-1) = 1)
                theta_bar += dPdt_n_m*sum(v̄p .* forwarddiff_deriv(L2B_vector_potential,args,10))
            end

            if VS
                #velocity += L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                args = (local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                for i=1:4
                    local_expansion_bar[1,i,index] += sum(v̄ .* forwarddiff_deriv(L2B_velocity,args,2*i-1))
                    local_expansion_bar[1,i,index] -= sum(v̄ .* forwarddiff_deriv(L2B_velocity,args,2*i))
                end
                #phi_bar += 0.0 # e^imϕ is still 1.0 for this case since m=0
                #r_bar += sum(v̄ .* forwarddiff_deriv(L2B_velocity,args,11))*(n-1)*rnm2 # n-1 = 0 for this case
                theta_bar += dPdt_n_m*sum(v̄ .* forwarddiff_deriv(L2B_velocity,args,12))

                # also need contributions from R, dPdt_n_m, and beta_n_m to theta/phi:
                # R contribution:
                # theta_bar += ((d L2B / dR)(dR/dtheta)) v̄
                # phi_bar += ((d L2B / dR)(dR/dphi)) v̄
                dL2BdR = forwarddiff_jacobian(L2B_velocity,args,16) # 3x3x3, squashed to 3x27.
                # dRdtheta and dRdphi only need to be calculated once per function call
                #dRdtheta = forwarddiff_deriv(rot,(theta,phi),1) # 3x3
                #dRdphi = forwarddiff_deriv(rot,(theta,phi),2) # 3x3
                for i=1:3
                    # j is actually two indices squished into one linear index.
                    for j=1:9
                        theta_bar += v̄[i]*dL2BdR[i,j]*dRdtheta[j]
                        phi_bar += v̄[i]*dL2BdR[i,j]*dRdphi[j]
                    end
                end
                # theta_bar += (d L2B / d(dPdt_n_m)) d2Pdt2_n_m v̄
                theta_bar += sum(forwarddiff_deriv(L2B_velocity,args,13) .* v̄) * d2Pdt2_n_m
                # theta_bar += (d L2B / d(beta_n_m)) d(beta_n_m)/dtheta v̄
                # theta_bar += 0.0
                # beta is zero for this case.
                
            end
            
            #--- m = 1, n = 1 ---#

            #####
            ##### recurse
            #####

            n = m = 1

            # beta
            beta_n_m = -1.0
            beta_nm1_m = 0.0
            beta_nm2_m = 0.0

            # we also need dbeta/dtheta:
            dbetadt_nm2_m = zero(TF)
            dbetadt_nm1_m = zero(TF)
            dbetadt_n_m = 3*st

            # Pnm
            Pnm = -st
            dPdt_n_m = -ct
            d2Pdt2_n_m = st

            # C_n_m
            n_plus_m = n + m
            C_n_m /= n_plus_m

            # eimp
            eimp_real, eimp_imag = complex_multiply(eimp_real, eimp_imag, cp, sp)

            # index
            index += 1

            # r is unchanged

            #####
            ##### evaluate
            #####

            # same things as before, but now phi derivatives matter. beta is still zero though.

            if PS
                #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                args  = (local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                local_expansion_bar[1,1,index] += 2*s̄p*forwarddiff_deriv(L2B_scalar_potential,args,1)
                local_expansion_bar[2,1,index] -= 2*s̄p*forwarddiff_deriv(L2B_scalar_potential,args,2)
                phi_bar += 2*s̄p*(-forwarddiff_deriv(L2B_scalar_potential,args,3)*eimp_imag + forwarddiff_deriv(L2B_scalar_potential,args,4)*eimp_real)*m
                r_bar += 2*s̄p*forwarddiff_deriv(L2B_scalar_potential,args,5) # n = 1, so d(rn)/dr = 1
                theta_bar += 2*s̄p*forwarddiff_deriv(L2B_scalar_potential,args,6) * dPdt_n_m
                
            end

            if VPS
                #vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                args = (local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                for i=1:3
                    local_expansion_bar[1,i,index] += 2*sum(v̄p.*forwarddiff_deriv(L2B_vector_potential,args,2*i-1))
                    local_expansion_bar[2,i,index] -= 2*sum(v̄p.*forwarddiff_deriv(L2B_vector_potential,args,2*i))
                end
                phi_bar += 2*sum(-v̄p.*(forwarddiff_deriv(L2B_vector_potential,args,7)*eimp_imag .+ forwarddiff_deriv(L2B_vector_potential,args,8)*eimp_real))*m
                r_bar += 2*sum(v̄p.*forwarddiff_deriv(L2B_vector_potential,args,9)) # n=1 -> r^n = r, so no extra chain rule terms for this case
                theta_bar += 2*sum(v̄p.*forwarddiff_deriv(L2B_vector_potential,args,10))*dPdt_n_m
            end

            if VS
                #velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                args = (local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                for i=1:4
                    local_expansion_bar[1,i,index] += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,2*i-1))
                    local_expansion_bar[2,i,index] -= 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,2*i))
                end
                phi_bar += 2*sum(-v̄.*(forwarddiff_deriv(L2B_vector_potential,args,7)*eimp_imag .+ forwarddiff_deriv(L2B_vector_potential,args,8)*eimp_real))*m
                #r_bar += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,11)) * 0 # rnm1 = 1 since m = 1, so there's no contribution from this term
                theta_bar += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,12))*dPdt_n_m
                theta_bar += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,13))*d2Pdt2_n_m
                # no beta contribution since beta is still constant
                # R contribution:
                dL2BdR = forwarddiff_jacobian(L2B_velocity,args,16) # 3x3x3, squashed to 3x27.
                # dRdtheta and dRdphi only need to be calculated once
                #dRdtheta = forwarddiff_deriv(rot,(theta,phi),1) # 3x3
                #dRdphi = forwarddiff_deriv(rot,(theta,phi),2) # 3x3
                for i=1:3
                    for j=1:9
                        theta_bar += 2*v̄[i]*dL2BdR[i,j]*dRdtheta[j]
                        phi_bar += 2*v̄[i]*dL2BdR[i,j]*dRdphi[j]
                    end
                end
            end
        end

        #--- m = 0, m = 1, n > 1 ---#

        if P > 1
            
            m = 1
            alpha_n_0 = zero(TF)        # n>1, m<2
            beta_n_0 = zero(TF)         # n>1, m<1
            P_nm2_0 = 0.0               # n=-1,m=0
            P_nm1_0 = 1.0               # n=0, m=0
            P_n_0 = ct                  # n=1, m=0
            dPdt_nm1_0 = zero(TF)       # n=0, m=0
            dPdt_n_0 = -st              # n=1, m=0
            d2Pdt2_nm1_0 = 0.0          # n=0, m=0
            d2Pdt2_n_0 = dPdt_n_m       # n=1, m=0
            d3Pdt3_nm2_0 = 0.0          # n=-1,m=0
            d3Pdt3_nm1_0 = 0.0          # n=0, m=0
            d3Pdt3_n_0 = st             # n=1, m=0
            dPdct_n_0 = 1.0             # n=1, m=0
            d2Pdct2_n_0 = 0.0           # n=1, m=0
            dPdct_nm1_0 = 0.0           # n=0, m=0
            d2Pdct2_nm1_0 = 0.0         # n=0, m=0
            P_nm1_m = 0.0               # n=0, m=1

            for n in 2:P

                #####
                ##### recurse
                #####
                
                # rn
                rnm2 = rnm1
                rnm1 = rn
                rn *= r
                
                # C_n_m
                n_plus_m += 1
                C_n_0 = -C_n_m
                C_n_m /= -n_plus_m
                
                # beta
                beta_nm2_m = beta_nm1_m
                beta_nm1_m = beta_n_m
                beta_n_m = ((2*n-1) * ct * beta_nm1_m - n * beta_nm2_m) / (n-1)

                # we also need dbeta/dtheta:
                dbetadt_nm2_m = dbetadt_nm1_m
                dbetadt_nm1_m = dbetadt_n_m
                dbetadt_n_m = ((2*n-1)*(-st*beta_nm1_m + ct*dbetadt_nm1_m) - n*dbetaedt_nm2_m)/(n-1)
                
                # Pnm, m=1
                P_nm1_m = Pnm
                Pnm = beta_n_m * st

                # Pnm, m=0
                P_nm2_0 = P_nm1_0
                P_nm1_0 = P_n_0
                P_n_0 = 1/n*((2*n-1) * ct * P_nm1_0 - (n-1) * P_nm2_0)
                
                # first derivatives, m=1
                dPdt_n_m = n * ct * beta_n_m - (n + m) * beta_nm1_m
                
                # first derivatives, m=0
                dPdt_nm1_0 = dPdt_n_0
                dPdt_n_0 = ct * dPdt_nm1_0 - n * st * P_nm1_0
                
                # second derivatives
                d2Pdt2_nm1_0 = d2Pdt2_n_0
                d2Pdt2_n_0 = dPdt_n_m
                d3Pdt3_nm2_0 = d3Pdt3_nm1_0
                d3Pdt3_nm1_0 = d3Pdt3_n_0
                d3Pdt3_n_0 = d3Pdt3_nm2_0 - (2*n-1) * (st*d2Pdt2_nm1_0 + 2*ct*P_nm1_m - st*P_nm1_0)
                d2Pdt2_n_m = d3Pdt3_n_0

                # second derivatives over sin(theta)
                dPdct_nm1_0 = dPdct_n_0
                d2Pdct2_nm1_0 = d2Pdct2_n_0
                dPdct_n_0 = n * P_nm1_0 + ct * dPdct_nm1_0
                d2Pdct2_n_0 = (n+1) * dPdct_nm1_0 + ct * d2Pdct2_nm1_0
                ddt_Pnm_st = st * d2Pdct2_n_0

                # index
                index = (n * (n + 1)) >> 1 + m + 1
                index_m0 = (n * (n + 1)) >> 1 + 0 + 1

                # eimp is unchanged

                #####
                ##### evaluate
                #####

                if PS
                    #scalar_potential += L2B_scalar_potential(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)
                    args = (local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)

                    local_expansion_bar[1,1,index_m0] += s̄*forwarddiff_deriv(L2B_scalar_potential,args,1)
                    local_expansion_bar[2,1,index_m0] -= s̄*forwarddiff_deriv(L2B_scalar_potential,args,2)
                    # next two entries are explicity constraints
                    # n is now > 1
                    r_bar += s̄*forwarddiff_deriv(scalar_potential,args,5) * n * rnm1
                    # check to make sure initializations are consistent
                    if abs(rnm1/r - rnm2) > 1e-12
                        error("rnm1/r is not equal to rnm2: $(rnm1/r) ≠ $rnm2")
                    end
                    theta_bar += s̄*forwarddiff_deriv(scalar_potential,args,6) * dPdt_n_0

                    #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                    args = (local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)

                    local_expansion_bar[1,1,index_m0] += 2*s̄*forwarddiff_deriv(L2B_scalar_potential,args,1)
                    local_expansion_bar[2,1,index_m0] -= 2*s̄*forwarddiff_deriv(L2B_scalar_potential,args,2)

                    phi_bar += 2*s̄*(-forwarddiff_deriv(L2B_scalar_potential,args,3)*eimp_imag + forwarddiff_deriv(L2B_scalar_potential,args,4)*eimp_real)*m
                    r_bar += 2*s̄*forwarddiff_deriv(scalar_potential,args,5) * n * rnm1
                    theta_bar += 2*s̄*forwarddiff_deriv(scalar_potential,args,6) * dPdt_n_m

                end
        
                if VPS
                    #vector_potential += L2B_vector_potential(local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)
                    args = (local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)
                    
                    for i=2:4
                        local_expansion_bar[1,i,index_m0] += sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,2*i-1))
                        local_expansion_bar[2,i,index_m0] -= sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,2*i))
                    end
                    # no phi contribution for this case
                    r_bar += sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,9))*n*rnm1
                    theta_bar += sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,10))*dPdt_n_0

                    #vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                    args = (local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)

                    for i=2:4
                        local_expansion_bar[1,i,index_m0] += 2*sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,2*i-1))
                        local_expansion_bar[2,i,index_m0] -= 2*sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,2*i))
                    end
                    phi_bar += 2*sum(v̄p.*(-forwarddiff_deriv(L2B_scalar_potential,args,3)*eimp_imag .+ forwarddiff_deriv(L2B_scalar_potential,args,4)*eimp_real))*m
                    r_bar += 2*sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,9))*n*rnm1
                    theta_bar += 2*sum(v̄p.*forwarddiff_deriv(L2B_velocity,args,10))*dPdt_n_m

                end
        
                if VS
                    #velocity += L2B_velocity(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm1, P_n_0, dPdt_n_0, beta_n_0, C_n_0, R, n, 0, derivatives_switch)
                    args = (local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm1, P_n_0, dPdt_n_0, beta_n_0, C_n_0, R, n, 0, derivatives_switch)

                    for i=1:4
                        local_expansion_bar[1,i,index_m0] += sum(v̄.*forwarddiff_deriv(L2B_velocity,args,2*i-1))
                        local_expansion_bar[2,i,index_m0] -= sum(v̄.*forwarddiff_deriv(L2B_velocity,args,2*i))
                    end
                    # no phi contribution in the normal spot... but it will contribute to R
                    r_bar += sum(v̄.*forwarddiff_deriv(L2B_velocity,args,11))*(n-1)*rnm2
                    if abs(rnm1/r - rnm2) > 1e-12
                        error("rnm1/r is not equal to rnm2: $(rnm1/r) ≠ $rnm2")
                    end
                    theta_bar += sum(v̄.*forwarddiff_deriv(L2B_velocity,args,12))*dPdt_n_0
                    theta_bar += sum(v̄.*forwarddiff_deriv(L2B_velocity,args,13))*d2Pdt2_n_0
                    #theta_bar += sum(v̄.*forwarddiff_deriv(L2B_velocity,args,14))*dbetadt_n_0 # need to calculuate dbetadt_n_m. However, beta_n_0 is zero so this term can be ignored the first time.
                    dL2BdR = forwarddiff_jacobian(L2B_velocity,args,16)
                    for i=1:3
                        for j=1:9
                            theta_bar += v̄[i]*dL2BdR[i,j]*dRdtheta[j]
                            phi_bar += v̄[i]*dL2BdR[i,j]*dRdphi[j]
                        end
                    end
                    
                    #velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                    args = (local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                    for i=1:4
                        local_expansion_bar[1,i,index_m0] += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,2*i-1))
                        local_expansion_bar[2,i,index_m0] -= 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,2*i))
                    end
                    phi_bar += 2*sum(v̄.*(-forwarddiff_deriv(L2B_scalar_potential,args,3)*eimp_imag .+ forwarddiff_deriv(L2B_scalar_potential,args,4)*eimp_real))*m
                    r_bar += 2*sum(v̄.(forwarddiff_deriv(L2B_velocity,args,11)))*(n-1)*rnm2
                    theta_bar += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,12))*dPdt_n_0
                    theta_bar += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,13))*d2Pdt2_n_0
                    theta_bar += 2*sum(v̄.*forwarddiff_deriv(L2B_velocity,args,14))*dbetadt_n_m

                    dL2BdR = forwarddiff_jacobian(L2B_velocity,args,16)
                    for i=1:3
                        for j=1:9
                            theta_bar += 2*v̄[i]*dL2BdR[i,j]*dRdtheta[j]
                            phi_bar += 2*v̄[i]*dL2BdR[i,j]*dRdphi[j]
                        end
                    end

                end

                if GS
                    # m = 0
                    ddt_Pn0_st = zero(TF)
                    #this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm2, P_n_0, dPdt_n_0, d2Pdt2_n_0, ddt_Pn0_st, alpha_n_0, beta_n_0, st, ct, sp, cp, C_n_0, R, n, 0, derivatives_switch)
                    args = (local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm2, P_n_0, dPdt_n_0, d2Pdt2_n_0, ddt_Pn0_st, alpha_n_0, beta_n_0, st, ct, sp, cp, C_n_0, R, n, 0, derivatives_switch)
                    _f(_args) = hcat(L2B_velocity_gradient(_args)...)
                    # for the local_expansion cotangents, the extra contribution from the rotation matrix disappears since dR/d(local expansion) = 0
                    for i=1:4
                        local_expansion_bar[1,i,index_m0] += sum(ḡ'*R*forwarddiff_deriv(_f,args,2*i-1))
                        local_expansion_bar[1,i,index_m0] -= sum(ḡ'*R*forwarddiff_deriv(_f,args,2*i))
                    end
                    # first ϕ contributions are zero due to fixed ϕ
                    dghatdr = forwarddiff_deriv(_f,args,11)
                    dghatdPn0 = forwarddiff_deriv(_f,args,12)
                    dghatddPdtn0 = forwarddiff_deriv(_f,args,13)
                    dghatd2Pdt2n0 = forwarddiff_deriv(_f,args,14)
                    #dghatdddtPn0 = forwarddiff_deriv(_f,args,15) # ddt_Pn0_st is zero
                    #dghatdalphan0 = forwarddiff_deriv(_f,args,16) # alpha_n_0 is zero
                    #dghatdbetan0 = forwarddiff_deriv(_f,args,17) # beta_n_0 is zer0
                    dghatdst = forwarddiff_deriv(_f,args,18)
                    dghatdct = forwarddiff_deriv(_f,args,19)
                    dghatdsp = forwarddiff_deriv(_f,args,20)
                    dghatdcp = forwarddiff_deriv(_f,args,21)
                    dghatdR = forwarddiff_jacobian(_f,args,23) # 81 entries in a 3x3x3x3 tensor, but indexed in a 9x9 array.
                    for i=1:3
                        for j=1:3
                            for k=1:3
                                r_bar += ḡ[k,i]*R[j,i]*dghatdr[j,k]*(n-2)*rnm2/r
                                theta_bar += ḡ[k,i]*R[j,i]*(dghatdPn0[j,k]*dPdt_n_0 + dghatddPdtn0[j,k]*d2Pdt2_n_0 + dghatd2Pdt2n0[j,k]*d3Pdt3_n_0 + dghatdst[j,k]*ct - dghatdct[j,k]*st + dghatdsp[j,k]*cp - dghatdcp[j,k]*sp + sum(dghatdR[27*(j-1)+9*(k-1)+1:27*(j-1)+9*(k-1)+9].*dRdtheta[1:end]))
                                phi_bar += ḡ[k,i]*R[j,i]*(sum(dghatdR[27*(j-1)+9*(k-1)+1:27*(j-1)+9*(k-1)+9].*dRdphi[1:end]))
                            end
                        end
                    end
                    
                    #dudr += this_dudr
                    #dudt_r += this_dudt_r
                    #dudp_r_st += this_dudp_r_st
                    
                    # m = 1

                    #this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, Pnm, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_0, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)
                    
                    #dudr += 2 * this_dudr
                    #dudt_r += 2 * this_dudt_r
                    #dudp_r_st += 2 * this_dudp_r_st

                    args = (local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, Pnm, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_0, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)
                    _f(_args) = hcat(L2B_velocity_gradient(_args)...)
                    # for the local_expansion cotangents, the extra contribution from the rotation matrix disappears since dR/d(local expansion) = 0
                    for i=1:4
                        local_expansion_bar[1,i,index_m0] += 2*sum(ḡ'*R*forwarddiff_deriv(_f,args,2*i-1))
                        local_expansion_bar[1,i,index_m0] -= 2*sum(ḡ'*R*forwarddiff_deriv(_f,args,2*i))
                    end
                    dghatdeimp_real = forwarddiff_deriv(_f,args,9)
                    dghatdeimp_imag = forwarddiff_deriv(_f,args,10)
                    dghatdr = forwarddiff_deriv(_f,args,11)
                    dghatdPnm = forwarddiff_deriv(_f,args,12)
                    dghatddPdtnm = forwarddiff_deriv(_f,args,13)
                    dghatd2Pdt2nm = forwarddiff_deriv(_f,args,14)
                    # alpha_n_1 = zero(TF) # this goes to zero, so the alpha contribution ends up not mattering. This makes the third-derivative term of Pnm also go to zero... I think. I'll have to double check that this is correct.
                    #d3Pdt3_n_m = 2*n^2*st*ct*alpha_n_m - (n*(n*st^2 + 1) - 1) * dalphadt_n_1
                    dghatdddtPnmst = forwarddiff_deriv(_f,args,15)
                    dddt_Pnm_stdt = (dPdt_n_m - ct/st*Pnm)/st
                    # dghatdalphan0 = forwarddiff_deriv(_f,args,16) # alpha_n_0 is zero
                    dghatdbetanm = forwarddiff_deriv(_f,args,17)
                    dghatdst = forwarddiff_deriv(_f,args,18)
                    dghatdct = forwarddiff_deriv(_f,args,19)
                    dghatdsp = forwarddiff_deriv(_f,args,20)
                    dghatdcp = forwarddiff_deriv(_f,args,21)
                    dghatdR = forwarddiff_jacobian(_f,args,23) # 81 entries in a 3x3x3x3 tensor, but indexed in a 9x27 array.
                    for i=1:3
                        for j=1:3
                            for k=1:3
                                r_bar += 2*ḡ[k,i]*R[j,i]*dghatdr[j,k]*(n-2)*rnm2/r
                                theta_bar += 2*ḡ[k,i]*R[j,i]*(dghatdPnm[j,k]*dPdt_n_m + dghatddPdtnm[j,k]*d2Pdt2_n_m + dghatdddtPnmst[j,k]*dddt_Pnm_stdt + dghatdbetanm[j,k]*dbetadt_n_m + dghatdst[j,k]*ct - dghatdct[j,k]*st + dghatdsp[j,k]*cp - dghatdcp[j,k]*sp + sum(dghatdR[27*(j-1)+9*(k-1)+1:27*(j-1)+9*(k-1)+9].*dRdtheta[1:end]))
                                phi_bar += 2*ḡ[k,i]*R[j,i]*(-m*forwarddiff_deriv(_f,args,9)*dghatdeimp_imag + m*forwarddiff_deriv(_f,args,10)*dghatdeimp_real + sum(dghatdR[27*(j-1)+9*(k-1)+1:27*(j-1)+9*(k-1)+9].*dRdphi[1:end]))
                            end
                        end
                    end
                end

            end

            #--- m > 1, n > 1 ---#

            # double factorial
            # double_factorial = 3

            # sin^(n-2)(theta)
            # st_nm2 = one(TF)

            # alpha_n_n
            alpha_n_n = TF(3.0) # sin(theta)^(n-2) (2*n-1)!! for n=2

            # rn
            rnm2 = zero(TF) # n=-1...
            rnm1 = one(TF)  # n=0
            rn = r          # n=1

            # C_n_m: n = m = 1
            n_plus_m = 2
            C_n_m = -0.5 # (-1)^1 / (1 + |1|)!

            for m in 2:P

                #####
                ##### recurse
                #####

                n = m

                # rn
                rnm2 = rnm1
                rnm1 = rn
                rn *= r
                
                # C_n_m: increment n
                n_plus_m += 1
                C_n_m /= -n_plus_m

                # increment m
                n_plus_m += 1
                C_n_m /= n_plus_m

                # alpha_n_n = (-1)^n * st^(n-2) * double_factorial(2*n-1)
                alpha_n_m = alpha_n_n
                alpha_nm1_m = 0.0
                alpha_nm2_m = 0.0

                # beta
                beta_n_m = alpha_n_m * st
                beta_nm1_m = 0.0

                # P_n_m
                P_n_m = beta_n_m * st

                # first derivative
                dPdt_n_m = n * ct * beta_n_m - (n + m) * beta_nm1_m

                # second derivative
                d2Pdt2_n_m = (n + m) * ct * alpha_nm1_m - (n * (n * st^2 + 1) - m^2) * alpha_n_m

                # derivative of Pnm/sin(theta)
                ddt_Pnm_st = (n-1) * ct * alpha_n_m - (n+m) * alpha_nm1_m

                # eimp
                eimp_real, eimp_imag = complex_multiply(eimp_real, eimp_imag, cp, sp)

                # index
                index = (n * (n + 1)) >> 1 + m + 1
                
                #####
                ##### evaluate
                #####

                if PS
                    #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
                    args = (local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
                    local_expansion_bar[1,1,index] += s̄*forwarddiff_deriv(L2B_scalar_potential,args,1)
                    local_expansion_bar[2,1,index] -= s̄*forwarddiff_deriv(L2B_scalar_potential,args,2)
                    phi_bar += s̄*(-forwarddiff_deriv(L2B_scalar_potential,args,3)*eimp_imag + forwarddiff_deriv(L2B_scalar_potential,args,4)*eimp_real)*m
                    r_bar += s̄*forwarddiff_deriv(L2B_scalar_potential,args,5)*m*rnm1 # rn is actually rm here -_-
                    theta_bar += s̄*forwarddiff_deriv(L2B_scalar_potential,args,6)*dPdt_n_m
                    
                end
        
                if VPS
                    vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
                    
                end
        
                if VS
                    velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, P_n_m, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                end

                if GS
                    this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)
                    
                    dudr += 2 * this_dudr
                    dudt_r += 2 * this_dudt_r
                    dudp_r_st += 2 * this_dudp_r_st
                end
                
                # prepare to recurse C_n_m
                C_n_m_loop = C_n_m
                n_plus_m_loop = n_plus_m
                rn_loop = rn
                rnm1_loop = rnm1
                rnm2_loop = rnm2

                for n in m+1:P

                    #####
                    ##### recurse
                    #####

                    # rn
                    rnm2_loop = rnm1_loop
                    rnm1_loop = rn_loop
                    rn_loop *= r

                    # alpha
                    alpha_nm2_m = alpha_nm1_m
                    alpha_nm1_m = alpha_n_m
                    alpha_n_m = 1/(n-m) * ( (2*n-1)*ct*alpha_nm1_m - (n+m-1)*alpha_nm2_m )

                    # beta
                    beta_nm1_m = beta_n_m
                    beta_n_m = alpha_n_m * st

                    # P_n_m
                    P_n_m = beta_n_m * st

                    # first derivative
                    dPdt_n_m = n * ct * beta_n_m - (n + m) * beta_nm1_m

                    # second derivative
                    d2Pdt2_n_m = (n + m) * ct * alpha_nm1_m - (n * (n * st^2 + 1) - m^2) * alpha_n_m

                    # derivative of Pnm/sin(theta)
                    ddt_Pnm_st = (n-1) * ct * alpha_n_m - (n+m) * alpha_nm1_m

                    # C_n_m
                    n_plus_m_loop += 1
                    C_n_m_loop /= -n_plus_m_loop

                    # index
                    index = (n * (n + 1)) >> 1 + m + 1

                    #####
                    ##### evaluate
                    #####
                    if PS
                        scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop)
                    end
            
                    if VPS
                        vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop)
                    end
            
                    if VS
                        velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1_loop, P_n_m, dPdt_n_m, beta_n_m, C_n_m_loop, R, n, m, derivatives_switch)
                    end
        
                    if GS
                        this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2_loop, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m_loop, R, n, m, derivatives_switch)
                        
                        dudr += 2 * this_dudr
                        dudt_r += 2 * this_dudt_r
                        dudp_r_st += 2 * this_dudp_r_st
                    end
                end

                # tail recursion
                alpha_n_n *= -st * (2*n+1)
                # double_factorial *= 2*n+1
                # st_nm2 *= -st

            end

        end

        # rotate to cartesian
        # R = SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
        duidrk = hcat(dudr, dudt_r, dudp_r_st)
        gradient = duidrk * R'
    

        # return input cotangents

        return NoTangent(),body_position_bar,expansion_center_bar,local_expansion_bar,NoTangent(),NoTangent()
    end

    return (scalar_potential, vector_potential, velocity, gradient), pullbacks
end
=#