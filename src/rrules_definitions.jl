# Defines rrules for ChainRules and registers them with ReverseDiff.
tracked_type = Union{ReverseDiff.TrackedReal,ReverseDiff.TrackedArray} # tracked type used for registering ChainRulesCore rrules with ReverseDiff.

ϵ(i::TI,j::TI,k::TI) where TI <: Int = (i == j || j == k || k == i) ? 0 : ((i + j) % 3 == 1 ? 1 : -1)

"""
function regular_harmonic!(harmonics, rho, theta, phi, P)
"""

# maybe a typo somewhere here?
function ChainRulesCore.rrule(::typeof(regular_harmonic!), harmonics, rho, theta, phi, P)
    #println("running regular harmonic in pullback:")
    #@time regular_harmonic!(harmonics,rho,theta,phi,P)
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
#ReverseDiff.@grad_from_chainrules regular_harmonic!(harmonics::AbstractArray{<:ReverseDiff.TrackedReal}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)

# nans show up here
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
#ReverseDiff.@grad_from_chainrules irregular_harmonic!(harmonics::AbstractArray{<:ReverseDiff.TrackedReal}, rho::tracked_type, theta::tracked_type, phi::tracked_type, P)

#=
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
=#

#=function ChainRulesCore.rrule(::typeof(M2L_loop!),LE,L,ME,h,expansion_order::Val{P}) where P

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
                                            expansion_order::Val)=#

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
    
end

#=ReverseDiff.@grad_from_chainrules L2L_loop!(CLE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            BLE::AbstractArray{<:ReverseDiff.TrackedReal},
                                            h::AbstractArray{<:ReverseDiff.TrackedReal},
                                            L::AbstractArray{<:ReverseDiff.TrackedReal},
                                            expansion_order::Val)=#

#=function rot(theta,phi)
    st, ct = sincos(theta)
    sp, cp = sincos(phi)
    return SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
end=#

# take derivative with respect to an input. Not necessary to define, but it makes the code a lot more readable.
@inline forwarddiff_deriv(f,args,idx) = ForwardDiff.derivative(_arg->f(args[1:idx-1]...,_arg,args[idx+1:length(args)]...),args[idx])
# same for the gradient
@inline forwarddiff_grad(f,args,idx) = ForwardDiff.gradient(_arg->f(args[1:idx-1]...,_arg,[args,idx+1:length(args)]...),args[idx])
# and for the jacobian
@inline forwarddiff_jacobian(f,args,idx) = ForwardDiff.jacobian(_arg->f(args[1:idx-1]...,_arg,args[idx+1:length(args)]...),args[idx])

function _L2B(args...)
    error("Dummy function _L2B was called. This function is only defined for the purposes of developing an L2B pullback and will be removed once that is done.")
end

# passes tests
function scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, s̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=false)
    
    local_expansion_bar[1,1,index] += eimp_real*rn*Pnm*C_n_m*s̄p*(scale_by_two ? 2 : 1)
    local_expansion_bar[2,1,index] -= eimp_imag*rn*Pnm*C_n_m*s̄p*(scale_by_two ? 2 : 1)
    r_theta_phi_bar[1] += n/r*s̄p*L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)*(scale_by_two ? 2 : 1)
    r_theta_phi_bar[2] += dPdt_n_m/Pnm*s̄p*L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)*(scale_by_two ? 2 : 1)
    r_theta_phi_bar[3] += -m*s̄p*(local_expansion[2,1,index]*eimp_real + local_expansion[1,1,index]*eimp_imag)*rn*Pnm*C_n_m*(scale_by_two ? 2 : 1)
    return nothing

end

# partially tested
function vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, v̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=false)
   
    for i=1:3
        local_expansion_bar[1,i+1,index] += eimp_real*rn*Pnm*C_n_m*v̄p[i]*(scale_by_two ? 2 : 1)
        local_expansion_bar[2,i+1,index] -= eimp_imag*rn*Pnm*C_n_m*v̄p[i]*(scale_by_two ? 2 : 1)
    end
    r_theta_phi_bar[1] += n / r * sum(v̄p.*L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m))*(scale_by_two ? 2 : 1)
    r_theta_phi_bar[2] += dPdt_n_m / Pnm * sum(v̄p.*L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m))*(scale_by_two ? 2 : 1)
    for i=1:3
        r_theta_phi_bar[3] -= m*rn*Pnm*C_n_m*(local_expansion[2,i+1,index]*eimp_real + local_expansion[1,i+1,index]*eimp_imag)*v̄p[i]*(scale_by_two ? 2 : 1)
    end
    return nothing
end

# something is wrong here...
function velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real::TF, eimp_imag, rnm1, Pnm, dPdt, beta_n_m, C_n_m, R, n, m, derivatives_switch::DerivativesSwitch{PS,VPS,<:Any,<:Any}, index, v̄, r, theta, phi, d2Pdt2, dbetadt_n_m, dRdt, dRdp, v; scale_by_two=false) where {TF, PS, VPS}
    #Extra inputs: d2Pdt2, dbetadt_n_m, dRdt, dRdp
    #TF = promote_type(TF1,TF2)

    # initialize
    #velocity = zero(SVector{3,TF})
    
    # intermediate quantities

    vx_real = n * Pnm
    vy_real = dPdt
    vz_imag = m * beta_n_m
    dvdt_x_real = n * dPdt
    dvdt_y_real = d2Pdt2
    dvdt_z_imag = m*dbetadt_n_m
    #rnm2 = rnm1/r
	
    # rotate to cartesian
    #=
	u_real_cartesian, u_imag_cartesian = complex_multiply(R, ux_real, 0, uy_real, 0, 0, uz_imag)
    dudt_real_cartesian, dudt_imag_cartesian = complex_multiply(R, dudt_x_real, 0, dudt_y_real, 0, 0, dudt_z_imag)
    dRdt_u_real, dRdt_u_imag = complex_multiply(dRdt, ux_real, 0, uy_real, 0, 0, uz_imag)
    dRdp_u_real, dRdp_u_imag = complex_multiply(dRdp, ux_real, 0, uy_real, 0, 0, uz_imag)
    =#
    # these should probably be pre-allocated.
    dRdt_v_real,dRdt_v_imag = complex_multiply(dRdt,vx_real,0,vy_real,0,0,vz_imag)
    R_dvdt_real, R_dvdt_imag = complex_multiply(R,dvdt_x_real,0,dvdt_y_real,0,0,dvdt_z_imag)
    R_v_real,R_v_imag = complex_multiply(R,vx_real,0,vy_real,0,0,vz_imag)
    dRdp_v_real,dRdp_v_imag = complex_multiply(dRdp,vx_real,0,vy_real,0,0,vz_imag)

    # takes care of r̄ contributions with one line.
    for i=1:3
        r_theta_phi_bar[1] += v̄[i]*(n-1)*v[i]/r#*(scale_by_two ? 2 : 1)
    end
    
    # velocity due to scalar potential
    if PS
        Lnm_eimp_real, Lnm_eimp_imag = complex_multiply(eimp_real, eimp_imag, local_expansion[1,1,index], local_expansion[2,1,index])
	
        #velocity -= SVector{3}(
		#	complex_multiply_real(u_real_cartesian[1], u_imag_cartesian[1], Lnm_eimp_real, Lnm_eimp_imag),
		#	complex_multiply_real(u_real_cartesian[2], u_imag_cartesian[2], Lnm_eimp_real, Lnm_eimp_imag),
		#	complex_multiply_real(u_real_cartesian[3], u_imag_cartesian[3], Lnm_eimp_real, Lnm_eimp_imag)
        #) * rnm1 * C_n_m

        for i=1:3
            for j=1:3
                r_theta_phi_bar[3] -= v̄[i]*rnm1*C_n_m*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdt_v_real[j] + R_dvdt_real[j], dRdt_v_imag[j] + R_dvdt_imag[j])*(scale_by_two ? 2 : 1)
                r_theta_phi_bar[2] -= v̄[i]*rnm1*C_n_m*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdp_v_real[j] - m*R_v_imag[j], dRdp_v_imag[j] + m*R_v_real[j])*(scale_by_two ? 2 : 1)
                local_expansion_bar[1,1,index] -= v̄[i]*rnm1*C_n_m*complex_multiply_real(R_v_real[j],R_v_imag[j],eimp_real,eimp_imag)*(scale_by_two ? 2 : 1)
                local_expansion_bar[2,1,index] += v̄[i]*rnm1*C_n_m*complex_multiply_imag(R_v_real[j],R_v_imag[j],eimp_real,eimp_imag)*(scale_by_two ? 2 : 1)
            end
        end
    end

    # velocity due to vector potential
    if VPS

		#Lnm_x_eimp_real, Lnm_x_eimp_imag = complex_multiply(Lnm_x_real, Lnm_x_imag, eimp_real, eimp_imag)
		#Lnm_y_eimp_real, Lnm_y_eimp_imag = complex_multiply(Lnm_y_real, Lnm_y_imag, eimp_real, eimp_imag)
		#Lnm_z_eimp_real, Lnm_z_eimp_imag = complex_multiply(Lnm_z_real, Lnm_z_imag, eimp_real, eimp_imag)

		#velocity += rnm1 * C_n_m * complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_eimp_real, Lnm_x_eimp_imag, Lnm_y_eimp_real, Lnm_y_eimp_imag, Lnm_z_eimp_real, Lnm_z_eimp_imag)

        for k=1:3
            Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(eimp_real,eimp_imag,local_expansion[1,k+1,index],local_expansion[2,k+1,index])
            for j=1:3
                for i=1:3
                    r_theta_phi_bar[3] += ϵ(i,j,k)*v̄[i]*rnm1*C_n_m*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdt_v_real[j] + R_dvdt_real[j], dRdt_v_imag[j] + R_dvdt_imag[j])*(scale_by_two ? 2 : 1)
                    r_theta_phi_bar[2] += ϵ(i,j,k)*v̄[i]*rnm1*C_n_m*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdp_v_real[j] - m*R_v_imag[j], dRdp_v_imag[j] + m*R_v_real[j])*(scale_by_two ? 2 : 1)
                    local_expansion_bar[1,k+1,index] += v̄[i]*rnm1*C_n_m*ϵ(i,j,k)*complex_multiply_real(eimp_real, eimp_imag, R_v_real[j],R_v_imag[j])*(scale_by_two ? 2 : 1)
                    local_expansion_bar[2,k+1,index] -= v̄[i]*rnm1*C_n_m*ϵ(i,j,k)*complex_multiply_imag(eimp_real, eimp_imag, R_v_real[j],R_v_imag[j])*(scale_by_two ? 2 : 1)
                end
            end
        end
    end

    return nothing
end

function gradient_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm2, Pnm, dPdt, d2Pdt2, ddt_Pnm_st, alpha, beta, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch::DerivativesSwitch{PS,VPS,<:Any,<:Any}, index, ḡ, r, dRdt, dRdp, dalphadt, dbetadt, d3Pdt3; scale_by_two=false) where {PS, VPS}
    # extra inputs to pass in: M, dRdt, dRdp, dbetadt_n_m
    # also duidrk = gradient*R

    TF = eltype(local_expansion)

    sb2 = (scale_by_two ? 2 : 1)

    #Lnm_real, Lnm_imag = local_expansion[1,1,index], local_expansion[2,1,index]
    #Lnm_x_real, Lnm_x_imag = local_expansion[1,2,index], local_expansion[2,2,index]
    #Lnm_y_real, Lnm_y_imag = local_expansion[1,3,index], local_expansion[2,3,index]
    #Lnm_z_real, Lnm_z_imag = local_expansion[1,4,index], local_expansion[2,4,index]

    # initialize outputs
    #dudr = zero(SVector{3,TF})
    #dudt_r = zero(SVector{3,TF})
    #dudp_r_st = zero(SVector{3,TF})

    # intermediate values
    rnm2Cnm = rnm2 * C_n_m
    nm1 = n - 1
    rnm3Cnm = rnm2Cnm/r
    nm2_r = (n-2)/r # needed 9 times, so this saves some calculations.
    
    #####
    ##### du/dr
    #####
    
    # intermediate quantities
    b_x_real = n*Pnm
    b_y_real = dPdt
    b_z_imag = m*beta

    dbdt_x_real = n*dPdt
    dbdt_y_real = d2Pdt2
    dbdt_z_imag = m*dbetadt

    R_b_real,R_b_imag = complex_multiply(R, b_x_real, zero(TF), b_y_real, zero(TF), zero(TF), b_z_imag)
    dRdt_b_real,dRdt_b_imag = complex_multiply(dRdt, b_x_real, zero(TF), b_y_real, zero(TF), zero(TF), b_z_imag)
    dRdp_b_real,dRdp_b_imag = complex_multiply(dRdp, b_x_real, zero(TF), b_y_real, zero(TF), zero(TF), b_z_imag)
    R_dbdt_real,R_dbdt_imag = complex_multiply(R, dbdt_x_real, zero(TF), dbdt_y_real, zero(TF), zero(TF), dbdt_z_imag)

    #u1x_real = n * Pnm
    #u1y_real = dPdt
    #u1z_imag = m * beta
    #dudt1x_real = n*dPdt
    #dudt1y_real = d2Pdt2
    #dudt1z_imag = m*dbetadt
    
    #vx_real, vx_imag = eimp_real * u1x_real, eimp_imag * u1x_real
    #vy_real, vy_imag = eimp_real * u1y_real, eimp_imag * u1y_real
    #vz_real, vz_imag = -eimp_imag * u1z_imag, eimp_real * u1z_imag

    #dvdt_x_real, dvdt_x_imag = eimp_real*dudt1x_real, eimp_imag*dudt1x_real 
    #dvdt_y_real, dvdt_y_imag = eimp_real*dudt1y_real, eimp_imag*dudt1y_real 
    #dvdt_z_real, dvdt_z_imag = eimp_real*dudt1z_imag, eimp_imag*dudt1z_imag

    #dvdp_x_real, dvdp_x_imag = -m * eimp_imag * u1x_real, m * eimp_real * u1x_real
    #dvdp_y_real, dvdp_y_imag = -m * eimp_imag * u1y_real, m * eimp_real * u1y_real
    #dvdp_z_real, dvdp_z_imag = -m * eimp_real * u1z_imag, -m * eimp_imag * u1z_imag
    
    # u_cartesian = R*eimp*u
    # dudt_cartesian = dRdt*eimp*u + R*eimp*dudt
    # dudp_cartesian = dRdp*eimp*u + R*im*eimp*u

	# transform to cartesian
	#u_real_cartesian, u_imag_cartesian = complex_multiply(R, vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag)
    #Rdvdt_real, Rdvdt_imag = complex_multiply(R, dvdt_x_real, dvdt_x_imag, dvdt_y_real, dvdt_y_imag, dvdt_z_real, dvdt_z_imag)
    #Rdvdp_real, Rdvdp_imag = complex_multiply(R, dvdp_x_real, dvdp_x_imag, dvdp_y_real, dvdp_y_imag, dvdp_z_real, dvdp_z_imag)
    #dRdt_v_real, dRdt_v_imag = complex_multiply(dRdt, u1x_real, 0, u1y_real, 0, 0, u1z_imag)
    #dRdp_v_real, dRdp_v_imag = complex_multiply(dRdp, u1x_real, 0, u1y_real, 0, 0, u1z_imag)
    #dudt_real_cartesian = dRdt_v_real + Rdvdt_real
    #dudt_imag_cartesian = dRdt_v_imag + Rdvdt_imag
    #dudp_real_cartesian = dRdp_v_real - Rdvdp_imag
    #dudp_imag_cartesian = dRdp_v_imag + Rdvdp_real
    # u_cartesian = Rv = eimp R u
    # dudt_cartesian = eimp (dRdt*u + R*dudt)
    # dudp_cartesian = (im*eimp*R + eimp*dRdp)*u
    # w = eimp*u3

    # r derivatives are the same for all cases, so just handle them all in one line.
    for i=1:3
        for j=1:3
            r_theta_phi_bar[j] += ḡ[j,i]*nm2_r*dudr[i]*R[i,j]
        end
    end

    #j=1
    
    # due to scalar potential
    # du/dr
    if PS
        Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(local_expansion[1,1,index],local_expansion[2,1,index],eimp_real,eimp_imag)
        #dudr -= nm1 * rnm2Cnm * SVector{3}(
		#	complex_multiply_real(u_real_cartesian[1], u_imag_cartesian[1], Lnm_real, Lnm_imag),
		#	complex_multiply_real(u_real_cartesian[2], u_imag_cartesian[2], Lnm_real, Lnm_imag),
		#	complex_multiply_real(u_real_cartesian[3], u_imag_cartesian[3], Lnm_real, Lnm_imag)
        #)
        for i=1:3
            r_theta_phi_bar[3] -= ḡ[1,i]*nm1*rnm2Cnm*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdt_b_real[i] + R_dbdt_real[i], dRdt_b_imag[i] + R_dbdt_imag[i])*sb2*R[i,1]
            r_theta_phi_bar[2] -= ḡ[1,i]*nm1*rnm2Cnm*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdp_b_real[i] - m*R_b_imag[i], dRdp_b_imag[i] + m*R_b_real[i])*sb2*R[i,1]
            local_expansion_bar[1,1,index] -= ḡ[1,i]*nm1*rnm2Cnm*complex_multiply_real(eimp_real,eimp_imag,R_b_real[i],R_b_imag[i])*sb2
            local_expansion_bar[2,1,index] += ḡ[1,i]*nm1*rnm2Cnm*complex_multiply_imag(eimp_real,eimp_imag,R_b_real[i],R_b_imag[i])*sb2
        end
    end
    
    # due to vector potential
    if VPS
		#dudr += nm1 * rnm2Cnm * complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
        for _p=1:3
            Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(local_expansion[1,_p+1,index],local_expansion[2,_p+1,index],eimp_real,eimp_imag)
            for k=1:3
                for i = 1:3
                    r_theta_phi_bar[3] += ḡ[1,i]*nm1*rnm2Cnm*ϵ(i,1,_p)*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdt_b_real[k] + R_dbdt_real[k], dRdt_b_imag[k] + R_dbdt_imag[k])*sb2*R[i,1]
                    r_theta_phi_bar[2] == ḡ[1,i]*nm1*rnm2Cnm*ϵ(i,1,_p)*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdp_b_real[k] - m*R_b_imag[k], dRdp_b_imag[k] + m*R_b_real[k])*sb2*R[i,1]
                    local_expansion_bar[1,_p+l,index] += ḡ[1,i]*nm1*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(eimp_real,eimp_imag,R_b_real[k],R_b_imag[k])*sb2
                    local_expansion_bar[2,_p+l,index] -= ḡ[1,i]*nm1*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_imag(eimp_real,eimp_imag,R_b_real[k],R_b_imag[k])*sb2
                end
            end
        end
    end

    #####
    ##### du/dphi / r
    #####
    
    # intermediate quantities

    c_x_real = nm1*dPdt
    c_y_real = d2Pdt2 + n*Pnm
    c_z_imag = m*ddt_Pnm_st

    dcdt_x_real = nm1*d2Pdt2
    dcdt_y_real = d3Pdt3 + n*dPdt
    dcdt_z_imag = m/st*(-ct/st*dPdt + d2Pdt2)

    R_c_real,R_c_imag = complex_multiply(R,c_x_real,zero(TF),c_y_real,zero(TF),zero(TF),c_z_imag)
    dRdt_c_real,dRdt_c_imag = complex_multiply(dRdt,c_x_real,zero(TF),c_y_real,zero(TF),zero(TF),c_z_imag)
    dRdp_c_real,dRdp_c_imag = complex_multiply(dRdp,c_x_real,zero(TF),c_y_real,zero(TF),zero(TF),c_z_imag)
    R_dcdt_real,R_dcdt_imag = complex_multiply(R,dcdt_x_real,zero(TF),dcdt_y_real,zero(TF),zero(TF),dcdt_z_imag)

    #u2x_real = nm1 * dPdt
    #u2y_real = d2Pdt2 + n * Pnm
    #u2z_imag = m * ddt_Pnm_st # double check that this is (1/sin θ)*dP/dt and not d/dt (P/sin(θ))

    #dudt2x_real = nm1*d2Pdt2
    #dudt2y_real = d3Pdt3 + n*dPdt
    #dudt2z_imag = m*(d2Pdt2/st - ct/st^2*dPdt)
    
    #vx_real, vx_imag = eimp_real * u2x_real, eimp_imag * u2x_real
    #vy_real, vy_imag = eimp_real * u2y_real, eimp_imag * u2y_real
    #vz_real, vz_imag = -eimp_imag * u2z_imag, eimp_real * u2z_imag

    #dvdt_x_real, dvdt_x_imag = eimp_real * dudt2x_real, eimp_imag * dudt2x_real
    #dvdt_y_real, dvdt_y_imag = eimp_real * dudt2y_real, eimp_imag * dudt2y_real
    #dvdt_z_real, dvdt_z_imag = -eimp_imag * dudt2z_imag, eimp_real * dudt2z_imag

    #dvdp_x_real, dvdp_x_imag = -m*eimp_imag * u2x_real, m*eimp_real * u2x_real
    #dvdp_y_real, dvdp_y_imag = -m*eimp_imag * u2y_real, m*eimp_real * u2y_real
    #dvdp_z_real, dvdp_z_imag = -m*eimp_real * u2z_imag, -m*eimp_imag * u2z_imag

	#u_real_cartesian, u_imag_cartesian = complex_multiply(R, vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag)
    #Rdvdt_real, Rdvdt_imag = complex_multiply(R, dvdt_x_real, dvdt_x_imag, dvdt_y_real, dvdt_y_imag, dvdt_z_real, dvdt_z_imag)
    #Rdvdp_real, Rdvdp_imag = complex_multiply(R, dvdp_x_real, dvdp_x_imag, dvdp_y_real, dvdp_y_imag, dvdp_z_real, dvdp_z_imag)
    #dRdt_v_real, dRdt_v_imag = complex_multiply(dRdt, u2x_real, 0, u2y_real, 0, 0, u2z_imag)
    #dRdp_v_real, dRdp_v_imag = complex_multiply(dRdp, u2x_real, 0, u2y_real, 0, 0, u2z_imag)
    #dudt_real_cartesian = dRdt_v_real + Rdvdt_real
    #dudt_imag_cartesian = dRdt_v_imag + Rdvdt_imag
    #dudp_real_cartesian = dRdp_v_real - Rdvdp_imag
    #dudp_imag_cartesian = dRdp_v_imag + Rdvdp_real

    #u_cross_L = complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
    #dudt_cross_L = complex_cross_real(dudt_real_cartesian[1], dudt_imag_cartesian[1], dudt_real_cartesian[2], dudt_imag_cartesian[2], dudt_real_cartesian[3], dudt_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
    #dudp_cross_L = complex_cross_real(dudp_real_cartesian[1], dudp_imag_cartesian[1], dudp_real_cartesian[2], dudp_imag_cartesian[2], dudp_real_cartesian[3], dudp_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)

    # du/dϕ / r

    # due to scalar potential
    if PS
        #dudt_r -= rnm2Cnm * SVector{3}(
		#	complex_multiply_real(u_real_cartesian[1], u_imag_cartesian[1], Lnm_real, Lnm_imag),
		#	complex_multiply_real(u_real_cartesian[2], u_imag_cartesian[2], Lnm_real, Lnm_imag),
		#	complex_multiply_real(u_real_cartesian[3], u_imag_cartesian[3], Lnm_real, Lnm_imag)
        #)
        Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(local_expansion[1,1,index],local_expansion[2,1,index],eimp_real,eimp_imag)
        for i=1:3
            r_theta_phi_bar[3] -= ḡ[2,i]*rnm2Cnm*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdt_c_real[i] + R_dcdt_real[i], dRdt_c_imag[i] + R_dcdt_imag[i])*sb2*R[i,2]
            r_theta_phi_bar[2] -= ḡ[2,i]*rnm2Cnm*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdp_c_real[i] - m*R_c_imag[i],dRdp_c_imag[i] + m*R_c_real[i])*sb2*R[i,2]
            local_expansion_bar[1,1,index] -= ḡ[2,i]*rnm2Cnm*complex_multiply_real(R_b_real[i],R_b_imag[i],eimp_real,eimp_imag)*sb2
            local_expansion_bar[2,1,index] += ḡ[2,i]*rnm2Cnm*complex_multiply_imag(R_b_real[i],R_b_imag[i],eimp_real,eimp_imag)*sb2
        end
    end
    
    # due to vector potential
    if VPS
		#dudt_r += rnm2Cnm * complex_cross_real(u_real_cartesian[1], u_imag_cartesian[1], u_real_cartesian[2], u_imag_cartesian[2], u_real_cartesian[3], u_imag_cartesian[3], Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
        for _p=1:3
            Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(local_expansion[1,_p+1,index],local_expansion[2,_p+1,index],eimp_real,eimp_imag)
            for k=1:3
                for i=1:3
                    r_theta_phi_bar[3] += ḡ[2,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdt_c_real[k] + R_dcdt_real[k],dRdt_c_imag[k] + R_dcdt_imag[k])*sb2*R[i,2]
                    r_theta_phi_bar[2] += ḡ[2,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dRdp_c_real[k] - m*R_c_imag[k],dRdp_c_imag[k] + m*R_c_real[k])*sb2*R[i,2]
                    local_expansion_bar[1,_p+1,index] += ḡ[2,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(eimp_real,eimp_imag,R_c_real[k],R_c_imag[k])*sb2
                    local_expansion_bar[2,_p+1,index] -= ḡ[2,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_imag(eimp_real,eimp_imag,R_c_real[k],R_c_imag[k])*sb2
                end
            end
        end
    end
    
    #####
    ##### du/dphi / r / sin(theta)
    #####
    
    # intermediate quantities

    Fnm = n * n * Pnm + d2Pdt2
    dFnmdt = n * n * dPdt + d3Pdt3
    m2m1 - m * m - 1

    #d_x_real = sp * Fnm
    #d_x_imag = m * cp * (alpha * m2m1 - Fnm)
    #d_y_real = -cp * Fnm
    #d_y_imag = m * sp * (alpha * m2m1 - Fnm)
    #d_z_real = 0.0
    #d_z_imag = m * (n * ct * beta - dPdt)

    d_real = [sp * Fnm, -cp * Fnm, zero(TF)]
    d_imag = [m * cp * (alpha * m2m1 - Fnm), m * sp * (alpha * m2m1 - Fnm), m * (n * ct * beta - dPdt)]

    #dddt_x_real = dFnmdt * st
    #dddt_x_imag = m * ct * (dalphadt * m2m1 - dFnmdt)
    #dddt_y_real = -dFnmdt * ct
    #dddt_y_imag = m * st * (dalphadt * m2m1 - dFnmdt)
    #dddt_z_real = 0.0
    #dddt_z_imag = m * (n * (dbetadt * ct - beta * st) - d2Pdt2)

    dddt_real = [dFnmdt * st, -dFnmdt * ct, zero(TF)]
    dddt_imag = [m * ct * (dalphadt * m2m1 - dFnmdt), m * st * (dalphadt * m2m1 - dFnmdt), m * (n * (dbetadt * ct - beta * st) - d2Pdt2)]

    #dddp_x_real = Fnm * ct
    #dddp_x_imag = -m * st * (alpha * m2m1 - Fnm)
    #dddp_y_real = Fnm * st
    #dddp_y_imag = m * ct * (alpha * m2m1 - Fnm)
    #dddp_z_real = 0.0
    #dddp_z_imag = 0.0

    dddp_real = [Fnm * ct, Fnm * st, zero(TF)]
    dddp_imag = [-m * st * (alpha * m2m1 - Fnm), m * ct * (alpha * m2m1 - Fnm), zero(TF)]

    #R_d_real,R_d_imag = complex_multiply(R,d_x_real,d_x_imag,d_y_real,d_y_imag,0.0,d_z_imag)
    #dRdt_d_real,dRdt_d_imag = complex_multiply(dRdt,d_x_real,d_x_imag,d_y_real,d_y_imag,0.0,d_z_imag)
    #dRdp_d_real,dRdp_d_imag = complex_multiply(dRdp,d_x_real,d_x_imag,d_y_real,d_y_imag,0.0,d_z_imag)
    #R_dddt_real,R_dddt_imag = complex_multiply(R,dddt_x_real,dddt_x_imag,dddt_y_real,dddt_y_imag,0.0,dddt_z_imag)
    #R_dddp_real,R_dddp_imag = complex_multiply(R,dddp_x_real,dddp_x_imag,dddp_y_real,dddp_y_imag,0.0,0.0)

    # du/dϕ / r / sin(θ)

    # due to scalar potential
    if PS
        #dudp_r_st -= SVector{3}(
        #    complex_multiply_real(vx_real, vx_imag, Lnm_real, Lnm_imag),
        #    complex_multiply_real(vy_real, vy_imag, Lnm_real, Lnm_imag),
		#	complex_multiply_real(vz_real, vz_imag, Lnm_real, Lnm_imag)
        #) * rnm2Cnm
        Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(local_expansion[1,1,index],local_expansion[2,1,index],eimp_real,eimp_imag)
        for i=1:3
            r_theta_phi_bar[3] -= ḡ[3,i]*rnm2Cnm*complex_multiply_real(Lnm_eimp_real, Lnm_eimp_imag, dddt_real[i],dddt_imag[i])*sb2*R[i,3]
            r_theta_phi_bar[2] -= ḡ[3,i]*rnm2Cnm*complex_multiply_real(Lnm_eimp_real, Lnm_eimp_imag, dddp_real[i] - m*d_imag[i],dddp_imag[i] + m*d_real[i])*sb2*R[i,3]
            local_expansion_bar[1,1,index] -= ḡ[3,i]*rnm2Cnm*complex_multiply_real(eimp_real,eimp_imag,eimp_imag,d_real[i],d_imag[i])*sb2
            local_expansion_bar[2,1,index] += ḡ[3,i]*rnm2Cnm*complex_multiply_imag(eimp_real,eimp_imag,eimp_imag,d_real[i],d_imag[i])*sb2
        end
    end

    # due to vector potential
    if VPS
        #dudp_r_st += rnm2Cnm * complex_cross_real(vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag, Lnm_x_real, Lnm_x_imag, Lnm_y_real, Lnm_y_imag, Lnm_z_real, Lnm_z_imag)
        for i=1:3
            for j=1:3
                r_theta_phi_bar[1] += ḡ[j,i]*(n-2)*rnm3Cnm*v_cross_L[i]*R[j,k]*sb2
                r_theta_phi_bar[2] += ḡ[j,i]*rnm2Cnm*dvdt_cross_L[i]*R[j,k]*sb2
                r_theta_phi_bar[3] += ḡ[j,i]*rnm2Cnm*dvdp_cross_L[i]*R[j,k]*sb2
            end
        end
        for _p=1:3
            Lnm_eimp_real,Lnm_eimp_imag = complex_multiply(local_expansion[1,_p+1,index],local_expansion[2,_p+1,index],eimp_real,eimp_imag)
            for k=1:3
                for i=1:3
                    r_theta_phi_bar[3] += ḡ[3,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dddt_real[k],dddt_imag[k])*sb2*R[i,3]
                    r_theta_phi_bar[2] += ḡ[3,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(Lnm_eimp_real,Lnm_eimp_imag,dddp_real[k] - m*d_imag[k], ddp_imag[k] + m*d_real[k])*sb2*R[i,3]
                    local_expansion_bar[1,_p+1,index] += ḡ[3,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_real(eimp_real,eimp_imag,d_real[k],d_imag[k])*sb2
                    local_expansion_bar[2,_p+1,index] -= ḡ[3,i]*rnm2Cnm*ϵ(i,k,_p)*complex_multiply_imag(eimp_real,eimp_imag,d_real[k],d_imag[k])*sb2
                end
            end
        end
    end

	#return dudr, dudt_r, dudp_r_st
    # updated in-place by the end: local_expansion_bar, r_theta_phi_bar
    return nothing
end

function C2S_pullback!(body_position_bar, expansion_center_bar, r_theta_phi_bar, dx, dy, dz, r, theta, phi)

    dC2S = [dx/r dy/r dz/r; dx*cot(theta)/r^2 dy*cot(theta)/r^2 -1/(r*sin(theta)) + dz*cot(theta)/r^2; -sin(phi)*(cos(phi)/dx) cos(phi)^2/dx 0]
    for j=1:3
        body_position_bar[j] = dC2S[1,j]*r_theta_phi_bar[1] + dC2S[2,j]*r_theta_phi_bar[2] + dC2S[3,j]*r_theta_phi_bar[3]
        expansion_center_bar[j] = -body_position_bar[j]
    end

    return nothing
end

# TODO: Make the rho/phi/theta order consistent, if possible. This might be a change in spherical.jl.
# TODO: Double-check the recursion relations; something might be wrong for higher expansion orders
function ChainRulesCore.rrule(::typeof(L2B),body_position, expansion_center::AbstractArray{TF}, local_expansion, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {TF,PS,VPS,VS,GS,P}

    scalar_potential, vector_potential, velocity, gradient = L2B(body_position, expansion_center, local_expansion, derivatives_switch, expansion_order)

    function pullbacks(ȳ)
        #println("running L2B pullback:")
        #@time begin
        #s̄p,v̄p,v̄,ḡ = ȳ # I don't actually want to allocate here... would a view be better?
        s̄p = ȳ[1] # scalar, so no view
        v̄p = view(ȳ[2],:)
        v̄ = view(ȳ[3],:)
        ḡ = reshape(view(ȳ[4],:),(3,3))

        # initialize cotangent vectors
        body_position_bar = zeros(length(body_position))
        expansion_center_bar = zeros(length(expansion_center))
        local_expansion_bar = zeros(size(local_expansion))

        # polar coordinate cotangents
        r_theta_phi_bar = zeros(TF,3) # [r, theta, phi] cotangents
        #duidrk_bar = zeros(TF,(3,3)) # spare memory for calculations

        #=
        R_n^m(\rho, \theta, \phi) &= (-1)^n \frac{\rho^n}{(n+\|m\|)!} P_n^{\|m\|}(\cos \theta) e^{i m \phi}\\
        I_n^m(\rho, \theta, \phi) &= (-1)^n \frac{(n-\|m\|)!}{\rho^{n+1}} P_n^{\|m\|}(\cos \theta) e^{i m \phi}\\
        phi = \sum \limits_{n=0}^P \sum \limits_{m=-n}^n L_n^m R_n^m
        =#
        # still need these intermediate polar coordinates
        dx, dy, dz = body_position - expansion_center
        r, theta, phi = cartesian_2_spherical(dx, dy, dz)
        
        #--- the following is inspired by Salloum and Lakkis (2020) with some improvements ---#

        # initialization # not needed - we already have scalar_potential, vector_potential, velocity, and gradient.
        #scalar_potential = zero(TF)
        #vector_potential = zero(SVector{3,TF})
        #velocity = zero(SVector{3,TF})
        #gradient = zero(SMatrix{3,3,TF,9})
        #dudr = zero(SVector{3,TF})
        #dudt_r = zero(SVector{3,TF})
        #dudp_r_st = zero(SVector{3,TF})

        # more intermediate values
        st, ct = sincos(theta)
        sp, cp = sincos(phi)

        if VS || GS
            R = SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
            # define two derivative matrices (dR/dθ and dR/dϕ) that will be needed later
            dRdt = SMatrix{3,3}(ct*cp,ct*sp,-st,-st*cp,-st*sp,-ct,0,0,0)
            dRdp = SMatrix{3,3}(-st*sp,st*cp,0,-ct*sp,ct*cp,0,-cp,-sp,0)
        end

        # note that by definition beta_n^m=0 for m<1 and alpha_n^m=0 for m<2
        
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
        dPdt_n_m = 0.0      # derivative w.r.t. theta of P_n_m
        d2Pdt2_n_m = 0.0    # second derivative w.r.t. theta of P_n_m
        index = 1           # index of the local expansion corresponding to n and m

        if PS
            #scalar_potential += L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            # I'm building the pullback functions with signatures of the form f!(cotangents,original inputs, other required inputs)
            scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, s̄p, r, theta, phi, n, m, dPdt_n_m)
            # scalar_potential += local_expansion[1,1,index]
        end

        if VPS
            #vector_potential += L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
            vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, v̄p, r, theta, phi, n, m, dPdt_n_m)
            # vector_potential += SVector{3}(
            #     local_expansion[1,2,1],
            #     local_expansion[1,3,1],
            #     local_expansion[1,4,1]
            # )
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
            dbetadt_n_m = zero(TF)
            
            # index
            index += 1

            #####
            ##### evaluate
            #####

            if PS
                #scalar_potential += L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, s̄p, r, theta, phi, n, m, dPdt_n_m)
            end

            if VPS
                #vector_potential += L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, v̄p, r, theta, phi, n, m, dPdt_n_m)
            end

            if VS
                #velocity += L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch, index, v̄, r, theta, phi, d2Pdt2_n_m, dbetadt_n_m, dRdt, dRdp, velocity)
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

            dbetadt_n_m = 0.0
            dbetadt_nm1_m = 0.0
            dbetadt_nm2_m = 0.0

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

            if PS
                #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, s̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
            end

            if VPS
                #vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, v̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
            end

            if VS
                #velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch, index, v̄, r, theta, phi, d2Pdt2_n_m, dbetadt_n_m, dRdt, dRdp, velocity; scale_by_two=true)
            end

        end
        
        #--- m = 0, m = 1, n > 1 ---#

        if P > 1
            
            m = 1
            alpha_n_0 = zero(TF)        # n>1, m<2
            beta_n_0 = zero(TF)         # n>1, m<1
            dbetadt_n_0 = zero(TF)      # n>1, m<1
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

            # additional tracked terms:
            alpha_nm2_m = 0.0
            alpha_nm1_m = 0.0
            alpha_n_m = 3.0 # alpha_2_2 = 3.0
            dalphadt_nm2_m = 0.0
            dalphadt_nm1_m = 0.0
            dalphadt_n_m = 0.0
            dalphadt_n_0 = 0.0 # alpha_0_0 is zero, so all alpha_n_0 are zero

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

                # dbetadt
                dbetadt_nm2_m = dbetadt_nm1_m
                dbetadt_nm1_m = dbetadt_n_m
                dbetadt_n_m = ((2*n-1)*(-st*beta_nm1_m + ct*dbetadt_nm1_m) - n*dbetadt_nm2_m)/(n-1)
                
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

                # third derivative for m=1
                alpha_nm2_m = alpha_nm1_m
                alpha_nm1_m = alpha_n_m
                alpha_n_m = ((2*n-1)*ct*alpha_nm1_m - (n+m-1)*alpha_nm2_m)/(n-m)
                dalphadt_nm2_m = dalphadt_nm1_m
                dalphadt_nm1_m = dalphadt_n_m
                dalphadt_n_m = ((2*n-1)*(st*alpha_nm1_m + ct*dalphadt_nm1_m) - (n+m-1)*dalphadt_nm2_m)/(n-m)
                d3Pdt3_n_m = (n+m)*(-st*alpha_nm1_m + ct*dalphadt_nm1_m) - 2*n*n*st*ct*alpha_n_m - (n*(n*st*st+1)-m*m)*dalphadt_n_m

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
                    #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                    scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, one(TF), zero(TF), rn, P_n_0, C_n_0, index_m0, s̄p, r, theta, phi, n, 0, dPdt_n_0)
                    scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, s̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
                end
        
                if VPS
                    #vector_potential += L2B_vector_potential(local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rn, P_n_0, C_n_0)
                    #vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, Pnm, C_n_m)
                    vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, one(TF), zero(TF), rn, P_n_0, C_n_0, index_m0, v̄p, r, theta, phi, n, 0, dPdt_n_0)
                    vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, Pnm, C_n_m, index, v̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
                end
        
                if VS
                    #velocity += L2B_velocity(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm1, P_n_0, dPdt_n_0, beta_n_0, C_n_0, R, n, 0, derivatives_switch)
                    #velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                    velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, one(TF), zero(TF), rnm1, P_n_0, dPdt_n_0, beta_n_0, C_n_0, R, n, 0, derivatives_switch, index_m0, v̄, r, theta, phi, d2Pdt2_n_0, dbetadt_n_0, dRdt, dRdp, velocity)
                    velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm1, Pnm, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch, index, v̄, r, theta, phi, d2Pdt2_n_m, dbetadt_n_m, dRdt, dRdp, velocity; scale_by_two=true)
                end

                if GS
                    # m = 0
                    ddt_Pn0_st = zero(TF)
                    #this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index_m0], local_expansion[2,1,index_m0], local_expansion[1,2,index_m0], local_expansion[2,2,index_m0], local_expansion[1,3,index_m0], local_expansion[2,3,index_m0], local_expansion[1,4,index_m0], local_expansion[2,4,index_m0], one(TF), zero(TF), rnm2, P_n_0, dPdt_n_0, d2Pdt2_n_0, ddt_Pn0_st, alpha_n_0, beta_n_0, st, ct, sp, cp, C_n_0, R, n, 0, derivatives_switch)
                    gradient_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, one(TF), zero(TF), rnm2, P_n_0, dPdt_n_0, d2Pdt2_n_0, ddt_Pn0_st, alpha_n_0, beta_n_0, st, ct, sp, cp, C_n_0, R, n, 0, derivatives_switch, index_m0, ḡ, r, dRdt, dRdp, dalphadt_n_0, dbetadt_n_0, d3Pdt3_n_0)
                    #dudr += this_dudr
                    #dudt_r += this_dudt_r
                    #dudp_r_st += this_dudp_r_st
                    
                    # m = 1

                    #this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, Pnm, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_0, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)
                    gradient_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm2, Pnm, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_0, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch, index, ḡ, r, dRdt, dRdp, dalphadt_n_m, dbetadt_n_m, d3Pdt3_n_m; scale_by_two=true)
                    #dudr += 2 * this_dudr
                    #dudt_r += 2 * this_dudt_r
                    #dudp_r_st += 2 * this_dudp_r_st
                end

            end

            #--- m > 1, n > 1 ---#

            # double factorial
            # double_factorial = 3

            # sin^(n-2)(theta)
            # st_nm2 = one(TF)

            # alpha_n_n
            alpha_n_n = TF(3.0) # sin(theta)^(n-2) (2*n-1)!! for n=2
            dalphadt_n_n = TF(0.0)

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

                # dalphadt
                dalphadt_n_m = dalphadt_n_n
                dalphadt_nm1_m = 0.0
                dalphadt_nm2_m = 0.0

                # dbetadt

                dbetadt_n_m = dalphadt_n_m*st + alpha_n_m*ct
                #dbetadt_nm1_m = 0.0 # not actually needed.

                # third derivative
                d3Pdt3_n_m = (n+m)*(-st*alpha_nm1_m + ct*dalphadt_nm1_m) - 2*n^2*st*ct*alpha_n_m - (n*(n*st^2+1) - m^2) * dalphadt_n_m
                
                #####
                ##### evaluate
                #####

                if PS
                    #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
                    scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, P_n_m, C_n_m, index, s̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
                end
        
                if VPS
                    #vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn, P_n_m, C_n_m)
                    vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn, P_n_m, C_n_m, index, v̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
                end
        
                if VS
                    #velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1, P_n_m, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch)
                    velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm1, P_n_m, dPdt_n_m, beta_n_m, C_n_m, R, n, m, derivatives_switch, index, v̄, r, theta, phi, d2Pdt2_n_m, dbetadt_n_m, dRdt, dRdp, velocity; scale_by_two=true)
                end

                if GS
                    #this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch)
                    gradient_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm2, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m, R, n, m, derivatives_switch, index, ḡ, r, dRdt, dRdp, dalphadt_n_m, dbetadt_n_m, d3Pdt3_n_m; scale_by_two=true)
                    #dudr += 2 * this_dudr
                    #dudt_r += 2 * this_dudt_r
                    #dudp_r_st += 2 * this_dudp_r_st
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

                    # dalphadt

                    dalphadt_nm2_m = dalphadt_nm1_m
                    dalphadt_nm1_m = dalphadt_n_m
                    dalphadt_n_m = ((2*n-1)*(st*alpha_nm1_m + ct*dalphadt_nm1_m) - (n+m-1)*dalphadt_nm2_m)/(n-m)

                    # dbetadt

                    dbetadt_n_m = dalphadt_n_m*st + alpha_n_m*ct
                    
                    # third derivative

                    d3Pdt3_n_m = (n+m)*(-st*alpha_nm1_m + ct*dalphadt_nm1_m) - 2 * n^2 * st * ct * alpha_n_m - (n*(n*st^2 + 1) - m^2) * dalphadt_n_m

                    #####
                    ##### evaluate
                    #####
                    if PS
                        #scalar_potential += 2 * L2B_scalar_potential(local_expansion[1,1,index], local_expansion[2,1,index], eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop)
                        scalar_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop, index, s̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
                    end
            
                    if VPS
                        #vector_potential += 2 * L2B_vector_potential(local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop)
                        vector_potential_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rn_loop, P_n_m, C_n_m_loop, index, v̄p, r, theta, phi, n, m, dPdt_n_m; scale_by_two=true)
                    end
            
                    if VS
                        #velocity += 2 * L2B_velocity(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm1_loop, P_n_m, dPdt_n_m, beta_n_m, C_n_m_loop, R, n, m, derivatives_switch)
                        velocity_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm1_loop, P_n_m, dPdt_n_m, beta_n_m, C_n_m_loop, R, n, m, derivatives_switch, index, v̄, r, theta, phi, d2Pdt2_n_m, dbetadt_n_m, dRdt, dRdp, velocity; scale_by_two=true)
                    end
        
                    if GS
                        #this_dudr, this_dudt_r, this_dudp_r_st = L2B_velocity_gradient(local_expansion[1,1,index], local_expansion[2,1,index], local_expansion[1,2,index], local_expansion[2,2,index], local_expansion[1,3,index], local_expansion[2,3,index], local_expansion[1,4,index], local_expansion[2,4,index], eimp_real, eimp_imag, rnm2_loop, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m_loop, R, n, m, derivatives_switch)
                        gradient_pullback!(local_expansion_bar, r_theta_phi_bar, local_expansion, eimp_real, eimp_imag, rnm2_loop, P_n_m, dPdt_n_m, d2Pdt2_n_m, ddt_Pnm_st, alpha_n_m, beta_n_m, st, ct, sp, cp, C_n_m_loop, R, n, m, derivatives_switch, index, ḡ, r, dRdt, dRdp, dalphadt_n_m, dbetadt_n_m, d3Pdt3_n_m; scale_by_two=true)
                        #dudr += 2 * this_dudr
                        #dudt_r += 2 * this_dudt_r
                        #dudp_r_st += 2 * this_dudp_r_st
                    end
                end

                # tail recursion
                dalphadt_n_n = -(2*n+1)*(ct*alpha_n_n + st*dalphadt_n_n)
                alpha_n_n *= -st * (2*n+1)
                # double_factorial *= 2*n+1
                # st_nm2 *= -st

            end

        end

        # rotate to cartesian
        # R = SMatrix{3,3}(st*cp,st*sp,ct,ct*cp,ct*sp,-st,-sp,cp,0)
        #duidrk = hcat(dudr, dudt_r, dudp_r_st)
        #gradient = duidrk * R'

        # extra term from product rule due to rotation:
        # r_theta_phi_bar[l] += duidrk[i,k]*dRd(r[l])[j,k]
        # also note that duidrk = gradient * R (since R' = R^-1)
        # so: r_theta_phi_bar[l] += gradient[i,ii]*R[ii,k]*dRd(r[l])[j,k]
        # manually iterate over l, since that corresponds to iterating over specific polar coordinates.
        # dRdt and dRdp are previously calculated. dRdr = 0.
        if VS || GS
            for i=1:3
                for j=1:3
                    # no contribution to r̄ since dRdr = 0
                    r_theta_phi_bar[3] += ḡ[j,i]*gradient[i,j]*R[j,i]*dRdt[i,j]
                    r_theta_phi_bar[2] += ḡ[j,i]*gradient[i,j]*R[j,i]*dRdp[i,j]
                end
            end
        end
        
        #return scalar_potential, vector_potential, velocity, gradient
        # account for cartesian coordinate inputs
        C2S_pullback!(body_position_bar, expansion_center_bar, r_theta_phi_bar, dx, dy, dz, r, theta, phi)
        #end
        #@show body_position_bar expansion_center_bar local_expansion_bar
        #@show r_theta_phi_bar
        return NoTangent(),body_position_bar,expansion_center_bar,local_expansion_bar,NoTangent(),NoTangent()
    end

    return (scalar_potential, vector_potential[1:3], velocity[1:3], gradient[1:3,1:3]), pullbacks
end
#=
@grad_from_chainrules_multiple_returns L2B(body_position::AbstractArray{<:ReverseDiff.TrackedReal},
                                   expansion_center::AbstractArray{<:ReverseDiff.TrackedReal},
                                   local_expansion::AbstractArray{<:ReverseDiff.TrackedReal},
                                   derivatives_switch::DerivativesSwitch,
                                   expansion_order::Val)=#