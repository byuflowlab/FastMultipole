function M2M!(branch, child, harmonics, M, expansion_order::Val{P}, ::Type{LambHelmholtz}) where P
    # get distance vector
    dx, dy, dz = branch.center - child.center
    r, theta, phi = cartesian_2_spherical(dx, dy, dz)
    regular_harmonic!(harmonics, r, theta, phi, P)

    for j in 0:P # iterate over new Multipole coefficients B_j^k
        for k in 0:j
            i_jk = ((j * (j+1)) >> 1) + k + 1 # current index
            M .= zero(eltype(M))

            for l in 0:j
                for m in max(-l,-j+k+l):min(k-1,l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) + k - m + 1
                    lm = l * l + l - m + 1
                    ipow = ipow2l(m)
                    oddeven = odd_or_even(l)
                    #C_tmp = harmonics[lm] * ipow * oddeven
                    C_tmp1 = harmonics[1,lm] * ipow * oddeven
                    C_tmp2 = harmonics[2,lm] * ipow * oddeven
                    for dim in 1:1
                        #@inbounds M[dim] += child.multipole_expansion[dim,jlkms] * C_tmp
                        @inbounds M[1,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp1 - child.multipole_expansion[2,dim,jlkms] * C_tmp2
                        @inbounds M[2,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp2 + child.multipole_expansion[2,dim,jlkms] * C_tmp1
                    end
                end

                for m in k:min(l,j+k-l)
                    jlkms = (((j-l) * (j-l+1)) >> 1) - k + m + 1
                    lm = l * l + l - m + 1
                    oddeven = odd_or_even(k + l + m)
                    #C_tmp = harmonics[lm] * oddeven
                    C_tmp1 = harmonics[1,lm] * oddeven
                    C_tmp2 = harmonics[2,lm] * oddeven
                    for dim in 1:1
                        #@inbounds M[dim] += conj(child.multipole_expansion[dim,jlkms]) * C_tmp
                        @inbounds M[1,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp1 + child.multipole_expansion[2,dim,jlkms] * C_tmp2
                        @inbounds M[2,dim] += child.multipole_expansion[1,dim,jlkms] * C_tmp2 - child.multipole_expansion[2,dim,jlkms] * C_tmp1
                    end
                end
            end

            for dim in 1:1
                #@inbounds branch.multipole_expansion[dim,i_jk] += M[dim]
                @inbounds branch.multipole_expansion[1,dim,i_jk] += M[1,dim]
                @inbounds branch.multipole_expansion[2,dim,i_jk] += M[2,dim]
            end
        end
    end
end
