# complex arithmetic functions to improve the efficiency of algorithmic differentiation with ReverseDiff.jl

@inline complex_add(z1_real, z1_imag, z2_real, z2_imag) = z1_real + z2_real, z1_imag + z2_imag

@inline complex_subtract(z1_real, z1_imag, z2_real, z2_imag) = z1_real - z2_real, z1_imag - z2_imag

@inline complex_multiply_real(z1_real, z1_imag, z2_real, z2_imag) = z1_real * z2_real - z1_imag * z2_imag
@inline complex_multiply_imag(z1_real, z1_imag, z2_real, z2_imag) = z1_real * z2_imag + z1_imag * z2_real
@inline complex_multiply(z1_real, z1_imag, z2_real, z2_imag) = complex_multiply_real(z1_real, z1_imag, z2_real, z2_imag), complex_multiply_imag(z1_real, z1_imag, z2_real, z2_imag)

@inline complex_multiply_real(z1r,z1i,z2r,z2i,z3r,z3i) = (z1r*z2r - z1i*z2i)*z3r - (z1r*z2i + z1i*z2r)*z3i
@inline complex_multiply_imag(z1r,z1i,z2r,z2i,z3r,z3i) = (z1r*z2i + z1i*z2r)*z3r + (z1r*z2r - z1i*z2i)*z3i
@inline complex_multiply(z1r,z1i,z2r,z2i,z3r,z3i) = complex_multiply_real(z1r,z1i,z2r,z2i,z3r,z3i), complex_multiply_imag(z1r,z1i,z2r,z2i,z3r,z3i)

@inline complex_divide_real(z1_real, z1_imag, z2_real, z2_imag) = (z1_real * z2_real + z1_imag * z2_imag) / (z2_real^2 + z2_imag^2)
@inline complex_divide_imag(z1_real, z1_imag, z2_real, z2_imag) = (z1_imag * z2_real - z1_real * z2_imag) / (z2_real^2 + z2_imag^2)
@inline function complex_divide(z1_real, z1_imag, z2_real, z2_imag)
    denom = z2_real^2 + z2_imag^2
    num_real = z1_real * z2_real + z1_imag * z2_imag
    num_imag = z1_imag * z2_real - z1_real * z2_imag
    return num_real / denom, num_imag / denom
end

# vector arithmetic

@inline function complex_cross_real(ux_real, ux_imag, uy_real, uy_imag, uz_real, uz_imag, vx_real, vx_imag, vy_real, vy_imag, vz_real, vz_imag)
    return SVector{3}(
        complex_multiply_real(uy_real, uy_imag, vz_real, vz_imag) - complex_multiply_real(vy_real, vy_imag, uz_real, uz_imag),
        complex_multiply_real(vx_real, vx_imag, uz_real, uz_imag) - complex_multiply_real(ux_real, ux_imag, vz_real, vz_imag),
        complex_multiply_real(ux_real, ux_imag, vy_real, vy_imag) - complex_multiply_real(vx_real, vx_imag, uy_real, uy_imag)
    )
end
