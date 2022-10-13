using LegendrePolynomials
import FLOWFMM as fmm
import PyPlot as plt
using LaTeXStrings

function Ylm(theta, phi, l, m)
    plm = Plm(cos(theta), l, abs(m))
    ylm = sqrt(factorial(big(l-abs(m)))/ factorial(big(l+abs(m)))) * plm * exp(im * m * phi)
    return ylm
end

function collect_ylm(rho, theta, phi, p)
    ylms = zeros(Complex{Float64}, (p+1)^2)
    i = 1
    for l in 0:p
        for m in -l:l
            ylms[i] = Ylm(theta, phi, l, m) * rho^l
            i += 1
        end
    end
    return ylms
end

function collect_local(rho, theta, phi, p)
    local_expansion = zeros(Complex{Float64}, (p+1)^2)
    i = 1
    for l in 0:p
        for m in -l:l
            # note: converges far from source locations
            local_expansion[i] = Ylm(theta, phi, l, m) / rho^(l+1)
            i += 1
        end
    end
    return local_expansion
end

function evaluate_biot_savart(x_source_old, x_target_old, P)
    x_source = fmm.cartesian_2_spherical(x_source_old)
    x_target = fmm.cartesian_2_spherical(x_target_old)
    ylms_target = collect_ylm(x_target[1], x_target[2], x_target[3], P)
    ylms_source = collect_local(x_source[1], x_source[2], -x_source[3], P) # get complex conjugate
    i = 1
    v = 0.0
    # r_ratio = x_target[1] / x_source[1]
    for l in 0:P
        for m in -l:l
            # v += real(-Ylm(x_source[2], x_source[3], l, -m) / x_source[1]^(l+1) * Ylm(x_target[2], x_target[3], l, m) * x_target[1]^l)
            # r = r_ratio^l / x_source[1]
            v += real(ylms_target[i] * ylms_source[i])
            # if imag(ylms_target[i] * ylms_source[i]) > 1e-12
            #     println("Warning: imaginary part of expansion is $(imag(ylms_target[i] * ylms_source[i]))")
            # end
            i += 1
        end
    end
    return v, ylms_target, ylms_source
end

function test_spherical_harmonics(x_target, x_source, P::Real)
    v_local = evaluate_biot_savart(x_source, x_target, P)
    return v_local
end

function test_spherical_harmonics(x_target, x_source, P::Vector{TR}) where TR <: Real
    dx = x_target - x_source
    v_true = 1/sqrt(dx' * dx)
    vs = zeros(length(P))
    for (i,p) in enumerate(P)
        vs[i] = test_spherical_harmonics(x_target, x_source, p)
    end
    return v_true, vs
end

# x_target = [0.5,0.2,0.1]
# x_source = 10*[0.5,0.2,0.1]
# Ps = collect(1:20)
# v_true, vs = test_spherical_harmonics(x_target, x_source, Ps)
# fig = plt.figure("test_local")
# fig.clear()
# ax = fig.add_subplot(111, xlabel="expansion order", ylabel="relative error")
# ax.plot(Ps, abs.(vs .- v_true) ./ v_true)#, label=L"r_t/r_s = " * "$(spacings[i])")
# ax.set_yscale("log")
# # ax.legend()
# fig.savefig("test_local.png")

x_target = [-0.1, 0.7, 0.4]
x_source = 10 * x_target
P = 3
v, ylms_target, ylms_source = evaluate_biot_savart(x_source, x_target, P)
