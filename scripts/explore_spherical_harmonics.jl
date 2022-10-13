using LegendrePolynomials
import FLOWFMM as fmm
import PyPlot as plt
using LaTeXStrings

function Ylm(theta, phi, l, m)
    ylm = (-1)^m * sqrt(factorial(big(l-abs(m)))/ factorial(big(l+abs(m)))) * Plm(cos(theta), l, abs(m)) * exp(im * m * phi)
    return ylm
end

function Ylms(theta, phi, p)
    ylms = Vector{Complex{Float64}}(undef,((p+1)*(p+2)) >> 1 + 1)
    i = 1
    for l in 0:p
        for m in -l:l
            ylms[i] = Ylm(theta, phi, l, m)
            i+=1
        end
    end
    return ylms
end

function evaluate_biot_savart(x_source, x_target, P; debug=true)
    if debug; (terms = zeros((P+1)^2)); end
    v = 0.0
    i = 1
    for l in 0:P
        for m in -l:l
            if debug; terms[i] = x_source[1]^l / x_target[1]^(l+1) * real(Ylm(x_target[2], x_target[3], l, m) * conj(Ylm(x_source[2], x_source[3], l, m))); end
            v += x_source[1]^l / x_target[1]^(l+1) * real(Ylm(x_target[2], x_target[3], l, m) * conj(Ylm(x_source[2], x_source[3], l, m)))
            i += 1
        end
    end
    debug && return v, terms
    return v
end

function test_spherical_harmonics(spacing; Ps = 0:1:20, x_source = [0.2,0.5,-0.1])
    x_target = [-10.1, -11.4, 3.5]
    x_target .*= sqrt(x_source' * x_source) / sqrt(x_target' * x_target) * spacing
    dx = sqrt(sum((x_source - x_target) .^2 ))

    one_over_r = 1/dx

    x_source_spherical = fmm.cartesian_2_spherical(x_source)
    x_target_spherical = fmm.cartesian_2_spherical(x_target)

    vals = [evaluate_biot_savart(x_source, x_target, P; debug=false) for P in Ps]
    return Ps, vals, one_over_r
end

function sweep_tsh(spacings; Ps=0:1:20, x_source = [0.2,0.5,-0.1])
    vals_vec = Vector{Vector{Float64}}(undef,length(spacings))
    for i in 1:length(spacings)
        Ps, vals, one_over_r = test_spherical_harmonics(spacings[i]; Ps, x_source)
        vals_vec[i] = abs.(vals .- one_over_r) ./ one_over_r
    end
    return Ps, vals_vec
end

spacings = 1:7
Ps, vals_vec = sweep_tsh(spacings)

fig = plt.figure("test")
fig.clear()
ax = fig.add_subplot(111, xlabel="expansion order", ylabel="relative error")
for (i,vals) in enumerate(vals_vec)
    ax.plot(Ps, vals, label=L"r_t/r_s = " * "$(spacings[i])")
end
# ax.set_yscale("log")
ax.legend()
fig.savefig("test_biot_savart_kernel.png")
