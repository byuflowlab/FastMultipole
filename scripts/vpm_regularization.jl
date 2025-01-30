using Roots
using SpecialFunctions: erf
using PythonPlot
using LaTeXStrings

# function myerf(x)
#     a1 = 0.278393
#     a2 = 0.230389
#     a3 = 2/sqrt(π)
#     return a3 * (x + a1 * x * x * x / (1 + a2 * x * x))
# end

function myerf(x)
    return tanh(x*pi/sqrt(6))
end

function upper_bound(σ, ω, ε)
    return ω / (8 * pi * ε * σ) * (sqrt(2/pi) + sqrt(2/(pi*σ*σ) + 16 * pi * ε / ω))
end

function residual(ρ_σ, σ, ω, ε)
    t1 = 4*pi*σ*σ*ε*ρ_σ*ρ_σ / ω
    # t2 = erf(ρ_σ / sqrt(2))
    t2 = myerf(ρ_σ / sqrt(2))
    t3 = sqrt(2/pi) * ρ_σ * exp(-ρ_σ*ρ_σ*0.5)
    return t1 + t2 - t3 - 1.0
end

function solve_ρ_over_σ(σ, ω, ε)
    return Roots.find_zero((x) -> residual(x, σ, ω, ε), (0.0, upper_bound(σ, ω, ε)), Roots.Brent())
end

ω = 1.0
σs = [10.0^(-i) for i in 0:5]
εs = [10.0^(-i) for i in range(0,stop=9,length=50)]

res = zeros(length(εs), length(σs))

for (iσ,σ) in enumerate(σs)
    for (iε,ε) in enumerate(εs)
        res[iε,iσ] = solve_ρ_over_σ(σ, ω, ε)
    end
end

fig = figure("regularization with myerf")
fig.clear()
fig.add_subplot(111, xlabel=L"\varepsilon_{tol}", ylabel=L"\rho / \sigma")
ax = fig.get_axes()[0]
ncolors = length(σs)
colormap = get_cmap("RdBu",ncolors)
for iσ in 1:length(σs)
    ax.plot(εs, res[:,iσ], label="σ = $(round(σs[iσ], sigdigits=1))", color=colormap(iσ-1))
end
ax.legend()
ax.set_xscale("log")
