using FastMultipole
using PythonPlot

include("rotor.jl")

fig = figure("test")
fig.clear()
for i in 1:5
    fig.add_subplot(510 + i, xlabel="2y/b", ylabel="circulation")
end
axs = fig.get_axes()

y2b = collect(range(0.0, stop=1, length=length(g5_cw)+1))

blade_xs = [cwblade1, cwblade2, cwblade3, cwblade4, cwblade5]
blades = [g1_cw, g2_cw, g3_cw, g4_cw, g5_cw]

gamma(x) = 27.78 * (x-0.06)^2 - 0.1

for i_blade in 1:5

    xs = blade_xs[i_blade]
    blade = blades[i_blade]

    dxs = [norm((xs[i,:] .+ xs[i+1,:])*0.5 - xs[1,:]) for i in 1:size(xs,1)-1]

    axs[i_blade-1].plot(dxs, blade)
    axs[i_blade-1].plot(dxs, gamma.(dxs))

    @show i_blade gamma.(dxs)

end

