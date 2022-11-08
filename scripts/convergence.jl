scripts_dir = @__DIR__
include(joinpath(scripts_dir,"..","test","gravitational.jl"))

# explore error with expansion order
xs = [
    1.2 1.1 0.8;
    0.8 0.9 0.2;
    0.1 0.2 0.9;
    0.1 0.3 0.2;
    0.2 0.25 0.4
]

ms = [
    0.8,
    1.1,
    2.2,
    0.5,
    1.9
]
function test_accuracy(exp_order)

    masses = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses[i] = Mass(x,mass,potential,force)
    end
    tree = fmm.Tree(masses, B2M!, P2P!; expansion_order=exp_order)

    masses_2 = Vector{Mass}(undef,length(ms))
    for i in 1:length(ms)
        x = xs[i,:]
        mass = [ms[i]]
        potential = zeros(1)
        force = zeros(3)
        masses_2[i] = Mass(x,mass,potential,force)
    end

    fmm.direct!(masses_2)
    potential_direct = [mass.potential[1] for mass in masses_2]

    theta = 4
    fmm.fmm!(tree, masses, theta; reset_tree=true)
    potential_fmm = [mass.potential[1] for mass in masses]

    return potential_direct, potential_fmm
end

function test_err(exp_order)
    potential_direct, potential_fmm = test_accuracy(exp_order)
    err = abs.(potential_direct .- potential_fmm)
    rel_err = err ./ abs.(potential_direct)
    return err, rel_err
end

os = [1,2,3,4]#,5,6,7]
errs = zeros(length(ms),length(os))
rel_errs = zeros(length(ms),length(os))
for (i,o) in enumerate(os)
    global this_rel_err, this_err
    this_rel_err, this_err = test_err(o)
    errs[:,i] .= this_err
    rel_errs[:,i] .= this_rel_err
end

conv_1 = plt.figure("convergence_1")
conv_1.clear()
conv_1.add_subplot(111,xlabel="order",ylabel="relative error")
ax = conv_1.get_axes()[1]
for i in 1:length(ms)
    ax.plot(os, rel_errs[i,:], label="mass $i")
end
# ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
conv_1.savefig("expansion_convergence.png")
