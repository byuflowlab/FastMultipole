#####
##### test sensitivity of the threshold parameter theta
#####
scripts_dir = @__DIR__
include(joinpath(scripts_dir,"..","test","gravitational.jl"))

import Statistics as S

n_bodies = 30
xs = rand(n_bodies,3)
ms = rand(n_bodies)
masses = Vector{Mass}(undef,length(ms))
for i in 1:length(ms)
    local x = xs[i,:]
    local mass = [ms[i]]
    potential = zeros(1)
    force = zeros(3)
    masses[i] = Mass(x,mass,potential,force)
end

function test_accuracy(theta, expansion_order, n_per_branch=1)
    for i in 1:length(masses)
        masses[i].potential .*= 0
    end
    time_fmm = @elapsed tree = fmm.fmm!(masses, derivatives, expansion_order, n_per_branch, theta)

    potential_fmm = [mass.potential[1] for mass in masses]
    for i in 1:length(masses)
        masses[i].potential .*= 0
    end
    time_direct = @elapsed fmm.direct!(masses)
    potential_direct = [mass.potential[1] for mass in masses]

    err = potential_fmm - potential_direct
    rel_err = err ./ potential_direct
    return err, rel_err, time_fmm, time_direct
end

function make_plots(expansion_order)
    thetas = [4,8,16,32,64]
    n = length(thetas)
    errs = Vector{Vector{Float64}}(undef,n)
    rel_errs = Vector{Vector{Float64}}(undef,n)
    times_fmm = Vector{Float64}(undef,n)
    times_direct = Vector{Float64}(undef,n)
    for (i,theta) in enumerate(thetas)
        println("THETA = $theta")
        local err
        err, rel_err, time_fmm, time_direct = test_accuracy(theta, expansion_order)
        errs[i] = err
        rel_errs[i] = rel_err
        times_fmm[i] = time_fmm
        times_direct[i] = time_direct
    end

    fig_theta = plt.figure("theta")
    fig_theta.clear()
    ax = fig_theta.add_subplot(111,xlabel="mass index", ylabel="relative error")
    for i in 1:length(thetas)
        ax.plot(1:length(masses), abs.(rel_errs[i]), label=L"\theta = "*string(thetas[i]))
    end
    ax.legend()
    ax.set_yscale("log")

    fig_theta_2 = plt.figure("theta2")
    fig_theta_2.clear()
    ax = fig_theta_2.add_subplot(111, xlabel="theta", ylabel="RMS error")
    rms_err = [sqrt(S.mean(rel_errs[i].^2)) for i in 1:n]
    ax.plot(thetas, rms_err)
    ax.set_yscale("log")
end

# make_plots(1)
# make_plots(2)
# make_plots(3)
make_plots(4)
