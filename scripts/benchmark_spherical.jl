import Statistics as S
import FLOWFMM as fmm
import PyPlot as plt
using LaTeXStrings
import JLD

include(joinpath("..", "test", "gravitational.jl"))

version = "221012"

"""
Builds mass list, performs naive n^2 computation and stores to JLD file. Returns JLD file names.
"""
function benchmark_direct(ns, is_direct, base_name = version)
    file_names = String[]
    for (i,n) in enumerate(ns)
        jld_name = base_name*"_direct_n$n.jld"
        push!(file_names, jld_name)
        println("n = $n elements...")
        if !isfile(jld_name)
            ms = rand(n)
            xs = rand(n,3)
            masses = [Mass(xs[i,:],[ms[i]],zeros(1),zeros(3)) for i in 1:length(ms)]
            potentials_direct = zeros(length(masses))

            if i in is_direct
                println("\tComputing direct...")
                time_direct = @elapsed fmm.direct!(masses; reflex=false)
                for ii in 1:n
                    potentials_direct[ii] = masses[ii].potential[1]
                end
                println("\t\tdirect time: $(time_direct) seconds")
            else
                println("\tSkipping direct computation...")
                time_direct = NaN
                potentials_direct = [NaN for _ in 1:length(masses)]
            end
            println("\tSaving JLD...")
            JLD.save(jld_name, "base_name", base_name, "masses", masses, "potentials", potentials_direct, "time", time_direct, "n", n)
        else
            println("\tFile "*jld_name*" already exists; skipping")
        end
    end
    return file_names
end

"Loads mass list and direct potential from jld file, builds trees, and computes errors. Saves to new JLD files and returns file names."
function benchmark_fmm(direct_files, ns, expansion_orders, ns_per_branch, thetas, bases, base_name=version)
    # initialize file names
    files = String[]

    # assertions
    @assert length(direct_files) == length(expansion_orders) "number of files = $(length(direct_files)) but number of expansion orders = $(length(expansion_orders))"
    @assert length(direct_files) == length(ns_per_branch) "number of files = $(length(direct_files)) but number of ns_per_branch = $(length(ns_per_branch))"
    @assert length(direct_files) == length(thetas) "number of files = $(length(direct_files)) but number of thetas = $(length(thetas))"

    # dummy run to reduce influence of compile time on @elapsed measurement
    i, n = 1, 5
    ms = rand(n)
    xs = rand(n,3)
    masses = [Mass(xs[i,:],[ms[i]],zeros(1),zeros(3)) for i in 1:length(ms)]
    basis = fmm.Cartesian()
    @elapsed fmm.fmm!(masses, expansion_order, n_per_branch, theta, basis)

    println("\nBegin FMM Test:")
    for (i,direct_file) in enumerate(direct_files)
        # extract run info
        expansion_order = expansion_orders[i] # these are floats somehow
        n_per_branch = ns_per_branch[i] # these are floats somehow
        theta = thetas[i]
        basis = bases[i]

        n = ns[i]

        # form names
        if basis == fmm.Cartesian()
            basis_tag = "cartesian"
        elseif basis == fmm.Spherical()
            basis_tag = "spherical"
        else
            @error "basis $basis not supported"
        end
        jld_name = base_name*"_"*basis_tag*"_t$(theta)_nmax$(n_per_branch)_p$(expansion_order)_n$n.jld"

        if !isfile(jld_name)
            # extract jld info
            direct_data = JLD.load(direct_file)
            masses = direct_data["masses"]
            for mass in masses; mass.potential .*= 0; end # reset potential
            # n = direct_data["n"]
            # base_name = direct_data["base_name"]
            potentials_direct = direct_data["potentials"]

            potentials_fmm = similar(potentials_direct)

            println("n = $n elements...")
            println("\tBuilding $basis_tag Tree...")
            time_tree = @elapsed tree = fmm.Tree(masses, basis; expansion_order, n_per_branch)
            println("\t\tcartesian tree time: $(time_tree) seconds")

            println("\tComputing $basis_tag FMM...")
            time_fmm = @elapsed fmm.fmm!(tree, masses, theta, basis)
            println("\t\t$basis_tag FMM time: $(time_fmm) seconds")
            for ii in 1:n
                potentials_fmm[ii] = masses[ii].potential[1]
            end

            # mean error
            mean_err = abs(S.mean(abs.(potentials_direct .- potentials_fmm)) / S.mean(potentials_direct))

            println("\tSaving JLD as "*jld_name*"...")
            time_jld = @elapsed JLD.save(jld_name, "n", n, "base_name", base_name, "expansion_order", expansion_order,
                "n_per_branch", n_per_branch, "theta", theta, "basis", basis, "potentials", potentials_fmm,
                "time_tree", time_tree, "time", time_fmm, "mean_err", mean_err)
            println("\t\tJLD save time: $time_jld seconds")
        else
            println("\tFile "*jld_name*" already exists; skipping")
        end
        push!(files, jld_name)
    end
    println("Done.\n")

    return files
end

function sweep_n(ns, is_direct, expansion_order, n_per_branch, theta, base_name=version)
    # compute direct
    direct_files = benchmark_direct(ns, is_direct)
    expansion_orders = ones(Int64,length(ns)) .* expansion_order
    ns_per_branch = ones(Int64,length(ns)) .* n_per_branch
    thetas = ones(length(ns)) .* theta
    cartesian_bases = [fmm.Cartesian() for _ in 1:length(ns)]
    spherical_bases = [fmm.Spherical() for _ in 1:length(ns)]

    # compute cartesian
    cartesian_files = benchmark_fmm(direct_files, ns, expansion_orders, ns_per_branch, thetas, cartesian_bases, base_name)
    # compute spherical
    spherical_files = benchmark_fmm(direct_files, ns, expansion_orders, ns_per_branch, thetas, spherical_bases, base_name)

    return direct_files, cartesian_files, spherical_files
end

function load_jlds(files)
    ns = zeros(Int64, length(files))
    base_names = Vector{String}(undef,length(files))
    potentials_list = Vector{Vector{Float64}}(undef,length(files))
    times = Vector{Float64}(undef,length(files))
    times_tree = Vector{Float64}(undef, length(files))
    expansion_orders = Vector{Float64}(undef,length(files))
    ns_per_branch = Vector{Int64}(undef,length(files))
    thetas = Vector{Float64}(undef,length(files))
    mean_errs = Vector{Float64}(undef,length(files))
    for (i,file) in enumerate(files)
        data = JLD.load(file)
        ns[i] = data["n"]
        base_names[i] = data["base_name"]
        potentials_list[i] = data["potentials"]
        times[i] = data["time"]
        times_tree[i] = haskey(data, "time_tree") ? data["time_tree"] : NaN
        expansion_orders[i] = haskey(data,"expansion_order") ? data["expansion_order"] : -1
        ns_per_branch[i] = haskey(data, "n_per_branch") ? data["n_per_branch"] : -1
        thetas[i] = haskey(data, "theta") ? data["theta"] : NaN
        mean_errs[i] = haskey(data, "mean_err") ? data["mean_err"] : NaN
    end
    return ns, base_names, potentials_list, times, times_tree, expansion_orders, ns_per_branch, thetas, mean_errs
end

function load_jlds(direct_files, cartesian_files, spherical_files)
    # direct files
    ns_direct, base_names_direct, _, times_direct, _,
        _, _, _, _ = load_jlds(direct_files)

    # cartesian files
    ns_cartesian, base_names_cartesian, _, times_cartesian, times_tree_cartesian,
        expansion_orders_cartesian, ns_per_branch_cartesian, thetas_cartesian, mean_errs_cartesian = load_jlds(cartesian_files)

    # spherical files
    ns_spherical, base_names_spherical, _, times_spherical, times_tree_spherical,
        expansion_orders_spherical, ns_per_branch_spherical, thetas_spherical, mean_errs_spherical = load_jlds(spherical_files)

    return ns_direct, times_direct, ns_cartesian, times_cartesian, ns_spherical, times_spherical, mean_errs_cartesian, mean_errs_spherical
end

function load_jld(file)
    data = JLD.load(file)
    theta = data["theta"]
    n_per_branch = data["n_per_branch"]
    expansion_order = data["expansion_order"]
    ns = data["ns"]
    is_direct = data["is_direct"]
    times_tree_spherical = data["times_tree_spherical"]
    times_fmm_spherical = data["times_fmm_spherical"]
    times_tree_cartesian = data["times_tree_cartesian"]
    times_fmm_cartesian = data["times_fmm_cartesian"]
    times_direct = data["times_direct"]
    mean_errs_spherical = data["mean_errs_spherical"]
    mean_errs_cartesian = data["mean_errs_cartesian"]
    return theta, n_per_branch, expansion_order, ns, is_direct, times_tree_spherical, times_fmm_spherical, times_tree_cartesian, times_fmm_cartesian, times_direct, mean_errs_spherical, mean_errs_cartesian
end

function loglog_slope(xs, ys)
    log_xs = log.(xs)
    log_ys = log.(ys)
    # minimize Ax=b where x is the line, A is [x 1], and b is ys
    A = hcat(log_xs, ones(length(xs)))
    X = A \ log_ys
    return X[1]
end

function plot_files(direct_files, cartesian_files, spherical_files, save_name, fig_title)
    ns_direct, times_direct, ns_cartesian, times_cartesian, ns_spherical,
        times_spherical, mean_errs_cartesian, mean_errs_spherical = load_jlds(direct_files, cartesian_files, spherical_files)

    fig = plt.figure("benchmark_fmm")
    fig.clear()
    fig.add_subplot(111,xlabel="elements",ylabel="time, seconds")
    # fig.suptitle(L"\theta="*"$theta, "*L"N_{max}="*"$(n_per_branch), "*L"P="*"$(expansion_order)")
    fig.suptitle(fig_title)
    ax = fig.get_axes()[1]
    ax.plot(ns_cartesian, times_cartesian, label="cartesian fmm: slope=$(loglog_slope(ns_cartesian,times_cartesian))")
    ax.plot(ns_spherical, times_spherical, label="spherical fmm: slope=$(loglog_slope(ns_spherical,times_spherical))")
    ax.plot(ns_direct, times_direct, label="direct: slope=$(loglog_slope(ns_direct[is_direct],times_direct[is_direct]))")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    fig.savefig(save_name*".png")

    fig = plt.figure("error")
    fig.clear()
    fig.add_subplot(111,xlabel="elements", ylabel="mean relative error")
    fig.suptitle(fig_title)
    ax = fig.get_axes()[1]
    ax.plot(ns_cartesian, mean_errs_cartesian, label="cartesian")
    ax.plot(ns_spherical, mean_errs_spherical, label="spherical")
    ax.legend()
    fig.savefig(save_name*"_error"*".png")

    return nothing
end

function benchmark_n(ns, is_direct, expansion_order, n_per_branch, theta, base_name=version)
    direct_files, cartesian_files, spherical_files = sweep_n(ns, is_direct, expansion_order, n_per_branch, theta, base_name)

    save_name = base_name*"_t$(theta)_nmax$(n_per_branch)_p$(expansion_order)"
    fig_title = L"\theta="*"$theta, "*L"N_{max}="*"$(n_per_branch), "*L"P="*"$(expansion_order)"

    plot_files(direct_files, cartesian_files, spherical_files, save_name, fig_title)
end

#####
##### run new
#####
# set up
ns = [10, 100, 1000, 10000]#, 100000, 1000000]# , 5000000]
is_direct = 1:4
theta = 4
n_per_branch = 50
expansion_order = 2
benchmark_n(ns, is_direct, expansion_order, n_per_branch, theta)
