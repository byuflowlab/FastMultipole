#=

# relative costs of various FMM tasks required for load balancing

# defaults
const TAU_B2M_DEFAULT = SVector{3}(1.4977347448647953e-7, 2.2143107802882807e-10, 3.1288829471775396e-9)
const ALLOC_M2M_L2L_DEFAULT = SVector{3}(6.607141609261922e-7, 9.642838049805685e-9, 3.9650702246920374e-10)
const TAU_M2M_DEFAULT = SVector{5}(2.5342296810051063e-8, 1.1457176734203204e-7, 1.7874812696411806e-9, 4.076474127582875e-9, 5.029318715939952e-10)
const ALLOC_M2L_DEFAULT = SVector{3}(1.233435958541507e-7, 1.8546354249173343e-8, 1.4835385731876622e-9)
const TAU_M2L_DEFAULT = SVector{5}(1.1527524574372772e-7, 7.267606919879385e-8, 5.275384866214009e-8, 5.435225827962213e-9, 2.613581164941739e-9)
const ALLOC_L2B_DEFAULT = SVector{3}(2.2602257309941528e-6, 1.360788334472431e-9, 1.6180982759930193e-9)
const TAU_L2L_DEFAULT = SVector{5}(2.5342296810051063e-8, 1.1457176734203204e-7, 1.7874812696411806e-9, 4.076474127582875e-9, 5.029318715939952e-10)
const TAU_L2B_DEFAULT = SVector{3}(6.251536970899492e-7, 7.584330604198042e-8, 3.705055312794977e-8)
const C_NEARFIELD_DEFAULT = 2.5e-9
const COST_PARAMETERS = ["ALLOC_M2M_L2L","TAU_M2M_L2L","ALLOC_M2L","TAU_M2L","ALLOC_L2B","TAU_L2B","ERROR","C_NEARFIELD","ERROR_NEARFIELD","TAU_B2M_1","TAU_B2M_2","TAU_B2M_3","ERROR_B2M"]

struct BenchmarkSystem{TF}
    val::MVector{3,TF}
    mat::MMatrix{3,3,TF,9}
end

BenchmarkSystem(TF) = BenchmarkSystem{TF}(rand(MVector{3,TF}),rand(MMatrix{3,3,TF,9}))

Base.getindex(::BenchmarkSystem,i_body,::Position) = rand(SVector{3,Float64})
Base.getindex(::BenchmarkSystem,i_body,::ScalarStrength) = rand()
Base.setindex(sys::BenchmarkSystem,val,i_body,::ScalarPotential) = sys.val[1] += val
Base.setindex(sys::BenchmarkSystem,val,i_body,::VectorPotential) = sys.val .+= val
Base.setindex(sys::BenchmarkSystem,val,i_body,::Velocity) = sys.val .+= val
Base.setindex(sys::BenchmarkSystem,val,i_body,::VelocityGradient) = sys.mat .+= val

function allocate_m2m_l2l(expansion_order, type=Float64)
    zeros(Complex{type}, (expansion_order+1)*(expansion_order+1)), zeros(type,4)
end

function allocate_m2l(expansion_order, type=Float64)
    zeros(Complex{type}, (expansion_order<<1 + 1)*(expansion_order<<1 + 1))
end

function allocate_l2b(expansion_order, type=Float64)
    vector_potential = zeros(type,3)
    potential_jacobian = zeros(type,3,4)
    potential_hessian = zeros(type,3,3,4)
    derivative_harmonics = zeros(Complex{type}, ((expansion_order+1) * (expansion_order+2)) >> 1)
    derivative_harmonics_theta = zeros(Complex{type}, ((expansion_order+1) * (expansion_order+2)) >> 1)
    derivative_harmonics_theta_2 = zeros(Complex{type}, ((expansion_order+1) * (expansion_order+2)) >> 1)
    workspace = zeros(type,3,4)
    return vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace
end

function get_tau_b2m(system, branch, b2m_n_bodies, harmonics_m2m_l2l, expansion_order)
    bodies_index = 1:min(b2m_n_bodies,length(system))
    tau_b2m = @elapsed B2M!(system, branch, bodies_index, harmonics_m2m_l2l, expansion_order)
    tau_b2m = @elapsed B2M!(system, branch, bodies_index, harmonics_m2m_l2l, expansion_order)
    return tau_b2m / length(bodies_index)
end

function get_tau_b2m(systems::Tuple, branch, b2m_n_bodies, harmonics_m2m_l2l, expansion_order)
    return SVector{length(systems),Float64}(get_tau_b2m(system, branch, b2m_n_bodies, harmonics_m2m_l2l, expansion_order) for system in systems)
end

function estimate_tau_fmm(expansion_order, type=Float64)
    branch1 = Branch(1:5, 1, 1:1, rand(SVector{3,type}), rand(), expansion_order)
    branch2 = Branch(1:5, 1, 1:1, rand(SVector{3,type}), rand(), expansion_order)
    harmonics_m2m_l2l, ML = allocate_m2m_l2l(expansion_order, type)
    alloc_m2m_l2l = @elapsed harmonics_m2m, M = allocate_m2m_l2l(expansion_order, type)
    alloc_m2m_l2l = @elapsed harmonics_m2m, M = allocate_m2m_l2l(expansion_order, type)
    tau_m2m_l2l = @elapsed M2M!(branch1, branch2, harmonics_m2m_l2l, ML, expansion_order)
    tau_m2m_l2l = @elapsed M2M!(branch1, branch2, harmonics_m2m_l2l, ML, expansion_order)
    harmonics_m2l = allocate_m2l(expansion_order, type)
    alloc_m2l = @elapsed harmonics_m2l = allocate_m2l(expansion_order, type)
    alloc_m2l = @elapsed harmonics_m2l = allocate_m2l(expansion_order, type)
    tau_m2l = @elapsed M2L!(branch1, branch2, harmonics_m2l, expansion_order)
    tau_m2l = @elapsed M2L!(branch1, branch2, harmonics_m2l, expansion_order)
    # harmonics_l2l, L = allocate_m2m_l2l(expansion_order, type)
    # alloc_l2l = @belapsed harmonics_l2l, L = allocate_m2m_l2l(expansion_order, type)
    # tau_l2l = @belapsed L2L!(branch1, branch2, harmonics_l2l, L, expansion_order)
    vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = allocate_l2b(expansion_order, type)
    alloc_l2b = @elapsed vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = allocate_l2b(expansion_order, type)
    alloc_l2b = @elapsed vector_potential, potential_jacobian, potential_hessian, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, workspace = allocate_l2b(expansion_order, type)
    body_position = rand(SVector{3,type})
    expansion_center = 10 * rand(SVector{3,type})
    tau_l2b = @elapsed L2B_loop!(vector_potential, potential_jacobian, potential_hessian, body_position, expansion_center, 
                branch2.local_expansion, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, expansion_order, workspace)
    tau_l2b = @elapsed L2B_loop!(vector_potential, potential_jacobian, potential_hessian, body_position, expansion_center, 
                branch2.local_expansion, derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, expansion_order, workspace)
    return alloc_m2m_l2l, tau_m2m_l2l, alloc_m2l, tau_m2l, alloc_l2b, tau_l2b
end

function get_error(param, mat, list)
    return maximum(abs.(mat * param - list))
end

function get_system_parameters(systems::Tuple, expansion_orders; n_bodies_max=1000, n_bodies_min=10, n_benchmarks=10)
    params = MVector{length(systems),Float64}(0.0 for _ in systems)
    errors = MVector{length(systems),Float64}(0.0 for _ in systems)
    b2m_params = SVector{length(systems),MVector{3,Float64}}(MVector{3}(zeros(3)) for _ in systems)
    b2m_errors = MVector{length(systems),Float64}(0.0 for _ in systems)
    for (i,system) in enumerate(systems)
        C, error, b2m_param, b2m_error = get_system_parameters(system, expansion_orders; n_bodies_max=n_bodies_max, n_bodies_min=n_bodies_min, n_benchmarks=n_benchmarks)
        params[i] = C
        errors[i] = error
        b2m_params[i] = b2m_param
        b2m_errors[i] = b2m_error
    end
    return SVector{length(systems)}(params), errors, b2m_params, b2m_errors
end

function get_system_parameters(system, expansion_orders; n_bodies_max=5000, n_bodies_min=100, n_benchmarks=5, epsilon=0.1)
    # set up simulation
    n_bodies_max = min(length(system),n_bodies_max)
    n_bodies_list = range(n_bodies_min, step=div(n_bodies_max-n_bodies_min,n_benchmarks), length=n_benchmarks)
    nearfield_benchmarks = MVector{length(n_bodies_list),Float64}(0.0 for _ in 1:length(n_bodies_list))
    b2m_benchmarks = MMatrix{length(expansion_orders),length(n_bodies_list),Float64,length(expansion_orders)*length(n_bodies_list)}(0.0 for _ in 1:length(expansion_orders)*length(n_bodies_list))
    branches = Tuple(Branch(1:1, 0, 1:0, rand(SVector{3,Float64}), rand(), expansion_order) for expansion_order in expansion_orders)
    harmonics_list = Tuple(zeros(Complex{eltype(system)}, (expansion_order+1)*(expansion_order+1)) for expansion_order in expansion_orders)

    # store original states
    original_states = [system[1]]
    resize!(original_states, n_bodies_list[end])
    for i in 2:n_bodies_list[end]
        original_states[i] = system[i]
    end

    # run benchmarks
    println("===== Running Nearfield Benchmarks =====")
    for (i,n_bodies) in enumerate(n_bodies_list)
        println("\tn bodies: n_bodies")
        bodies_index = 1:n_bodies
        nearfield_benchmarks[i] = @elapsed direct!(system, bodies_index, system, bodies_index)
        nearfield_benchmarks[i] = @elapsed direct!(system, bodies_index, system, bodies_index)
        for (j,(expansion_order,branch,harmonics)) in enumerate(zip(expansion_orders,branches,harmonics_list))
            b2m_benchmarks[j,i] = @elapsed B2M!(system, branch, bodies_index, harmonics, expansion_order)
            b2m_benchmarks[j,i] = @elapsed B2M!(system, branch, bodies_index, harmonics, expansion_order)
            b2m_benchmarks[j,i] /= n_bodies
        end
    end
    println("===== Done. =====")

    # restore original states
    for i in 1:n_bodies_list[end]
        system[i] = original_states[i]
    end

    # estimate nearfield parameters
    Cs = [bm/n_bodies^2 for (bm,n_bodies) in zip(nearfield_benchmarks, n_bodies_list)]
    C = sum(Cs)/length(Cs)
    this_error = 0.0
    for this_C in Cs
        this_error = max(abs(this_C-C),this_error)
    end
    this_error /= C # relative error
    this_error > epsilon && (@warn "relative error in nearfield cost estimate greater than epsilon; load balancing may have unexpected behavior; Cs = Cs")

    # estimate b2m parameters
    b2m_mean = MVector(sum(b2m_benchmarks,dims=2)) ./ size(b2m_benchmarks,2)
    b2m_error = maximum(abs.(b2m_benchmarks .- b2m_mean) ./ b2m_mean)
    b2m_error > epsilon && (@warn "relative error in b2m cost estimate greater than epsilon; load balancing may have unexpected behavior; b2m = (b2m_benchmarks)")
    mat_nx3 = MMatrix{length(expansion_orders),3,Float64}(p^i for p in expansion_orders, i in 0:2)
    b2m_params = mat_nx3 \ b2m_mean

    return C, this_error, SVector(b2m_params), b2m_error
end

function read_param(row,i)
    @assert row[1] == COST_PARAMETERS[i] "cost parameters file formatted incorrectly"
    i_end = findfirst((val) -> !(typeof(val) <: Number), view(row,2:length(row)))
    isnothing(i_end) && (i_end = length(row))
    if i_end > 2
        return SVector{i_end-1,Float64}(view(row,2:i_end))
    else # single value
        return row[2]
    end
end

function read_param_b2m(b2m_c1::AbstractVector{<:AbstractVector}, b2m_c2, b2m_c3)
    b2m_params = SVector{length(b2m_c1)}(SVector{3,Float64}(c1,c2,c3) for (c1,c2,c3) in zip(b2m_c1, b2m_c2, b2m_c3))
end

function read_param_b2m(b2m_c1, b2m_c2, b2m_c3)
    b2m_params = SVector{3,Float64}(b2m_c1,b2m_c2,b2m_c3)
end

function read_cost_parameters(cost_file_path)
    file_contents = readdlm(cost_file_path,',')
    params = Tuple(read_param(view(file_contents,i,:),i) for i in 1:6)
    errors = read_param(view(file_contents,7,:),7)
    nearfield_params = read_param(view(file_contents,8,:),8)
    nearfield_error = read_param(view(file_contents,9,:),9)
    b2m_c1 = read_param(view(file_contents,10,:),10)
    b2m_c2 = read_param(view(file_contents,11,:),11)
    b2m_c3 = read_param(view(file_contents,12,:),12)
    b2m_params = read_param_b2m(b2m_c1, b2m_c2, b2m_c3)
    
    b2m_error = read_param(view(file_contents,13,:),13)
    return params, errors, nearfield_params, nearfield_error, b2m_params, b2m_error
end

function unpack_b2m!(file_contents, b2m_params::AbstractVector{<:AbstractVector})
    for (i,b2m_param) in enumerate(b2m_params)
        file_contents[1:3,i+1] .= b2m_param
    end
end

function unpack_b2m!(file_contents, b2m_params)
    file_contents[1:3,2] .= b2m_params
end

function write_cost_parameters(cost_file_path, params, errors, nearfield_params, nearfield_error, b2m_params, b2m_error)
    @assert length(params) + 7 == length(COST_PARAMETERS) "received the wrong number of cost parameters for file write; expected (length(COST_PARAMETERS)); got (length(params) + 7)"
    max_length = 0
    for param in (params..., errors, nearfield_params, nearfield_error, b2m_params, b2m_error)
        max_length = max(max_length,length(param))
    end
    file_contents = Matrix{Any}(undef,length(COST_PARAMETERS),1 + max_length)
    file_contents .= ""
    file_contents[:,1] .= COST_PARAMETERS
    for (i,(param_name,param)) in enumerate(zip(COST_PARAMETERS[1:9], (params..., errors, nearfield_params, nearfield_error)))
        file_contents[i,2:length(param)+1] .= param
    end
    file_contents[end,2] = b2m_error
    unpack_b2m!(view(file_contents,10:12,:), b2m_params)

    open(cost_file_path, "w") do io
        writedlm(io, file_contents, ',')
    end
    return nothing
end

initialize_tau_b2m(systems::Tuple) = SVector{length(systems),MVector{3,Float64}}(MVector{3,Float64}(zeros(3)) for _ in 1:length(systems))

initialize_tau_b2m(systems) = MVector{3,Float64}(zeros(3))

function initialize_cost_parameters(systems, param_name)
    if param_name == "TAU_B2M"
        return initialize_tau_b2m(systems)
    elseif param_name in ("TAU_M2M_L2L", "TAU_M2L")
        return MVector{5,Float64}(zeros(5))
    else
        return MVector{3,Float64}(zeros(3))
    end
end

function estimate_tau(systems, type=Float64; expansion_orders = 1:3:20, epsilon=0.1, cost_file_read=true, cost_file_write=true, cost_file_path="cost_parameters_(typeof(systems)).csv")
    read_file = false
    params_index = 1:6
    if cost_file_read && isfile(cost_file_path)
        read_file = true
        println("Reading cost parameter file (cost_file_path)...")
        params, errors, nearfield_params, nearfield_error = read_cost_parameters(cost_file_path)
    else
        # preallocate
        lists = [Vector{Float64}(undef, length(expansion_orders)) for _ in params_index]
        params = Tuple(
            initialize_cost_parameters(systems, param_name) for param_name in COST_PARAMETERS[params_index]
        )
        errors = zeros(length(params))
        mat_nx3 = Float64[p^i for p in expansion_orders, i in 0:2]
        mat_nx5 = Float64[p^i for p in expansion_orders, i in 0:4]

        # get fmm benchmarks
        println("===== BEGIN FMM BENCHMARKS =====")
        for (ip,expansion_order) in enumerate(expansion_orders)
            println("\tp = expansion_order")
            benchmarks = estimate_tau_fmm(expansion_order, type)
            for (list, benchmark) in zip(lists, benchmarks)
                list[ip] = benchmark
            end
        end
        println("===== Done. =====")
        
        # estimate fmm parameters and errors
        for (i,(list,param)) in enumerate(zip(lists,params))
            if length(param) == 3
                param .= mat_nx3 \ list
                errors[i] = get_error(param, mat_nx3, list)
            elseif length(param) == 5
                param .= mat_nx5 \ list
                errors[i] = get_error(param, mat_nx5, list)
            else
                error("load balance parameter has the incorrect size; expected 3 or 5, got length(param)")
            end
        end
        
        for (i,list) in enumerate(lists)
            errors[i] /= list[end]
        end
        
        # get nearfield benchmarks
        nearfield_params, nearfield_error, b2m_params, b2m_error = get_system_parameters(systems, expansion_orders)
    end

    maximum(errors) > epsilon && (@warn "error in load-balancing- inefficient multithreading behavior may occur; relative error=(maximum(errors))")

    if cost_file_write && !read_file
        println("Writing cost parameter file (cost_file_path)...")
        write_cost_parameters(cost_file_path, params, errors, nearfield_params, nearfield_error, b2m_params, b2m_error)
    end

    # return result
    return params, errors, nearfield_params, nearfield_error, b2m_params, b2m_error
end

function get_t(n_times, coefficients, expansion_order)
    t = get_t(coefficients, expansion_order)
    return n_times * t
end

function get_t(n_times, coefficients::AbstractVector{<:AbstractVector}, expansion_order)
    t = get_t(coefficients, expansion_order)
    return sum(n_times .* t)
end

function get_t(coefficients, expansion_order)
    t = 0.0
    for (i,coefficient) in enumerate(coefficients)
        t += coefficient * expansion_order^(i-1)
    end
    return t
end

function get_t(coefficients::AbstractVector{<:AbstractVector}, expansion_order)
    t = SVector{length(coefficients)}(zeros(length(coefficients)))
    for (i,coefficient) in enumerate(coefficients)
        t = t + coefficient * expansion_order^(i-1)
    end
    return t
end

=#
get_n_bodies(branch::SingleBranch) = length(branch.bodies_index)
get_n_bodies(branch::MultiBranch) = SVector{length(branch.bodies_index)}(length(index) for index in branch.bodies_index)
#=
function get_t_nearfield(direct_list, branches::Vector{<:SingleBranch}, C_nearfield)
    t = 0.0
    for (i_target, j_source) in direct_list
        n_bodies_source = length(branches[j_source].bodies_index)
        n_bodies_target = length(branches[i_target].bodies_index)
        t += C_nearfield * n_bodies_source * n_bodies_target
    end
    return t
end

function get_t_nearfield(direct_list, branches::Vector{<:MultiBranch}, C_nearfield)
    t = 0.0
    n_systems = length(C_nearfield)
    for (i_target, j_source) in direct_list
        n_bodies_target = 0
        for i_system in 1:n_systems
            n_bodies_target += length(branches[i_target].bodies_index[i_system])
        end
        for j_source_system in 1:n_systems
            t += C_nearfield[i_system] * length(branches[j_source].bodies_index[j_source_system]) * n_bodies_target
        end
    end
    return t
end

function estimate_cost(direct_list, m2l_list, branches, expansion_order, cost_parameters)
    # total number of bodies (svector if a multisystem)
    n_bodies_total = get_n_bodies(branches[1])
        
    # cost of b2m for all bodies across all systems
    t_B2M = get_t(n_bodies_total, cost_parameters.tau_B2M, expansion_order)

    # cost of m2m and l2l for all branches, including allocations per thread TODO: perhaps neglect this?
    t_M2M_L2L = 2 * get_t(length(branches)-1, cost_parameters.tau_M2M_L2L, expansion_order)
    a_M2M_L2L = 2 * get_t(cost_parameters.alloc_M2M_L2L, expansion_order)

    # cost of m2l for all branches, including allocations per thread
    t_M2L = get_t(length(m2l_list), cost_parameters.tau_M2L, expansion_order)
    a_M2L = get_t(cost_parameters.alloc_M2L, expansion_order)

    # cost of l2b for all bodies, including allocations per thread
    t_L2B = get_t(n_bodies_total, cost_parameters.tau_L2B, expansion_order)
    a_L2B = get_t(cost_parameters.alloc_L2B, expansion_order)

    # cost of nearfield interactions for all bodies across all systems
    t_nearfield = get_t_nearfield(direct_list, branches, cost_parameters.C_nearfield)

    return t_B2M, t_M2M_L2L, a_M2M_L2L, t_M2L, a_M2L, t_L2B, a_L2B, t_nearfield
end

function load_balance_residual(n_threads_fmm, direct_list, m2l_list, branches, expansion_order, cost_parameters, n_threads=Threads.nthreads(); tau_spawn = 0.0)
    # total number of bodies (svector if a multisystem)
    n_bodies_total = get_n_bodies(branches[1])
    
    # cost of b2m for all bodies across all systems
    t_b2m = get_t(n_bodies_total, cost_parameters.tau_b2m, expansion_order)
    
    # cost of m2m and l2l for all branches, including allocations per thread TODO: perhaps neglect this?
    t_m2m_l2l = 2 * get_t(length(branches)-1, cost_parameters.tau_m2m_l2l, expansion_order)
    a_m2m_l2l = 2 * get_t(cost_parameters.alloc_m2m_l2l, expansion_order)
    
    # cost of m2l for all branches, including allocations per thread
    t_m2l = get_t(length(m2l_list), cost_parameters.tau_m2l, expansion_order)
    a_m2l = get_t(cost_parameters.alloc_m2l, expansion_order)

    # cost of l2b for all bodies, including allocations per thread
    t_l2b = get_t(length(), cost_parameters.tau_l2b, expansion_order)

    # cost of nearfield interactions for all bodies across all systems
    t_nearfield = get_t_nearfield(direct_list, branches, cost_parameters.C_nearfield)

    # residual function for determining optimal thread designation
    residual = 1/n_threads_fmm * t_b2m + tau_spawn * n_threads_fmm +
                t_m2l / n_threads_fmm + a_m2l * n_threads_fmm + tau_spawn * n_threads_fmm +
                t_l2b / n_threads_fmm + a_l2b * n_threads_fmm + tau_spawn * n_threads_fmm +
                2 * (t_m2m_l2l / n_threads_fmm + a_m2m_l2l * n_threads_fmm) -
                1/(n_threads - n_threads_fmm) * t_nearfield - tau_spawn * (n_threads - n_threads_fmm)
    
    return residual
end

# function load_balance_nearfield(direct_list, m2l_list, branches, cost_parameters)

# end
=#

function direct_cost_estimate(system, n_per_branch; n_iter=10)
    # store original states
    original_states = [system[1]]
    resize!(original_states, n_per_branch)
    for i in 2:n_per_branch
        original_states[i] = system[i]
    end

    # benchmark
    t = 0.0
    for i in 1:n_iter+1 # one for precompilation
        i > 1 && (t += @elapsed direct!(system, 1:n_per_branch, system, 1:n_per_branch))
    end
    t /= n_iter # mean time per iteration
    t /= n_per_branch^2 # mean time per interaction

    # restore original states
    for i in 1:n_per_branch
        system[i] = original_states[i]
    end

    return t
end

function dummy_direct_cost_estimate(system, n_per_branch)
    return NaN
end

function direct_cost_estimate(systems::Tuple, n_per_branch; n_iter=10)
    return SVector{length(systems),Float64}(direct_cost_estimate(system, n_per_branch; n_iter=n_iter) for system in systems)
end

function dummy_direct_cost_estimate(systems::Tuple, n_per_branch)
    return SVector{length(systems),Float64}(dummy_direct_cost_estimate(system, n_per_branch) for system in systems)
end

