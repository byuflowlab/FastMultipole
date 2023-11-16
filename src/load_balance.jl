# relative costs of various FMM tasks required for load balancing

const TAU_B2M_DEFAULT = 1.0
const ALLOC_M2M_DEFAULT = 1.0
const TAU_M2M_DEFAULT = 1.0
const ALLOC_M2L_DEFAULT = 1.0
const TAU_M2L_DEFAULT = 1.0
const TAU_L2L_DEFAULT = 1.0
const ALLOC_L2L_DEFAULT = 1.0
const TAU_L2B_DEFAULT = 1.0
const ALLOC_L2B_DEFAULT = 1.0

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

function estimate_tau(branch1, branch2, dummy_sys, expansion_order, type=Float64)
    alloc_m2m = @belapsed harmonics_m2m, M = allocate_m2m_l2l($expansion_order, $type)
    dummy_index = 1:1
    tau_b2m = @belapsed B2M!_sourcepoint($dummy_sys, $branch1, $dummy_index, $harmonics_m2m, $expansion_order)
    tau_m2m = @belapsed M2M!($branch1, $branch2, $harmonics_m2m, $M, $expansion_order)
    alloc_m2l = @belapsed harmonics_m2l = allocate_m2l($expansion_order, $type)
    tau_m2l = @belapsed M2L!($branch1, $branch2, $harmonics_m2l, $expansion_order)
    alloc_l2l = @belapsed harmonics_l2l, L = allocate_m2m_l2l($expansion_order, $type)
    tau_l2l = @belapsed L2L!($branch1, $branch2, $harmonics_l2l, $L, $expansion_order)
    alloc_l2b = @belapsed vector_potential, potential_jacobian, potential_hessian, 
                derivative_harmonics, derivative_harmonics_theta, derivative_harmonics_theta_2, 
                workspace = allocate_l2b($expansion_order, $type)
    body_position = rand(SVector{3,type})
    tau_l2b = @belapsed L2B_loop!($vector_potential, $potential_jacobian, $potential_hessian, $body_position, $expansion_center, 
                $local_expansion, $harmonics, $harmonics_theta, $harmonics_theta_2, $expansion_order, $workspace)
    return alloc_m2m, tau_b2m, tau_m2m, alloc_m2l, tau_m2l, tau_l2l, alloc_l2b, tau_l2b
end

function get_error(param, mat, list)
    return maximum(abs.(mat * param - list))
end

function estimate_tau(type=Float64; expansion_orders = 1:20)
    # preallocate
    branch1 = Branch(1:5, 1, 1:1, rand(SVector{3,type}), rand(), expansion_order)
    branch2 = Branch(1:5, 1, 1:1, rand(SVector{3,type}), rand(), expansion_order)
    dummy_sys = BenchmarkSystem(type)
    lists = [Vector{Float64}(undef, length(expansion_orders)) for _ in 1:9]
    params = Tuple(
        zeros(3), zeros(3), zeros(3), zeros(3), zeros(5), zeros(3), zeros(3), zeros(3)
    )
    errors = zeros(length(params))
    mat_nx3 = SMatrix{length(expansion_orders),3,Float64}(p^i for p in expansion_orders, i in 0:2)
    mat_nx5 = SMatrix{length(expansion_orders),3,Float64}(p^i for p in expansion_orders, i in 0:4)

    # get benchmarks
    for (ip,p) in enumerate(expansion_orders)
        benchmarks = estimate_tau(branch1, branch2, dummy_sy, p, type)
        for (list, benchmark) in zip(lists, benchmarks)
            list[ip] = benchmark
        end
    end

    # estimate parameters and errors
    for (i,(list,param)) in enumerate(zip(lists,params))
        if length(param) == 3
            param .= mat_nx3 \ list
            errors[i] = get_error(param, mat_nx3, list)
        else
            param .= mat_nx5 \ list
            errors[i] = get_error(param, mat_nx5, list)
        end
    end

    # return result
    return params
end
