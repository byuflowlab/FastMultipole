using Statistics

function spherical_to_cartesian(ρ,θ,ϕ)
    z = ρ * cos(θ)
    x = ρ * sin(θ) * cos(ϕ)
    y = ρ * sin(θ) * sin(ϕ)
    return x, y, z
end

function build_system(center, ρs, θs, ϕs, radii, qs)
    source_bodies = zeros(8,length(ρs))
    for (i,(ρ,θ,ϕ,r,q)) in enumerate(zip(ρs, θs, ϕs, radii, qs))
        x, y, z = spherical_to_cartesian(ρ,θ,ϕ)
        source_bodies[1:3,i] .= center + SVector{3}(x,y,z)
        source_bodies[4,i] = r
        source_bodies[5,i] = q
    end
    system = Gravitational(source_bodies)
    return system
end

function generate_multipole(center, ρs, θs, ϕs, radii, qs, expansion_order, shrink)
    # build system
    system = build_system(center, ρs, θs, ϕs, radii, qs)

    # build branch
    branch = Branch(1:length(ρs), 0, 1:0, 0, 1, center, sqrt(3.0), sqrt(3.0), SVector{6}(-1.0,1.0,-1.0,1.0,-1.0,1.0), SVector{3}(1.0,1.0,1.0), expansion_order)

    # shrink branch
    branches = [branch]
    shrink && FastMultipole.shrink_recenter!(branches, [1:1], system)
    branch = branches[1]

    # multipole coefficients
    body_to_multipole!(branch, system, branch.harmonics, Val(expansion_order))

    return branch, system
end

function evaluate_multipole(xt, branch::Branch, expansion_order)
    ϕ_m2b, v_m2b, g_m2b = evaluate_multipole(xt, branch.center, branch.multipole_expansion, DerivativesSwitch(), Val(expansion_order))
    return ϕ_m2b
end

function evaluate_direct(xt, system::Gravitational)
    u = 0.0
    for body in system.bodies
        xs = body.position
        q = body.strength
        u += q / norm(xt-xs)
    end
    return u / 4 / pi
end

function multipole_to_local!(local_branch, multipole_branch, expansion_order)
    # reset expansion
    local_branch.local_expansion .= zero(eltype(local_branch.local_expansion))

    # preallocate containers
    lamb_helmholtz = Val(false)
    Hs_π2 = [1.0]
    FastMultipole.update_Hs_π2!(Hs_π2, Val(expansion_order))
    Ts = zeros(FastMultipole.length_Ts(expansion_order))
    eimϕs = zeros(2, expansion_order+1)
    weights_tmp_1 = initialize_expansion(expansion_order, eltype(Ts))
    weights_tmp_2 = initialize_expansion(expansion_order, eltype(Ts))

    # normalization
    ζs_mag = zeros(FastMultipole.length_ζs(expansion_order))
    FastMultipole.update_ζs_mag!(ζs_mag, 0, expansion_order)
    ηs_mag = zeros(FastMultipole.length_ηs(expansion_order))
    FastMultipole.update_ηs_mag!(ηs_mag, 0, expansion_order)

    # local coefficients
    FastMultipole.multipole_to_local!(local_branch, multipole_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ζs_mag, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)

    return nothing
end

function evaluate_local(xt, branch::Branch, expansion_order)
    Δx = xt - branch.center
    velocity_n_m = zeros(2,3,size(branch.multipole_expansion,3))
    lamb_helmholtz = Val(false)
    u, velocity, gradient = FastMultipole.evaluate_local(Δx, branch.harmonics, velocity_n_m, branch.local_expansion, Val(expansion_order), lamb_helmholtz, DerivativesSwitch())
    return u
end

#--- RUN TESTS ---#

function solve_p_old(ε, ρ, am, pmax)
    c = ρ/am - 1.0
    p_old = ceil(-log(c,ε))
    if sign(p_old) < 0 || isinf(p_old)
        p_old = maxintfloat()
    end
    return min(Integer(p_old),pmax)
end

function get_body_positions(x_set, y_set, z_set, r_range, q_range)
    ρs = []
    θs = []
    ϕs = []
    radii = []
    qs = []
    Random.seed!(123)
    for x in x_set
        for y in y_set
            for z in z_set
                r = r_range[1] + rand() * (r_range[2]-r_range[1])
                q = q_range[1] + rand() * (q_range[2]-q_range[1])
                ρ, θ, ϕ = FastMultipole.cartesian_to_spherical(x,y,z)
                push!(ρs, ρ)
                push!(θs, θ)
                push!(ϕs, ϕ)
                push!(radii, r)
                push!(qs, q)
            end
        end
    end
    return ρs, θs, ϕs, radii, qs
end

function relative_error(u_true, u_test)
    return (u_test .- u_true) ./ u_true
end

function test_p(epsilon::Number, pmax=50;
        multipole_center = SVector{3}(0.0,0,0),
        local_center = SVector{3}(4.0,0,0),
        source_xs = [1.0],
        source_ys = [0.0],
        source_zs = [0.0],
        source_q_range = [1.0, 1.0],
        source_r_range = [0.0, 0.0],
        target_xs = [-1.0],
        target_ys = [0.0],
        target_zs = [0.0],
        target_q_range = [1.0, 1.0],
        target_r_range = [0.0, 0.0],
        shrink=true,
    )
    # define multipole (sources)
    ρs, θs, ϕs, rs, qs = get_body_positions(source_xs, source_ys, source_zs, source_r_range, source_q_range)
    multipole, multipole_system = generate_multipole(multipole_center, ρs, θs, ϕs, rs, qs, pmax, shrink)

    # define local expansion (targets)
    ρs, θs, ϕs, rs, qs = get_body_positions(target_xs, target_ys, target_zs, target_r_range, target_q_range)
    local_branch, local_system = generate_multipole(local_center, ρs, θs, ϕs, rs, qs, pmax, shrink)

    # translate for local expansion
    multipole_to_local!(local_branch, multipole, pmax)

    # determine expansion order using Pringle method
    ρ = norm(local_branch.center - multipole.center)
    ρ2 = ρ * ρ
    am = multipole.source_radius
    p_equal_spheres_check = solve_p_old(epsilon, ρ, am, pmax)

    # build expansion order
    expansion_order = FastMultipole.Dynamic(pmax, epsilon)

    # FastMultipole old error method (Equal spheres)
    rs = FastMultipole.get_r_ρ(local_branch, multipole, ρ2, FastMultipole.EqualSpheres())
    p_equal_spheres = FastMultipole.get_P(rs..., ρ, expansion_order, FastMultipole.EqualSpheres())

    # adapted for unequal spheres
    rs = FastMultipole.get_r_ρ(local_branch, multipole, ρ2, FastMultipole.UnequalSpheres())
    p_unequal_spheres = FastMultipole.get_P(rs..., ρ, expansion_order, FastMultipole.UnequalSpheres())

    if !(shrink && local_branch.source_radius > multipole.source_radius) # can break if target branch has a larger radius than source branch
        @test p_unequal_spheres <= p_equal_spheres
    end

    # adapted for unequal boxes
    rs = FastMultipole.get_r_ρ(local_branch, multipole, ρ2, FastMultipole.UnequalBoxes())
    p_unequal_boxes = FastMultipole.get_P(rs..., ρ, expansion_order, FastMultipole.UnequalBoxes())

    @test p_unequal_boxes <= p_unequal_spheres

    #--- evaluate potential ---#

    # equal spheres
    xts = (body.position for body in local_system.bodies)
    u_equal_spheres_multipole = [evaluate_multipole(xt, multipole, p_equal_spheres) for xt in xts]
    u_equal_spheres_local = [evaluate_local(xt, local_branch, p_equal_spheres) for xt in xts]

    # unequal spheres
    u_unequal_spheres_multipole = [evaluate_multipole(xt, multipole, p_unequal_spheres) for xt in xts]
    u_unequal_spheres_local = [evaluate_local(xt, local_branch, p_unequal_spheres) for xt in xts]

    # unequal boxes
    u_unequal_boxes_multipole = [evaluate_multipole(xt, multipole, p_unequal_boxes) for xt in xts]
    u_unequal_boxes_local = [evaluate_local(xt, local_branch, p_unequal_boxes) for xt in xts]

    # direct
    u_direct = [evaluate_direct(xt, multipole_system) for xt in xts]

    # prepare return values
    multipole_error_equal_spheres = relative_error(u_direct, u_equal_spheres_multipole)
    total_error_equal_spheres = relative_error(u_direct, u_equal_spheres_local)

    if !(shrink && local_branch.source_radius > multipole.source_radius) # can break if target branch has a larger radius than source branch
        @test maximum(abs.(total_error_equal_spheres)) <= epsilon
    end

    multipole_error_unequal_spheres = relative_error(u_direct, u_unequal_spheres_multipole)
    total_error_unequal_spheres = relative_error(u_direct, u_unequal_spheres_local)

    @test maximum(abs.(total_error_unequal_spheres)) <= epsilon

    multipole_error_unequal_boxes = relative_error(u_direct, u_unequal_boxes_multipole)
    total_error_unequal_boxes = relative_error(u_direct, u_unequal_boxes_local)

    @test maximum(abs.(total_error_unequal_boxes)) <= epsilon

    return p_equal_spheres, multipole_error_equal_spheres, total_error_equal_spheres, p_unequal_spheres, multipole_error_unequal_spheres, total_error_unequal_spheres, p_unequal_boxes, multipole_error_unequal_boxes, total_error_unequal_boxes
end

function test_p_set(epsilon, shrink; pmax=50,
        multipole_center = SVector{3}(0.0,0,0),
        local_center = SVector{3}(4.0,0,0),
    )


    # single source, single target: highest error
    source_xs = [1.0]
    source_ys = [0.0]
    source_zs = [0.0]
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = [-1.0]
    target_ys = [0.0]
    target_zs = [0.0]
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    # single source, single target: lowest error
    source_xs = [-1.0]
    source_ys = [0.0]
    source_zs = [0.0]
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = [1.0]
    target_ys = [0.0]
    target_zs = [0.0]
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    # single source, many targets: highest error
    source_xs = [1.0]
    source_ys = [0.0]
    source_zs = [0.0]
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = range(-1,stop=1,length=10)
    target_ys = range(-1,stop=1.0,length=10)
    target_zs = range(-1,stop=1.0,length=10)
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    # single source, many targets: lowest error
    source_xs = [-1.0]
    source_ys = [0.0]
    source_zs = [0.0]
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = range(-1,stop=1,length=10)
    target_ys = range(-1,stop=1.0,length=10)
    target_zs = range(-1,stop=1.0,length=10)
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    # many sources, single target: highest error
    source_xs = range(-1,stop=1.0,length=10)
    source_ys = range(-1,stop=1.0,length=10)
    source_zs = range(-1,stop=1.0,length=10)
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = [-1.0]
    target_ys = [0.0]
    target_zs = [0.0]
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    # many sources, single target: lowest error
    source_xs = range(-1,stop=1.0,length=10)
    source_ys = range(-1,stop=1.0,length=10)
    source_zs = range(-1,stop=1.0,length=10)
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = [1.0]
    target_ys = [0.0]
    target_zs = [0.0]
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    # many sources, many targets
    source_xs = range(-1,stop=1.0,length=10)
    source_ys = range(-1,stop=1.0,length=10)
    source_zs = range(-1,stop=1.0,length=10)
    source_q_range = [1.0, 1.0]
    source_r_range = [0.0, 0.0]
    target_xs = range(-1.0, stop=1,length=10)
    target_ys = range(-1,stop=1.0,length=10)
    target_zs = range(-1,stop=1.0,length=10)
    target_q_range = [1.0, 1.0]
    target_r_range = [0.0, 0.0]
    pes, mes, tes, pus, mus, tus, pub, mub, tub = test_p(epsilon, pmax;
        multipole_center, local_center, source_xs,  source_ys,  source_zs,  source_q_range,  source_r_range,  target_xs,  target_ys,  target_zs,  target_q_range,  target_r_range,  shrink)

    return nothing
end

@testset "dynamic expansion order" begin

test_p_set(1e-2, false)
test_p_set(1e-2, true)
test_p_set(1e-4, false)
test_p_set(1e-4, true)
test_p_set(1e-11, true)

end

function recursive_l2l(unsorted_system, tree, i_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, Pmax)
    branches = tree.branches
    branch = branches[i_branch]
    velocity_n_m = zeros(2,3,size(weights_tmp_1,3))
    if length(branch.branch_index) == 0 # leaf
        for i_sorted in branch.bodies_index
			i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
			x_target = unsorted_system[i_unsorted, FastMultipole.Position()]
			dx = x_target - branch.center
			phi, v, vg = FastMultipole.evaluate_local(dx, branch.harmonics, velocity_n_m, branch.local_expansion, Pmax, lamb_helmholtz, DerivativesSwitch())
            unsorted_system.potential[1,i_unsorted] = phi
        end
    else
        for i_child in branch.branch_index
            FastMultipole.local_to_local!(branches[i_child], branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz)
            recursive_l2l(unsorted_system, tree, i_child, weights_tmp_1, weights_tmp_2, Ts, eimϕs, ηs_mag, Hs_π2, expansion_order, lamb_helmholtz, Pmax)
        end
    end
end

function check_m2l(unsorted_system, tree, m2l_list, direct_list;
	Ts = zeros(FastMultipole.length_Ts(30)),
	eimϕs = zeros(2, 30+1),
	weights_tmp_1 = initialize_expansion(30),
	weights_tmp_2 = initialize_expansion(30),
	lamb_helmholtz = Val(false),
	harmonics = FastMultipole.initialize_harmonics(30),
)
	velocity_n_m = zeros(2,3,size(harmonics,3))

	# check all m2l_list interactions
	errors = []
	for (j_target, i_source, P) in m2l_list
		source_branch = tree.branches[i_source]
		target_branch = tree.branches[j_target]
		target_branch.local_expansion .= zero(eltype(target_branch.local_expansion))
		FastMultipole.multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ζs_mag, FastMultipole.ηs_mag, FastMultipole.Hs_π2, P, lamb_helmholtz)

		# loop over target bodies
		errs = Float64[]
		for i_sorted in target_branch.bodies_index
			i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
			x_target = unsorted_system[i_unsorted, FastMultipole.Position()]
			dx = x_target - target_branch.center
			phi, v, vg = FastMultipole.evaluate_local(dx, harmonics, velocity_n_m, target_branch.local_expansion, Val(P), lamb_helmholtz, DerivativesSwitch())

			# calculate potential directly
			phi_direct = 0.0
			for i_source_sorted in source_branch.bodies_index
				i_source_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_source_sorted, tree)
				x_source = unsorted_system[i_source_unsorted, FastMultipole.Position()]
				r = norm(x_target - x_source)
				strength = unsorted_system[i_source_unsorted, FastMultipole.Strength()]
				phi_direct += strength / r / 4 / pi
			end
			push!(errs, phi_direct - phi)
		end
		push!(errors, errs)
	end

    # repeat, but include L2L operations
	errors_l2l = []
	for (j_target, i_source, P) in m2l_list
		source_branch = tree.branches[i_source]
		target_branch = tree.branches[j_target]
		target_branch.local_expansion .= zero(eltype(target_branch.local_expansion))
		FastMultipole.multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ζs_mag, FastMultipole.ηs_mag, FastMultipole.Hs_π2, P, lamb_helmholtz)

        # downward pass manually
        unsorted_system.potential .= 0.0
        recursive_l2l(unsorted_system, tree, j_target, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ηs_mag, FastMultipole.Hs_π2, Val(P), lamb_helmholtz, tree.expansion_order)

		# loop over target bodies
		errs = Float64[]
		for i_sorted in target_branch.bodies_index
			i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
            x_target = unsorted_system[i_unsorted, FastMultipole.Position()]
            phi = unsorted_system.potential[1,i_unsorted]

			# calculate potential directly
			phi_direct = 0.0
			for i_source_sorted in source_branch.bodies_index
				i_source_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_source_sorted, tree)
				x_source = unsorted_system[i_source_unsorted, FastMultipole.Position()]
				r = norm(x_target - x_source)
				strength = unsorted_system[i_source_unsorted, FastMultipole.Strength()]
				phi_direct += strength / r / 4 / pi
			end
			push!(errs, phi_direct - phi)

		end
		push!(errors_l2l, errs)
	end

	# ensure all interactions are accounted for
	n_bodies = FastMultipole.get_n_bodies(unsorted_system)
	m2l_interactions_unsorted = zeros(Int64, n_bodies, n_bodies)

	# loop over m2l_list and check interactions
	for (i_target, i_source, P) in m2l_list
		target_branch = tree.branches[i_target]
		source_branch = tree.branches[i_source]
		for i_target_sorted in target_branch.bodies_index
			i_target_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_target_sorted, tree)
			for i_source_sorted in source_branch.bodies_index
				i_source_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_source_sorted, tree)
				m2l_interactions_unsorted[i_target_unsorted, i_source_unsorted] += 1
			end
		end
	end

	# direct interactions now
	direct_interactions_unsorted = zeros(Int64, n_bodies, n_bodies)

	# loop over direct_list and check interactions
	for (i_target, i_source) in direct_list
		target_branch = tree.branches[i_target]
		source_branch = tree.branches[i_source]
		for i_target_sorted in target_branch.bodies_index
			i_target_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_target_sorted, tree)
			for i_source_sorted in source_branch.bodies_index
				i_source_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_source_sorted, tree)
				direct_interactions_unsorted[i_target_unsorted, i_source_unsorted] += 1
			end
		end
	end

	return errors, errors_l2l, m2l_interactions_unsorted, direct_interactions_unsorted
end

# check L2L
function check_l2l(unsorted_system, tree::Tree{<:Any,P}, m2l_source, m2l_target, l2l_target) where P
    # initialize containers
	Ts = zeros(FastMultipole.length_Ts(30))
	eimϕs = zeros(2, 30+1)
	weights_tmp_1 = initialize_expansion(30)
	weights_tmp_2 = initialize_expansion(30)
	lamb_helmholtz = Val(false)
	harmonics = FastMultipole.initialize_harmonics(30)
	velocity_n_m = zeros(2,3,size(harmonics,3))

    # reset potential/expansions
    unsorted_system.potential .= 0.0
    for branch in tree.branches
        branch.local_expansion .= 0.0
    end

    # m2l
    target_branch = tree.branches[m2l_target]
    source_branch = tree.branches[m2l_source]
    l2l_target_branch = tree.branches[l2l_target]
    FastMultipole.multipole_to_local!(target_branch, source_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ζs_mag, FastMultipole.ηs_mag, FastMultipole.Hs_π2, P, lamb_helmholtz)

    # check influence
	m2l_errs = Float64[]
	for i_sorted in l2l_target_branch.bodies_index
		i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
		x_target = unsorted_system[i_unsorted, FastMultipole.Position()]
		dx = x_target - target_branch.center
		phi, v, vg = FastMultipole.evaluate_local(dx, harmonics, velocity_n_m, target_branch.local_expansion, Val(P), lamb_helmholtz, DerivativesSwitch())

		# calculate potential directly
		phi_direct = 0.0
		for i_source_sorted in source_branch.bodies_index
			i_source_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_source_sorted, tree)
			x_source = unsorted_system[i_source_unsorted, FastMultipole.Position()]
			r = norm(x_target - x_source)
			strength = unsorted_system[i_source_unsorted, FastMultipole.Strength()]
			phi_direct += strength / r / 4 / pi
		end
		push!(m2l_errs, phi_direct - phi)
	end

    # l2l
    FastMultipole.local_to_local!(l2l_target_branch, target_branch, weights_tmp_1, weights_tmp_2, Ts, eimϕs, FastMultipole.ηs_mag, FastMultipole.Hs_π2, Val(P), lamb_helmholtz)

    # check influence
	l2l_errs = Float64[]
	for i_sorted in l2l_target_branch.bodies_index
		i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
		x_target = unsorted_system[i_unsorted, FastMultipole.Position()]
		dx = x_target - l2l_target_branch.center
		phi, v, vg = FastMultipole.evaluate_local(dx, harmonics, velocity_n_m, l2l_target_branch.local_expansion, Val(P), lamb_helmholtz, DerivativesSwitch())

		# calculate potential directly
		phi_direct = 0.0
		for i_source_sorted in source_branch.bodies_index
			i_source_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_source_sorted, tree)
			x_source = unsorted_system[i_source_unsorted, FastMultipole.Position()]
			r = norm(x_target - x_source)
			strength = unsorted_system[i_source_unsorted, FastMultipole.Strength()]
			phi_direct += strength / r / 4 / pi
		end
		push!(l2l_errs, phi_direct - phi)
	end

    # return result
    return m2l_errs, l2l_errs
end

function check_m2l_again(unsorted_system, tree, i_source, i_target, expansion_order; redo_b2m=false, check_b2m=false)
    # get branches
    source_branch = tree.branches[i_source]
    target_branch = tree.branches[i_target]

    # b2m (if requested)
    if redo_b2m
        old_multipole = deepcopy(source_branch.multipole_expansion)
        source_branch.multipole_expansion .= 0.0
        for i_sorted in source_branch.bodies_index
            i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)

            # relative body position
            x = unsorted_system[i_unsorted,Position()]
            Δx = x - source_branch.center

            # update values
            FastMultipole.body_to_multipole_point!(Point{Source}, source_branch.multipole_expansion, source_branch.harmonics, Δx, unsorted_system[i_unsorted, Strength()], Val(expansion_order))
        end
        new_multipole = deepcopy(source_branch.multipole_expansion)
    end

    if check_b2m
        for i_sorted in target_branch.bodies_index
            i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
            x_target = unsorted_system[i_unsorted,Position()]

            phi_multipole, _, _ = evaluate_multipole(x_target, source_branch.center, source_branch.multipole_expansion, DerivativesSwitch(), Val(expansion_order))

            phi_direct = 0.0
            for i_sorted_source in source_branch.bodies_index
                i_unsorted_source = FastMultipole.sorted_index_2_unsorted_index(i_sorted_source, tree)
                x_source = unsorted_system[i_unsorted_source, Position()]
                strength = unsorted_system[i_unsorted_source, Strength()]
                r = norm(x_target - x_source)
                phi_direct += strength / 4 / pi / r
            end

        end
    end

    # reset local expansion
    target_branch.local_expansion .= 0.0

    # m2l
    multipole_to_local!(target_branch, source_branch, expansion_order)

    # evaluate local
    phis = Float64[]
    phis_direct = Float64[]
    for i_sorted in target_branch.bodies_index
        i_unsorted = FastMultipole.sorted_index_2_unsorted_index(i_sorted, tree)
        x_target = unsorted_system[i_unsorted, Position()]

        phi = evaluate_local(x_target, target_branch, expansion_order)
        push!(phis, phi)

        phi_direct = 0.0
        for i_sorted_source in source_branch.bodies_index
            i_unsorted_source = FastMultipole.sorted_index_2_unsorted_index(i_sorted_source, tree)
            x_source = unsorted_system[i_unsorted_source, Position()]
            strength = unsorted_system[i_unsorted_source, Strength()]
            r = norm(x_target - x_source)
            phi_direct += strength / 4 / pi / r
        end
        push!(phis_direct, phi_direct)
    end

    return phis, phis_direct
end

function analyze_m2l_list(m2l_list)
    Ps = [m2l[3] for m2l in m2l_list]
    return length(Ps), mean(Ps), std(Ps)
end

@testset "dynamic expansion order: UnequalSpheres without shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UnequalSpheres()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
#FastMultipole.DEBUG[] = true
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=false, unsort_bodies=true)
#FastMultipole.DEBUG[] = false
FastMultipole.visualize("bad_error", masses_fmm, tree; toggle_branches=true, toggle_bodies=true)
u_fmm = masses_fmm.potential[1,:]

Ps = analyze_m2l_list(m2l_list)
@show Ps

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

@test maximum(abs.(err)) < rtol

@show maximum(abs.(err))
end

@testset "dynamic expansion order: UnequalSpheres with shrinking" begin

n_bodies = 10000
Pmax = 18
rtol = 1e-4
error_method = FastMultipole.UnequalSpheres()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
FastMultipole.DEBUG[] = true
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=true, unsort_bodies=true)
FastMultipole.DEBUG[] = false
u_fmm = masses_fmm.potential[1,:]

Ps = analyze_m2l_list(m2l_list)
@show Ps

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

@test maximum(abs.(err)) < rtol
@show maximum(abs.(err))

end

#=
@testset "dynamic expansion order: UniformUnequalSpheres without shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UniformUnequalSpheres()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=false, unsort_bodies=true)
u_fmm = masses_fmm.potential[1,:]

Ps = analyze_m2l_list(m2l_list)
@show Ps

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

@test maximum(abs.(err)) < rtol

@show maximum(abs.(err))

end
=#

@testset "dynamic expansion order: UniformUnequalSpheres with shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UniformUnequalSpheres()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=true, unsort_bodies=true)
u_fmm = masses_fmm.potential[1,:]

Ps = analyze_m2l_list(m2l_list)
@show Ps

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

@test maximum(abs.(err)) < rtol
@show maximum(abs.(err))

end

@testset "dynamic expansion order: UnequalBoxes without shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UnequalBoxes()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=false, unsort_bodies=true)
u_fmm = masses_fmm.potential[1,:]

Ps = analyze_m2l_list(m2l_list)
@show Ps

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

@test maximum(abs.(err)) < rtol
@show maximum(abs.(err))

end

@testset "dynamic expansion order: UnequalBoxes with shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UnequalBoxes()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
#FastMultipole.DEBUG[] = true
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=true, unsort_bodies=true)

Ps = analyze_m2l_list(m2l_list)
@show Ps
#FastMultipole.DEBUG[] = false
#FastMultipole.visualize("bad_error", masses_fmm, tree; toggle_branches=true, toggle_bodies=true)
u_fmm = masses_fmm.potential[1,:]

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

#=
errs, errs_l2l, m2l, direct = check_m2l(masses_fmm, tree, m2l_list, direct_list)


m2l_source, m2l_target, l2l_target = 11, 18, 18
m2l_errs, l2l_errs = check_l2l(masses_fmm, tree, m2l_source, m2l_target, l2l_target)

# check m2l from 11 to 18 manually
i_source, i_target = 11, 18
#i_source, i_target = 18, 11
phis_again, phis_direct_again = check_m2l_again(masses_fmm, tree, i_source, i_target, Pmax; redo_b2m=false, check_b2m=true)
=#

@test maximum(abs.(err)) < rtol
@show maximum(abs.(err))

end

#=
@testset "dynamic expansion order: UniformUnequalBoxes without shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UniformUnequalBoxes()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=false, unsort_bodies=true)
u_fmm = masses_fmm.potential[1,:]

Ps = analyze_m2l_list(m2l_list)
@show Ps

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

@test maximum(abs.(err)) < rtol
@show maximum(abs.(err))

end
=#

#=
@testset "dynamic expansion order: UniformUnequalBoxes with shrinking" begin

n_bodies = 10000
Pmax = 10
rtol = 1e-4
error_method = FastMultipole.UniformUnequalBoxes()
leaf_size = 100
multipole_threshold=0.5
seed = 123

masses_fmm = generate_gravitational(seed, n_bodies; radius_factor=0.1)
masses_direct = generate_gravitational(seed, n_bodies; radius_factor=0.1)
expansion_order = Dynamic(Pmax,rtol)
#FastMultipole.DEBUG[] = true
tree, m2l_list, direct_list, derivatives_switch = fmm.fmm!(masses_fmm; expansion_order, error_method, leaf_size, multipole_threshold, shrink_recenter=true, unsort_bodies=true)

Ps = analyze_m2l_list(m2l_list)
@show Ps
#FastMultipole.DEBUG[] = false
#FastMultipole.visualize("bad_error", masses_fmm, tree; toggle_branches=true, toggle_bodies=true)
u_fmm = masses_fmm.potential[1,:]

# direct
masses_fmm.potential .= 0.0
fmm.direct!(masses_direct)
u_direct = masses_direct.potential[1,:]

err = relative_error(u_direct, u_fmm)

#=
errs, errs_l2l, m2l, direct = check_m2l(masses_fmm, tree, m2l_list, direct_list)


m2l_source, m2l_target, l2l_target = 11, 18, 18
m2l_errs, l2l_errs = check_l2l(masses_fmm, tree, m2l_source, m2l_target, l2l_target)

# check m2l from 11 to 18 manually
i_source, i_target = 11, 18
#i_source, i_target = 18, 11
phis_again, phis_direct_again = check_m2l_again(masses_fmm, tree, i_source, i_target, Pmax; redo_b2m=false, check_b2m=true)
=#

@test maximum(abs.(err)) < rtol
@show maximum(abs.(err))

end
=#
