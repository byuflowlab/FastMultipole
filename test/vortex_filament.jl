struct VortexFilaments{TF}
    x::Matrix{SVector{3,TF}}
    strength::Vector{SVector{3,TF}}
    potential::Vector{TF}
    force::Vector{SVector{3,TF}}
    gradient::Vector{SMatrix{3,3,TF,9}}
end

function VortexFilaments(x, strength::Vector{SVector{3,TF}};
        potential = zeros(size(x,2)),
        force = zeros(SVector{3,TF},size(x,2)),
        gradient = zeros(SMatrix{3,3,TF,9},size(x,2))
    ) where TF

    return VortexFilaments(x, strength, potential, force, gradient)
end

function generate_vortex_filaments(ntheta; nrings=2, r=1.0, dz=1.0, strength=1e-2)
    pts = zeros(SVector{3,Float64}, 2, (ntheta-1)*nrings)
    strength_vec = zeros(SVector{3,Float64}, size(pts,2))

    # create rings
    i_filament = 1
    for z in range(0, stop=dz, length=nrings)
        x1 = SVector{3}(1.0, 0.0, z)
        for theta in range(0,stop=2*pi, length=ntheta)
            if theta > 10*eps()
                # filament
                x2 = SVector{3}(cos(theta)*r, sin(theta)*r, z)
                pts[1,i_filament] = x1
                pts[2,i_filament] = x2

                # determine strength
                Γ = x2 - x1
                Γ = Γ * strength / norm(Γ)
                strength_vec[i_filament] = Γ

                # recurse
                x1 = x2
                i_filament += 1
            end
        end
    end

    # create filaments
    potential = zeros(length(strength_vec))
    force = zeros(SVector{3,Float64}, length(strength_vec))
    gradient = zeros(SMatrix{3,3,Float64,9}, length(strength_vec))

    return VortexFilaments(pts, strength_vec, potential, force, gradient)
end

function vortex_filament(x1,x2,xt,q)
    r1 = xt - x1
    r2 = xt - x2

    nr1 = norm(r1)
    nr2 = norm(r2)

    num = cross(r1, r2)
    norm_num = norm(num)
    if norm_num < eps() * 10
        return zero(SVector{3,eltype(x1)})
    end

    f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
    f2 = (1/nr1 + 1/nr2)

    V = (f1*f2)/(4*pi) * norm(q)

    return V
end

function reset!(filaments::VortexFilaments)
    filaments.potential .= zero(eltype(filaments.potential))
    for i in eachindex(filaments.potential)
        filaments.force[i] = SVector{3,Float64}(0.0,0,0)
        filaments.gradient[i] = zero(SMatrix{3,3,Float64,9})
    end
end

#------- compatibility functions -------#

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::VortexFilaments, i_body)
    buffer[1:3,i_buffer] .= (system.x[1,i_body]  + system.x[2,i_body]) * 0.5
    buffer[4,i_buffer] = norm(system.x[2,i_body] - system.x[1,i_body]) * 0.5
    buffer[5:7,i_buffer] .= system.strength[i_body]
    buffer[8:10,i_buffer] .= system.x[1,i_body]
    buffer[11:13,i_buffer] .= system.x[2,i_body]
end

function FastMultipole.strength_dims(system::VortexFilaments)
    return 3
end

FastMultipole.get_n_bodies(system::VortexFilaments) = length(system.strength)

FastMultipole.data_per_body(::VortexFilaments) = 13

Base.eltype(::VortexFilaments{TF}) where TF = TF

FastMultipole.get_position(system::VortexFilaments, i) = (system.x[1,i] + system.x[2,i]) * 0.5

function FastMultipole.reset!(system::VortexFilaments)
    system.potential .= zero(eltype(system.potential))
    system.force .= zero(eltype(system.force))
    system.gradient .= zero(eltype(system.gradient))
end

function FastMultipole.direct!(target_system, target_index, derivatives_switch::DerivativesSwitch{PS,VS,GS}, source_system::VortexFilaments, source_buffer, source_index) where {PS,VS,GS}
    for i_source in source_index
        x1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
        x2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
        q = FastMultipole.get_strength(source_buffer, source_system, i_source)

        for i_target in target_index
            xt = FastMultipole.get_position(target_system, i_target)

            if VS
                v = vortex_filament(x1,x2,xt,q)
                FastMultipole.set_velocity!(target_system, i_target, v)
            end
        end
    end
end

FastMultipole.body_to_multipole!(system::VortexFilaments, args...) = FastMultipole.body_to_multipole!(Filament{Vortex}, system, args...)

function FastMultipole.buffer_to_target_system!(target_system::VortexFilaments, i_target, derivatives_switch, target_buffer, i_buffer)

    # extract from buffer
    velocity = FastMultipole.get_velocity(target_buffer, i_buffer)

    # load into system
    target_system.force[i_target] += velocity

end

#------- test -------#
#=
filaments = generate_vortex_filaments(100; nrings=50, r=1.0, dz=10.0, strength=1e-2)

direct!(filaments)
v_direct = deepcopy(filaments.force)

reset!(filaments)
optargs, target_tree, source_tree, m2l_list, direct_list, derivatives_switches, error_success = fmm!(filaments; leaf_size_source=10, multipole_threshold=0.5, expansion_order=10)

v_fmm = deepcopy(filaments.force)

@show maximum(norm.(v_fmm - v_direct))
=#
