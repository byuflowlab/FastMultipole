# @testset "3D vortex particles" begin
using LinearAlgebra
include("../test/vortex.jl")

function vortex_ring(origin; orientation=[1.0 0 0;0 1 0; 0 0 1],
        azimuth_radius=1.0, torroid_radius=1/7,
        n_azimuth=5, n_rho=1,
        zeta_azimuth=1.3, zeta_torroid=1.3,
        circulation=1.0
    )
    sigma = 2 * pi * azimuth_radius / n_azimuth * zeta_azimuth

    zhat = [0,0,1.0]
    n_per_azimuth = 0
    if n_rho == 1
        torroid_radius *= 0
        n_per_azimuth += 1
    else
        for rho in range(0, stop=torroid_radius, length=n_rho)
            for _ in 1:ceil(2*pi*rho*zeta_torroid/sigma)
                n_per_azimuth += 1
            end
        end
    end

    n_particles = n_per_azimuth * (n_azimuth-1)
    particle_strength = circulation / n_particles#/ (4/3*pi*sigma^3 * n_particles)

    azimuth_locations = zeros(3,n_per_azimuth)
    azimuth_transformed = zeros(3,n_per_azimuth)
    particles = zeros(i_STRENGTH_vortex[end],n_particles)
    particles_transformed = zeros(i_STRENGTH_vortex[end],n_particles)
    for (itheta,theta) in enumerate(range(0,stop=2*pi, length=n_azimuth)[1:end-1])
        azimuth_origin = azimuth_radius .* [cos(theta), sin(theta), 0.0]
        xhat = azimuth_origin / azimuth_radius
        yhat = cross(zhat,xhat)
        azimuth_orientation = hcat(xhat, yhat, zhat)
        i_tape = 1
        for rho in range(0, stop=torroid_radius, length=n_rho)
            for phi in range(0,stop=2*pi,length=Int(ceil(2*pi*rho*zeta_torroid/sigma)))
                azimuth_locations[1,i_tape] = rho * cos(phi)
                azimuth_locations[3,i_tape] = rho * sin(phi)
                i_tape += 1
            end
        end
        mul!(azimuth_transformed, azimuth_orientation, azimuth_locations)
        azimuth_transformed .+= azimuth_origin
        particles[i_POSITION_vortex,(itheta-1)*n_per_azimuth+1:itheta*n_per_azimuth] .+= azimuth_transformed
        particles_transformed[i_STRENGTH_vortex,(itheta-1)*n_per_azimuth+1:itheta*n_per_azimuth] .+= yhat * particle_strength
    end

    # rotate and translate to new origin/orientation
    mul!(view(particles_transformed,i_POSITION_vortex,:), orientation, particles[i_POSITION_vortex,:])
    particles_transformed[i_POSITION_vortex,:] .+= origin
    return particles_transformed
end

# function add_vortex_ring!(vortex_particles, origin; orientation=[1.0 0 0;0 1 0; 0 0 1],
#     azimuth_radius=1.0, torroid_radius=1/7,
#     n_azimuth=14, n_rho=1,
#     zeta_azimuth=1.3, zeta_torroid=1.3,
#     circulation=1.0
# )

# end

# test
# function test_leapfrog()
circulation = 1
n_azimuth = 13
ring_1 = vortex_ring(zeros(3); circulation, n_azimuth)
ring_2 = vortex_ring([0.0,0.0,1.0]; circulation, n_azimuth)

# vortex_particles = VortexParticles(ring_1)
vortex_particles = VortexParticles(hcat(ring_1,ring_2))

expansion_order = 20
n_per_branch = 5
theta = 4
# # fmm.direct!(vortex_particles)#, expansion_order, n_per_branch, theta)
# tree = fmm.fmm!(vortex_particles, fmm.Options(expansion_order, n_per_branch, theta))

# update_velocity_stretching!(vortex_particles)
# save_vtk("vortex_rings.$i", vortex_particles)

# vortex_particles.potential .*= 0
# vortex_particles.velocity_stretching .*= 0

dt = 1.0
fmm_options = fmm.Options(expansion_order, n_per_branch, theta)
integrate! = Euler(dt)
nsteps = 10

convect!(vortex_particles, nsteps;
        # integration options
        integrate!,
        # save options
        save=true, filename="leaping_vortex_rings", compress=false,
        # fmm options
        fmm_options
    )
