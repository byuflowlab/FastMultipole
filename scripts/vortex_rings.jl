# @testset "3D vortex particles" begin
using LinearAlgebra
include("../test/vortex.jl")
include("vtk.jl")

function vortex_ring(origin; orientation=[1.0 0 0;0 1 0; 0 0 1],
        azimuth_radius=1.0, torroid_radius=1/7,
        n_azimuth=14, n_rho=1,
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

    n_particles = n_per_azimuth * n_azimuth
    particle_strength = circulation / (4/3*pi*sigma^3 * n_particles)

    azimuth_locations = zeros(3,n_per_azimuth)
    azimuth_transformed = zeros(3,n_per_azimuth)
    particles = zeros(7,n_particles)
    particles_transformed = similar(particles)
    for (itheta,theta) in enumerate(range(0,stop=2*pi, length=n_azimuth))
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
        particles[1:3,(itheta-1)*n_per_azimuth+1:itheta*n_per_azimuth] .+= azimuth_transformed
        particles_transformed[5:7,(itheta-1)*n_per_azimuth+1:itheta*n_per_azimuth] .+= yhat * particle_strength
    end

    # rotate and translate to new origin/orientation
    mul!(view(particles_transformed,1:3,:), orientation, view(particles,1:3,:))
    particles_transformed[1:3,:] .+= origin
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
ring_1 = vortex_ring(zeros(3))
ring_2 = vortex_ring([0.0,0.0,1.0])

vortex_particles = VortexParticles(hcat(ring_1,ring_2))

save_vtk("test_vortex_rings", vortex_particles)
