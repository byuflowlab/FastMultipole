# @testset "3D vortex particles" begin
using LinearAlgebra

function vortex_ring(origin, orientation, centerline_radius, torroid_radius, n_azimuth, n_rho, sigma, zeta)
    # zhat = [0,0,1.0]
    for theta in range(0,stop=2*pi, length=n_azimuth)
        azimuth_origin = centerline_radius .* [cos(theta), sin(theta), 0.0]
        # azimuth_y = cross(zhat, azimuth_origin)
        # azimuth_orientation = hcat(azimuth_orientation, azimuth_y, zhat)
        for rho in range(0, stop=torroid_radius, length=n_rho)
            for phi in range(0,stop=2*pi,length=ceil(2*pi*rho*zeta/sigma))
                # TODO: finish vortex ring
            end
        end
    end

end
