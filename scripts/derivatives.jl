using Test

@testset "derivatives" begin
"""
dr_k/dx_idx_j
"""
function d2rdx2(r, theta, phi)
    derivatives = zeros(3,3,3)
    derivatives[:,:,1] .= [
        (1-cos(phi)^2 * sin(theta)^2)/r -sin(theta)^2*cos(phi)*sin(phi)/r -sin(theta)*cos(phi)*cos(theta)/r;
        (-sin(theta)^2*cos(phi)*sin(phi))/r (1-sin(theta)^2*sin(phi)^2)/r -sin(theta)*sin(phi)*cos(theta)/r;
        -sin(theta)*cos(phi)*cos(theta)/r -sin(theta)*sin(phi)*cos(theta)/r sin(theta)^2/r
    ]
    derivatives[:,:,2] .= [
        cos(theta)/sin(theta)*(1-cos(phi)^2*(1+2*sin(theta)^2))/r^2 -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/r^2 cos(phi)*(1-2*cos(theta)^2)/r^2;
        -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/r^2 cos(theta)/sin(theta)*(1-sin(phi)^2*(1+2*sin(theta)^2))/r^2 (2*sin(theta)^2-1)/r^2*sin(phi);
        cos(phi)*(1-2*cos(theta)^2)/r^2 (2*sin(theta)^2-1)/r^2*sin(phi) 2*sin(theta)*cos(theta)/r^2
    ]
    derivatives[:,:,3] .= [
        2*cos(phi)*sin(phi)/r^2/sin(theta)^2 (2*sin(phi)^2-1)/r^2/sin(theta)^2 0;
        (2*sin(phi)^2-1)/r^2/sin(theta)^2 -2*sin(phi)*cos(phi)/r^2/sin(theta)^2 0;
        0 0 0
    ]
    return derivatives
end

function cartesian_2_spherical(x,y,z; vec=false)
    r = sqrt(x^2+y^2+z^2)
    theta = acos(z/r)
    phi = atan(y,x)
    vec && return [r,theta,phi]
    return r, theta, phi
end

# function spherical_2_cartesian(r,theta,phi)
#     x = r*sin(theta)*cos(phi)
#     y = r*sin(theta)*sin(phi)
#     z = r*cos(theta)
#     return x,y,z
# end

function d2rdx2_cart(x,y,z)
    r, theta, phi = cartesian_2_spherical(x,y,z)
    return d2rdx2(r,theta,phi)
end

"""
dr_j/dx_i
"""
function drdx(r,theta,phi)
    derivatives = [
        sin(theta)*cos(phi) cos(theta)*cos(phi)/r -sin(phi)/r/sin(theta);
        sin(theta)*sin(phi) cos(theta)*sin(phi)/r cos(phi)/r/sin(theta);
        cos(theta) -sin(theta)/r 0
    ]
    # derivatives = [
    #     sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
    #     cos(theta)*cos(phi)/r cos(theta)*sin(phi)/r -sin(theta)/r;
    #     -sin(phi)/r/sin(theta) cos(phi)/r/sin(theta) 0
    # ]
    return derivatives
end

function drdx_cart(x,y,z)
    r, theta, phi = cartesian_2_spherical(x,y,z)
    return drdx(r, theta, phi)
end

function d2rdx2_fd(x,y,z)
    derivatives = zeros(3,3,3)
    derivatives[1,:,:] .= (drdx_cart(x+1e-6,y,z) - drdx_cart(x,y,z))/1e-6
    derivatives[2,:,:] .= (drdx_cart(x,y+1e-6,z) - drdx_cart(x,y,z))/1e-6
    derivatives[3,:,:] .= (drdx_cart(x,y,z+1e-6) - drdx_cart(x,y,z))/1e-6
    return derivatives
end

x,y,z = rand(3)

function drdx_fd(x,y,z)
    derivatives = zeros(3,3)
    derivatives[1,:] .= (cartesian_2_spherical(x+1e-6,y,z;vec=true) - cartesian_2_spherical(x,y,z;vec=true))/1e-6
    derivatives[2,:] .= (cartesian_2_spherical(x,y+1e-6,z;vec=true) - cartesian_2_spherical(x,y,z;vec=true))/1e-6
    derivatives[3,:] .= (cartesian_2_spherical(x,y,z+1e-6;vec=true) - cartesian_2_spherical(x,y,z;vec=true))/1e-6
    return derivatives
end

fd_1 = drdx_fd(x,y,z)
anal_1 = drdx_cart(x,y,z)

for i in 1:length(fd_1)
    @test isapprox(fd_1[i], anal_1[i]; atol=1e-4)
end

fd = d2rdx2_fd(x,y,z)
anal = d2rdx2_cart(x,y,z)

for i in 1:length(fd)
    @test isapprox(fd[i], anal[i]; atol=1e-4)
end
end

# chain rule
"""
dr_k/dx_idx_j
"""
function d2rdx2(r, theta, phi)
    derivatives = zeros(3,3,3)
    derivatives[:,:,1] .= [
        (1-cos(phi)^2 * sin(theta)^2)/r -sin(theta)^2*cos(phi)*sin(phi)/r -sin(theta)*cos(phi)*cos(theta)/r;
        (-sin(theta)^2*cos(phi)*sin(phi))/r (1-sin(theta)^2*sin(phi)^2)/r -sin(theta)*sin(phi)*cos(theta)/r;
        -sin(theta)*cos(phi)*cos(theta)/r -sin(theta)*sin(phi)*cos(theta)/r sin(theta)^2/r
    ]
    derivatives[:,:,2] .= [
        cos(theta)/sin(theta)*(1-cos(phi)^2*(1+2*sin(theta)^2))/r^2 -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/r^2 cos(phi)*(1-2*cos(theta)^2)/r^2;
        -cos(theta)/sin(theta)*sin(phi)*cos(phi)*(1+2*sin(theta)^2)/r^2 cos(theta)/sin(theta)*(1-sin(phi)^2*(1+2*sin(theta)^2))/r^2 (2*sin(theta)^2-1)/r^2*sin(phi);
        cos(phi)*(1-2*cos(theta)^2)/r^2 (2*sin(theta)^2-1)/r^2*sin(phi) 2*sin(theta)*cos(theta)/r^2
    ]
    derivatives[:,:,3] .= [
        2*cos(phi)*sin(phi)/r^2/sin(theta)^2 (2*sin(phi)^2-1)/r^2/sin(theta)^2 0;
        (2*sin(phi)^2-1)/r^2/sin(theta)^2 -2*sin(phi)*cos(phi)/r^2/sin(theta)^2 0;
        0 0 0
    ]
    return derivatives
end
"""
dr_j/dx_i
"""
function drdx(r,theta,phi)
    derivatives = [
        sin(theta)*cos(phi) cos(theta)*cos(phi)/r -sin(phi)/r/sin(theta);
        sin(theta)*sin(phi) cos(theta)*sin(phi)/r cos(phi)/r/sin(theta);
        cos(theta) -sin(theta)/r 0
    ]
    # derivatives = [
    #     sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
    #     cos(theta)*cos(phi)/r cos(theta)*sin(phi)/r -sin(theta)/r;
    #     -sin(phi)/r/sin(theta) cos(phi)/r/sin(theta) 0
    # ]
    return derivatives
end

# """
# dr_i/dx_j
# """
# function drdx(r,theta,phi)
#     # derivatives = [
#     #     sin(theta)*cos(phi) cos(theta)*cos(phi)/r -sin(phi)/r/sin(theta);
#     #     sin(theta)*sin(phi) cos(theta)*sin(theta)/r cos(phi)/r/sin(theta);
#     #     cos(theta) -sin(theta)/r 0
#     # ]
#     derivatives = [
#         sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta);
#         cos(theta)*cos(phi)/r cos(theta)*sin(phi)/r -sin(theta)/r;
#         -sin(phi)/r/sin(theta) cos(phi)/r/sin(theta) 0
#     ]
#     return derivatives
# end

function derivatives_2_cartesian(first_derivatives, second_derivatives, r, theta, phi)
    @assert size(spherical_derivatives) == (3,3)
    d2_unit_vec = d2rdx2(r, theta, phi)
    d_unit_vec = drdx(r,theta,phi)
    cartesian_derivatives = zeros(3,3)
    for k in 1:3
        cartesian_derivatives .+= d2_unit_vec[:,:,k] .* first_derivatives[k]
    end
    cartesian_derivatives .+= d_unit_vec * second_derivatives * transpose(d_unit_vec)
    return cartesian_derivatives
end

scalar_func(r,theta,phi) = 1/r*sin(theta)*cos(theta) + sin(phi)^2

function fd_derivative(func, r, theta, phi, i; step=1e-8)
    mask = zeros(3)
    mask[i] = step
    return (func.([r,theta,phi] + mask...) - func.(r,theta,phi))/step
end

simple(r,theta,phi) = 3*r

vec_func(r,theta,phi) = [3*r, 4*theta, 5*phi]

fd_derivative(simple, 1.0,1.0,1.0,1)
fd_derivative(vec_func, 1.0,1.0,1.0,2)

function first_derivatives(scalar_func, r, theta, phi)
    first_derivatives = [
        fd_derivative(scalar_func, r, theta, phi, i) for i in 1:3
    ]
    return first_derivatives
end


function second_derivatives(scalar_func, r, theta, phi)
    second_derivatives = zeros(3,3)
    for i in 1:3
        second_derivatives
    end
end