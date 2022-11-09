using Test
import FLOWFMM
fmm = FLOWFMM

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
    @test isapprox(fd[i], anal[i]; atol=1e-3)
end
end

@testset "chain rule" begin

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

function cs_derivative(func, i, args; step=1e-25)
    args[i] += step*im
    return imag(func(args))/step
end

function simple(X)
    r = X[1]
    theta = X[2]
    phi = X[3]
    return 3*r^4*cos(theta) + 2*theta - sin(phi)
end

function cartesian_2_spherical(x,y,z; vec=false)
    r = sqrt(x^2+y^2+z^2)
    theta = acos(z/r)
    phi = atan(y/x)
    vec && return [r,theta,phi]
    return r, theta, phi
end

function spherical_2_cartesian(r,theta,phi; vec=false)
    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
    vec && return [x,y,z]
    return x,y,z
end

function simple_cart(X)
    args = cartesian_2_spherical(X...; vec=true)
    return simple(args)
end

r0, theta0, phi0 = rand(3)

dsdr(r,theta,phi) = 12 * r^3 * cos(theta)
dsdt(r,theta,phi) = 3*r^4*-sin(theta) + 2
dsdp(r,theta,phi) = -cos(phi)
d2sdr2(r,theta,phi) = 36 * r^2 * cos(theta)
d2sdrdt(r,theta,phi) = -12*r^3*sin(theta)
d2sdrdp(r,theta,phi) = 0.0
d2sdt2(r,theta,phi) = -3*r^4*cos(theta)
d2sdtdp(r,theta,phi) = 0.0
d2sdp2(r,theta,phi) = sin(phi)

testd = cs_derivative(simple,1,Complex{Float64}[r0,theta0,phi0])
@show testd - dsdr(r0,theta0,phi0)

function second_derivative(func,i,j,args; step=1e-8)
    first_func(args) = cs_derivative(func,i,args)
    mask = zeros(length(args))
    mask[j] += step
    sec = (first_func(args + mask) - first_func(args))/step
    return sec
end

function fd_hessian(func, args)
    simple_jacobian = zeros(3,3)
    for j in 1:3
        for i in 1:3
            simple_jacobian[j,i] = second_derivative(func,i,j,convert(Vector{Complex{Float64}},args))
        end
    end
    return simple_jacobian
end
simple_jacobian = fd_hessian(simple,[r0,theta0,phi0])
simple_jacobian_anal = [
    d2sdr2(r0,theta0,phi0) d2sdrdt(r0,theta0,phi0) d2sdrdp(r0,theta0,phi0);
    d2sdrdt(r0,theta0,phi0) d2sdt2(r0,theta0,phi0) d2sdtdp(r0,theta0,phi0);
    d2sdrdp(r0,theta0,phi0) d2sdtdp(r0,theta0,phi0) d2sdp2(r0,theta0,phi0);
]

for i in 1:length(simple_jacobian)
    @test isapprox(simple_jacobian[i], simple_jacobian_anal[i]; atol=1e-6)
end

# test chain rule
args = [r0,theta0,phi0]
spherical_grad = [cs_derivative(simple,i,convert(Vector{Complex{Float64}},args)) for i in 1:3]
spherical_hessian = fd_hessian(simple, args)
d2rkdxidxj = d2rdx2(args...)
drjdxi = drdx(args...)

cartesian_hessian = zeros(3,3)
for k in 1:3
    cartesian_hessian .+= d2rkdxidxj[:,:,k] * spherical_grad[k]
end

cartesian_hessian .+= drjdxi * spherical_hessian * transpose(drjdxi)

cartesian_hessian_fd = fd_hessian(simple_cart, spherical_2_cartesian(args...;vec=true))

for i in 1:length(cartesian_hessian)
    @test isapprox(cartesian_hessian[i], cartesian_hessian_fd[i]; atol=1e-4)
end

# now test FMM function
potential_hessian = zeros(3,3,4)
potential_jacobian = zeros(3,4)
for i in 1:4
    potential_hessian[:,:,i] .= spherical_hessian
    potential_jacobian[:,i] .= spherical_grad
end
workspace = zeros(3,4)

fmm.spherical_2_cartesian!(potential_jacobian, potential_hessian, workspace, r0, theta0, phi0)

for i in 1:3
    for j in 1:3
        @test isapprox(cartesian_hessian[i,j], potential_hessian[i,j,1]; atol=1e-12)
    end
end
end
