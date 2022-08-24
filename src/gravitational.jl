#####
##### gravitational elements
#####
function get_x end

struct Mass{TF} <: Element{TF}
    X::Vector{TF} # location of each element
    m::Vector{TF} # mass of the element
    V::Vector{TF} # potential induced by other elements
end

function Mass(X::Vector{TF}, m::TF, V::TF=0.0) where TF
    Mass(X, [m], [V])
end

function Mass{TF}(dims) where TF
    Mass(zeros(dims), zeros(1), zeros(1))
end

function get_X(mass::Mass)
    return mass.X
end

function get_q(mass::Mass)
    return mass.m[1]
end

function get_V(mass::Mass)
    return mass.V[1]
end

function set_X(mass::Mass, new_X)
    mass.X .= new_X
end

function set_q(mass::Mass, new_q)
    mass.m[1] = new_q
end

function add_q(mass::Mass, summand)
    set_q(mass, get_q(mass) + summand)
end

function set_V(mass::Mass, new_V)
    mass.V[1] = new_V
end

function add_V(mass::Mass, summand)
    set_V(mass, get_V(mass) + summand)
end

#####
##### gravitational kernel
#####
# function kernel!(target::Mass, source::Mass)

# end


struct Gravitational{TF,dims} <: Kernel{TF,dims}
    parameters::TF # Universal gravitational constant
    potential_derivatives::Array{Function,dims}
    order::UInt8
end

function Gravitational(dims=3; G=1.0, order=4)
    @assert dims == 3 || dims == 2 "requested Gravitational kernel of dimensions $dims which is not yet supported"
    potential_derivatives = Array{Function, dims}(undef, (order+1) * ones(UInt8,dims)...) # default through 4th derivatives
    if dims == 3
        # raw potential
        potential_derivatives[1,1,1] = gravitational_potential_3D
        # first derivatives
        potential_derivatives[2,1,1] = gravitational_dx
        potential_derivatives[1,2,1] = gravitational_dy
        potential_derivatives[1,1,2] = gravitational_dz
        # second derivatives
        potential_derivatives[3,1,1] = gravitational_dx2
        potential_derivatives[1,3,1] = gravitational_dy2
        potential_derivatives[1,1,3] = gravitational_dz2
        potential_derivatives[2,2,1] = gravitational_dxdy
        potential_derivatives[1,2,2] = gravitational_dydz
        potential_derivatives[2,1,2] = gravitational_dxdz
        # third derivatives
        potential_derivatives[4,1,1] = gravitational_dx3
        potential_derivatives[1,4,1] = gravitational_dy3
        potential_derivatives[1,1,4] = gravitational_dz3
        potential_derivatives[3,2,1] = gravitational_dx2dy
        potential_derivatives[3,1,2] = gravitational_dx2dz
        potential_derivatives[1,3,2] = gravitational_dy2dz
        potential_derivatives[2,3,1] = gravitational_dxdy2
        potential_derivatives[2,1,3] = gravitational_dxdz2
        potential_derivatives[1,2,3] = gravitational_dydz2
        potential_derivatives[2,2,2] = gravitational_dxdydz
        # fourth derivatives
        potential_derivatives[5,1,1] = gravitational_dx4
        potential_derivatives[1,5,1] = gravitational_dy4
        potential_derivatives[1,1,5] = gravitational_dz4
        potential_derivatives[4,2,1] = gravitational_dx3dy
        potential_derivatives[4,1,2] = gravitational_dx3dz
        potential_derivatives[1,4,2] = gravitational_dy3dz
        potential_derivatives[2,4,1] = gravitational_dxdy3
        potential_derivatives[2,1,4] = gravitational_dxdz3
        potential_derivatives[1,2,4] = gravitational_dydz3
        potential_derivatives[3,3,1] = gravitational_dx2dy2
        potential_derivatives[1,3,3] = gravitational_dy2dz2
        potential_derivatives[3,1,3] = gravitational_dx2dz2
        potential_derivatives[3,2,2] = gravitational_dx2dydz
        potential_derivatives[2,3,2] = gravitational_dxdy2dz
        potential_derivatives[2,2,3] = gravitational_dxdydz2
    elseif dims == 2
        # potential
        potential_derivatives[1,1] = gravitational_potential_2D
        # first derivatives
        potential_derivatives[2,1] = gravitational_2d_ddx
        potential_derivatives[1,2] = gravitational_2d_ddy
        # second derivatives
        potential_derivatives[3,1] = gravitational_2d_d2dx2
        potential_derivatives[2,2] = gravitational_2d_dxdy
        potential_derivatives[1,3] = gravitational_2d_d2dy2
        # third derivatives
        potential_derivatives[4,1] = gravitational_2d_d3dx3
        potential_derivatives[3,2] = gravitational_2d_d3dx2dy
        potential_derivatives[2,3] = gravitational_2d_d3dxdy2
        potential_derivatives[1,4] = gravitational_2d_d3dy3
        # fourth derivatives
        potential_derivatives[5,1] = gravitational_2d_d4dx4
        potential_derivatives[4,2] = gravitational_2d_d4dx3dy
        potential_derivatives[3,3] = gravitational_2d_d4dx2dy2
        potential_derivatives[2,4] = gravitational_2d_d4dxdy3
        potential_derivatives[1,5] = gravitational_2d_d4dy4
        # fifth derivatives
        potential_derivatives[6,1] = gravitational_2d_d5dx5
        potential_derivatives[5,2] = gravitational_2d_d5dx4dy
        potential_derivatives[4,3] = gravitational_2d_d5dx3dy2
        potential_derivatives[3,4] = gravitational_2d_d5dx2dy3
        potential_derivatives[2,5] = gravitational_2d_d5dxdy4
        potential_derivatives[1,6] = gravitational_2d_d5dy5
    end

    return Gravitational(G, potential_derivatives, UInt8(order))
end

function gravitational_potential_3D(Rho, m_source, G)
    rho_sq = Rho' * Rho
    rho = sqrt(rho_sq)
    return -G * m_source / rho
end

function gravitational_dx(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dVdx = G * m_source * Rho[1] / rho_sq^1.5
end

function gravitational_dy(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dVdx = G * m_source * Rho[2] / rho_sq^1.5
end

function gravitational_dz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dVdx = G * m_source * Rho[3] / rho_sq^1.5
end

function gravitational_dx2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d2Vdx2 = G * m_source * (-2 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^2.5
end

function gravitational_dy2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d2Vdy2 = G * m_source * (Rho[1]^2 - 2*Rho[2]^2 + Rho[3]^2) / rho_sq^2.5
end

function gravitational_dz2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d2Vdz2 = G * m_source * (Rho[1]^2 + Rho[2]^2 - 2*Rho[3]^2) / rho_sq^2.5
end

function gravitational_dxdy(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d2Vdxdy = -3 * G * m_source * Rho[1] * Rho[2] / rho_sq^2.5
end

function gravitational_dydz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d2Vdydz = -3 * G * m_source * Rho[2] * Rho[3] / rho_sq^2.5
end

function gravitational_dxdz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d2Vdxdz = -3 * G * m_source * Rho[3] * Rho[1] / rho_sq^2.5
end

function gravitational_dx3(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdx3 = 3 * G * m_source * Rho[1] * (2 * Rho[1]^2 - 3 * (Rho[2]^2 + Rho[3]^2)) / rho_sq^3.5
end

function gravitational_dy3(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdy3 = 3 * G * m_source * Rho[2] * (2 * Rho[2]^2 - 3 * (Rho[3]^2 + Rho[1]^2)) / rho_sq^3.5
end

function gravitational_dz3(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdz3 = 3 * G * m_source * Rho[3] * (2 * Rho[3]^2 - 3 * (Rho[1]^2 + Rho[2]^2)) / rho_sq^3.5
end

function gravitational_dx2dy(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdx2dy = -3 * G * m_source * Rho[2] * (-4 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^3.5
end

function gravitational_dx2dz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdx2dy = -3 * G * m_source * Rho[3] * (-4 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^3.5
end

function gravitational_dy2dz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdy2dz = -3 * G * m_source * Rho[3] * (-4 * Rho[2]^2 + Rho[3]^2 + Rho[1]^2) / rho_sq^3.5
end

function gravitational_dxdy2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdy2dz = -3 * G * m_source * Rho[1] * (-4 * Rho[2]^2 + Rho[3]^2 + Rho[1]^2) / rho_sq^3.5
end

function gravitational_dxdz2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdz2dx = -3 * G * m_source * Rho[1] * (-4 * Rho[3]^2 + Rho[1]^2 + Rho[2]^2) / rho_sq^3.5
end

function gravitational_dydz2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdz2dx = -3 * G * m_source * Rho[2] * (-4 * Rho[3]^2 + Rho[1]^2 + Rho[2]^2) / rho_sq^3.5
end

function gravitational_dxdydz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    d3Vdxdydz = 15 * G * m_source * Rho[1] * Rho[2] * Rho[3] / rho_sq^3.5
end

function gravitational_dx4(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dx4 = -3 * G * m_source * (8 * Rho[1]^4 - 24 * Rho[1]^2 * (Rho[2]^2 + Rho[3]^2) + 3 * (Rho[2]^2 + Rho[3]^2)^2) / rho_sq^4.5
end

function gravitational_dy4(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dy4 = -3 * G * m_source * (8 * Rho[2]^4 - 24 * Rho[2]^2 * (Rho[3]^2 + Rho[1]^2) + 3 * (Rho[3]^2 + Rho[1]^2)^2) / rho_sq^4.5
end

function gravitational_dz4(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dz4 = -3 * G * m_source * (8 * Rho[3]^4 - 24 * Rho[3]^2 * (Rho[1]^2 + Rho[2]^2) + 3 * (Rho[1]^2 + Rho[2]^2)^2) / rho_sq^4.5
end

function gravitational_dx3dy(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dx3dy = 15 * G * m_source * Rho[2] * (3 * Rho[1] * (Rho[2]^2 + Rho[3]^2) - 4 * Rho[1]^3) / rho_sq^4.5
end

function gravitational_dx3dz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dx3dz = 15 * G * m_source * Rho[3] * (3 * Rho[1] * (Rho[3]^2 + Rho[2]^2) - 4 * Rho[1]^3) / rho_sq^4.5
end

function gravitational_dy3dz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dy3dz = 15 * G * m_source * Rho[3] * (3 * Rho[2] * (Rho[3]^2 + Rho[1]^2) - 4 * Rho[2]^3) / rho_sq^4.5
end

function gravitational_dxdy3(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dxdy3 = 15 * G * m_source * Rho[1] * (3 * Rho[2] * (Rho[1]^2 + Rho[3]^2) - 4 * Rho[2]^3) / rho_sq^4.5
end

function gravitational_dxdz3(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dxdz3 = 15 * G * m_source * Rho[1] * (3 * Rho[3] * (Rho[1]^2 + Rho[2]^2) - 4 * Rho[3]^3) / rho_sq^4.5
end

function gravitational_dydz3(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dydz3 = 15 * G * m_source * Rho[2] * (3 * Rho[3] * (Rho[2]^2 + Rho[1]^2) - 4 * Rho[3]^3) / rho_sq^4.5
end

function gravitational_dx2dy2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dx2dy2 = 3 * G * m_source * (4 * Rho[1]^4 + 3 * Rho[1]^2 * (Rho[3]^2 - 9 * Rho[2]^2) +
             4 * Rho[2]^4 + 3 * Rho[2]^2 * Rho[3]^2 - Rho[3]^4) / rho_sq^4.5
end

function gravitational_dy2dz2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dy2dz2 = 3 * G * m_source * (4 * Rho[3]^4 + 3 * Rho[3]^2 * (Rho[1]^2 - 9 * Rho[2]^2) +
            4 * Rho[2]^4 + 3 * Rho[2]^2 * Rho[1]^2 - Rho[1]^4) / rho_sq^4.5
end

function gravitational_dx2dz2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dx2dz2 = 3 * G * m_source * (4 * Rho[1]^4 + 3 * Rho[1]^2 * (Rho[2]^2 - 9 * Rho[3]^2) +
            4 * Rho[3]^4 + 3 * Rho[3]^2 * Rho[2]^2 - Rho[2]^4) / rho_sq^4.5
end

function gravitational_dx2dydz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dx2dydz = 15 * G * m_source * Rho[2] * Rho[3] * (-6 * Rho[1]^2 + Rho[2]^2 + Rho[3]^2) / rho_sq^4.5
end

function gravitational_dxdy2dz(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dxdy2dz = 15 * G * m_source * Rho[1] * Rho[3] * (-6 * Rho[2]^2 + Rho[1]^2 + Rho[3]^2) / rho_sq^4.5
end

function gravitational_dxdydz2(Rho, m_source, G)
    rho_sq = Rho' * Rho
    dxdydz2 = 15 * G * m_source * Rho[2] * Rho[1] * (-6 * Rho[3]^2 + Rho[2]^2 + Rho[1]^2) / rho_sq^4.5
end

##
## 2 dimensions
##
function gravitational_potential_2D(Rho, m_source, G)
    rho_sq = Rho' * Rho
    rho = sqrt(rho_sq)
    return G * m_source * log(rho)
end

function gravitational_2d_ddx(Rho, m_source, G)
    x, y = Rho
    return x / (x^2 + y^2) * m_source * G
end

function gravitational_2d_ddy(Rho, m_source, G)
    x, y = Rho
    return y / (x^2 + y^2) * m_source * G
end

function gravitational_2d_d2dx2(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return (y2 - x2) / (x2 + y2)^2 * m_source * G
end

function gravitational_2d_dxdy(Rho, m_source, G)
    x, y = Rho
    return -2 * x * y / (x^2 + y^2)^2 * m_source * G
end

function gravitational_2d_d2dy2(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return (x2 - y2) / (x2 + y2)^2 * m_source * G
end

function gravitational_2d_d3dx3(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 2 * x * (x2 - 3*y2) / (x2 + y2)^3 * m_source * G
end

function gravitational_2d_d3dx2dy(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return -2 * y * (y2 - 3*x2) / (x2 + y2)^3 * m_source * G
end

function gravitational_2d_d3dxdy2(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return -2 * x * (x2 - 3*y2) / (x2 + y2)^3 * m_source * G
end

function gravitational_2d_d3dy3(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    x2_p_y2 = x2 + y2
    return 2 * y * (y^2 - 3 * x^2) / (x2_p_y2)^3 * m_source * G
end

function gravitational_2d_d4dx4(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    x2_p_y2 = x2 + y2
    return -6 * (x^4 - 6*x2 * y2 + y^4) / x2_p_y2^4 * m_source * G
end

function gravitational_2d_d4dx3dy(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return -24 * x * y * (x2 - y2) / (x2 + y2)^4 * m_source * G
end

function gravitational_2d_d4dx2dy2(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 6 * (x^4 - 6 * x2 * y2 + y^4) / (x2 + y2)^4 * m_source * G
end

function gravitational_2d_d4dxdy3(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 24 * x * y * (x2 - y2) / (x2 + y2)^4
end

function gravitational_2d_d4dy4(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    x2_p_y2 = x2 + y2
    return -6 * (x^4 - 6*x2 * y2 + y^4) / x2_p_y2^4 * m_source * G
end

function gravitational_2d_d5dx5(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 24 * x * (x^4 - 10*x2*y2 + 5*y^4) / (x2 + y2)^5
end

function gravitational_2d_d5dx4dy(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 24 * y * (5 * x^4 - 10 * x2 * y2 + y^4) / (x2 + y2)^5
end

function gravitational_2d_d5dx3dy2(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return -24 * x * (x^4 - 10 * x2 * y2 + 5 * y^4) / (x2 + y2)^5
end

function gravitational_2d_d5dx2dy3(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return -24 * y * (5 * x^4 - 10 * x2 * y2 + y^4) / (x2 + y2)^5
end

function gravitational_2d_d5dxdy4(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 24 * x * (x^4 - 10 * x2 * y2 + 5 * y^4) / (x2 + y2)^5
end

function gravitational_2d_d5dy5(Rho, m_source, G)
    x, y = Rho
    x2 = x^2
    y2 = y^2
    return 24 * y * (5 * x^4 - 10 * x2 * y2 + y^4) / (x2 + y2)^5
end
