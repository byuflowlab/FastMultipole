import FLOWVPM
vpm = FLOWVPM
import Elliptic
import Roots
import Cubature

include(joinpath(dirname(pathof(vpm)), "..", "examples", "vortexrings", "vortexrings_functions.jl"))

pfield = vpm.ParticleField(300)
circulation = 1.0/10
R = 1.0
AR = 2.5
Rcross = 0.4
Nphi = 10
nc = 1
sigma = 0.5
addvortexring(pfield, circulation, R, AR, Rcross, Nphi, nc, sigma;
    Oaxis=[1.0 0 0;0 1.0 0;0 0 1.0])

# try original setup
vpm.UJ_direct(pfield)
u_direct = zeros(3,pfield.np)
for i in 1:pfield.np
    u_direct[:,i] .= pfield.particles[i].U
end

new_pfield = vpm.ParticleField(300)
addvortexring(new_pfield, circulation, R, AR, Rcross, Nphi, nc, sigma;
    Oaxis=[1.0 0 0;0 1.0 0;0 0 1.0])

# try my setup
vpm.UJ_fmm(new_pfield, vpm.Fmm())
u_myfmm = zeros(3,new_pfield.np)
for i in 1:new_pfield.np
    u_myfmm[:,i] .= new_pfield.particles[i].U
end

@show u_direct-u_myfmm