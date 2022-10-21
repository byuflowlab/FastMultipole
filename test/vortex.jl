import FLOWFMM
fmm = FLOWFMM

#####
##### classic vortex particle method
#####

struct Vorton <: fmm.VectorSource
    position
    strength
    velocity
    potential
    J_potential # jacobian of the vector potential
    H_potential # hessian of the vector potential
end

function Vorton(position, strength;
    velocity = zeros(3),
    potential = zeros(3),
    J_potential = zeros(3,3),
    H_potential = zeros(3,3,3)
)
    Vorton(position, strength, velocity, J_potential, H_potential)
end

