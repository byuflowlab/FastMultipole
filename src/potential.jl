function update_velocity!(velocity, body, potential)
    # velocity is the curl of the vector potential
    # minus the gradient of the scalar potential
    jacobian = reshape(potential[i_POTENTIAL_JACOBIAN],3,4)
    velocity[1] = (-jacobian[1,1] + jacobian[2,4] - jacobian[3,3])
    velocity[2] = (-jacobian[2,1] + jacobian[3,2] - jacobian[1,4])
    velocity[3] = (-jacobian[3,1] + jacobian[1,3] - jacobian[2,2])

    return nothing
end

function update_stretching!(stretching, body, potential)
    # jacobian of the velocity
    hessian = reshape(potential[i_POTENTIAL_HESSIAN],3,3,4)
    duidxj = @SMatrix [
        -hessian[1,1,1]+hessian[2,1,4]-hessian[3,1,3] -hessian[1,2,1]+hessian[2,2,4]-hessian[3,2,3] -hessian[1,3,1]+hessian[2,3,4]-hessian[3,3,3];
        -hessian[2,1,1]+hessian[3,1,2]-hessian[1,1,4] -hessian[2,2,1]+hessian[3,2,2]-hessian[1,2,4] -hessian[2,3,1]+hessian[3,3,2]-hessian[1,3,4];
        -hessian[3,1,1]+hessian[1,1,3]-hessian[2,1,2] -hessian[3,2,1]+hessian[1,2,3]-hessian[2,2,2] -hessian[3,3,1]+hessian[1,3,3]-hessian[2,3,2];
    ]

    # stretching term (omega dot nabla)
    fmm.mul!(stretching, duidxj, body[i_STRENGTH_vortex])

    return nothing
end
