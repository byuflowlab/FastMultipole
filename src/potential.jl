function update_velocity!(elements)
    for i in 1:size(elements.bodies)[2]
        # velocity is the curl of the vector potential
        # minus the gradient of the scalar potential
        jacobian = reshape(elements.potential[i_POTENTIAL_JACOBIAN,i],3,4)
        @show elements.velocity[:,i]
        elements.velocity[1,i] = (-jacobian[1,1] + jacobian[2,4] - jacobian[3,3]) * ONE_OVER_4PI
        elements.velocity[2,i] = (-jacobian[2,1] + jacobian[3,2] - jacobian[1,4]) * ONE_OVER_4PI
        elements.velocity[3,i] = (-jacobian[3,1] + jacobian[1,3] - jacobian[2,2]) * ONE_OVER_4PI
        @show elements.velocity[:,i]
    end
    @show elements.potential[i_POTENTIAL_JACOBIAN]
end