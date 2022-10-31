function update_velocity!(elements)
    for element in elements
        # velocity is the curl of the vector potential
        element.velocity[1] = (-element.J_potential[1,1] + element.J_potential[2,4] - element.J_potential[3,3]) * ONE_OVER_4PI
        element.velocity[2] = (-element.J_potential[2,2] + element.J_potential[3,2] - element.J_potential[1,4]) * ONE_OVER_4PI
        element.velocity[3] = (-element.J_potential[3,3] + element.J_potential[1,3] - element.J_potential[2,2]) * ONE_OVER_4PI
    end
end