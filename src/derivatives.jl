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