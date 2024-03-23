"""
dr_k/dx_idx_j
"""
function d2rdx2(r, multipole_acceptance_criterion, phi)
    derivatives = zeros(3,3,3)
    derivatives[:,:,1] .= [
        (1-cos(phi)^2 * sin(multipole_acceptance_criterion)^2)/r -sin(multipole_acceptance_criterion)^2*cos(phi)*sin(phi)/r -sin(multipole_acceptance_criterion)*cos(phi)*cos(multipole_acceptance_criterion)/r;
        (-sin(multipole_acceptance_criterion)^2*cos(phi)*sin(phi))/r (1-sin(multipole_acceptance_criterion)^2*sin(phi)^2)/r -sin(multipole_acceptance_criterion)*sin(phi)*cos(multipole_acceptance_criterion)/r;
        -sin(multipole_acceptance_criterion)*cos(phi)*cos(multipole_acceptance_criterion)/r -sin(multipole_acceptance_criterion)*sin(phi)*cos(multipole_acceptance_criterion)/r sin(multipole_acceptance_criterion)^2/r
    ]
    derivatives[:,:,2] .= [
        cos(multipole_acceptance_criterion)/sin(multipole_acceptance_criterion)*(1-cos(phi)^2*(1+2*sin(multipole_acceptance_criterion)^2))/r^2 -cos(multipole_acceptance_criterion)/sin(multipole_acceptance_criterion)*sin(phi)*cos(phi)*(1+2*sin(multipole_acceptance_criterion)^2)/r^2 cos(phi)*(1-2*cos(multipole_acceptance_criterion)^2)/r^2;
        -cos(multipole_acceptance_criterion)/sin(multipole_acceptance_criterion)*sin(phi)*cos(phi)*(1+2*sin(multipole_acceptance_criterion)^2)/r^2 cos(multipole_acceptance_criterion)/sin(multipole_acceptance_criterion)*(1-sin(phi)^2*(1+2*sin(multipole_acceptance_criterion)^2))/r^2 (2*sin(multipole_acceptance_criterion)^2-1)/r^2*sin(phi);
        cos(phi)*(1-2*cos(multipole_acceptance_criterion)^2)/r^2 (2*sin(multipole_acceptance_criterion)^2-1)/r^2*sin(phi) 2*sin(multipole_acceptance_criterion)*cos(multipole_acceptance_criterion)/r^2
    ]
    derivatives[:,:,3] .= [
        2*cos(phi)*sin(phi)/r^2/sin(multipole_acceptance_criterion)^2 (2*sin(phi)^2-1)/r^2/sin(multipole_acceptance_criterion)^2 0;
        (2*sin(phi)^2-1)/r^2/sin(multipole_acceptance_criterion)^2 -2*sin(phi)*cos(phi)/r^2/sin(multipole_acceptance_criterion)^2 0;
        0 0 0
    ]
    return derivatives
end

"""
dr_j/dx_i
"""
function drdx(r,multipole_acceptance_criterion,phi)
    derivatives = [
        sin(multipole_acceptance_criterion)*cos(phi) cos(multipole_acceptance_criterion)*cos(phi)/r -sin(phi)/r/sin(multipole_acceptance_criterion);
        sin(multipole_acceptance_criterion)*sin(phi) cos(multipole_acceptance_criterion)*sin(phi)/r cos(phi)/r/sin(multipole_acceptance_criterion);
        cos(multipole_acceptance_criterion) -sin(multipole_acceptance_criterion)/r 0
    ]
    # derivatives = [
    #     sin(multipole_acceptance_criterion)*cos(phi) sin(multipole_acceptance_criterion)*sin(phi) cos(multipole_acceptance_criterion);
    #     cos(multipole_acceptance_criterion)*cos(phi)/r cos(multipole_acceptance_criterion)*sin(phi)/r -sin(multipole_acceptance_criterion)/r;
    #     -sin(phi)/r/sin(multipole_acceptance_criterion) cos(phi)/r/sin(multipole_acceptance_criterion) 0
    # ]
    return derivatives
end