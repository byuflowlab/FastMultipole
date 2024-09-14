function evaluate_local!(systems, branch::MultiBranch, derivatives_switches, expansion_order)
    for i in eachindex(systems)
        evaluate_local!(systems[i], branch.bodies_index[i], branch.local_expansion, branch.center, derivatives_switches[i], expansion_order)
    end
end

function evaluate_local!(system, branch::SingleBranch, derivatives_switch, expansion_order)
    evaluate_local!(system, branch.bodies_index, branch.local_expansion, branch.center, derivatives_switch, expansion_order)
end

function evaluate_local!(system, bodies_index, local_expansion, expansion_center, derivatives_switch::DerivativesSwitch{PS,VPS,VS,GS}, expansion_order::Val{P}) where {PS,VPS,VS,GS,P}
    for i_body in bodies_index

    end
end
