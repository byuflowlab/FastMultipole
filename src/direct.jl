function direct(sources::Vector{e}, x_target::Vector{TF}, kernel) where {e<:Element, TF}
    V = 0.0
    for source in sources
        x_source = get_X(source)
        q_source = get_q(source)
        Rho = x_target - x_source
        rho_squared = Rho' * Rho
        if rho_squared > 0.0
            V += kernel(x_source, q_source, x_target)
        end
    end
    return V
end

function direct(sources::Vector{e1}, targets::Vector{e2}, kernel) where {e1<:Element, e2<:Element}
    V = zeros(length(targets))
    for (i,target) in enumerate(targets)
        x_target = get_X(target)
        V[i] = direct(sources, x_target, kernel)
    end
    return V
end

function direct(sources::Vector{e}, kernel) where e <: Element
    return direct(sources, sources, kernel)
end
