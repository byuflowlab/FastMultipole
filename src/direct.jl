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

function direct!(sources::Vector{e1}, targets::Vector{e2}, kernel) where {e1<:Element, e2<:Element}
    for (i,target) in enumerate(targets)
        x_target = get_X(target)
        V = direct(sources, x_target, kernel)
        add_V(target, V)
    end
    return nothing
end

function direct(sources::Vector{e1}, targets::Vector{e2}, kernel) where {e1<:Element, e2<:Element}
    Vs = zeros(length(targets))
    for (i,target) in enumerate(targets)
        x_target = get_X(target)
        Vs[i] = direct(sources, x_target, kernel)
    end
    return Vs
end

function direct!(sources::Vector{e}, kernel) where e <: Element
    direct!(sources, sources, kernel)
end

function direct(sources::Vector{e}, kernel) where e <: Element
    return direct(sources, sources, kernel)
end

function direct(source::e, target::Vector{TF}, kernel) where {e <: Element, TF}
    x_source = get_X(source)
    q_source = get_q(source)
    Rho = x_target - x_source
    rho_squared = Rho' * Rho
    V = rho_squared > 0.0 ? kernel(x_source, q_source, x_target) : 0.0
end

function direct!(source::e1, target::e2, kernel) where {e1 <: Element, e2 <: Element}
    x_source = get_X(source)
    q_source = get_q(source)
    x_target = get_X(target)
    Rho = x_target - x_source
    rho_squared = Rho' * Rho
    V = rho_squared > 0.0 ? kernel(x_source, q_source, x_target) : 0.0
    add_V(target, V)
    return nothing
end
