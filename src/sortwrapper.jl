"""
    SortWrapper(system)

Convenience wrapper for systems whose elements cannot be sorted in-place (e.g. structured grids). The resulting object is treated like any other `system`.
"""
function SortWrapper(system)
    return SortWrapper(system,collect(1:get_n_bodies(system)))
end

@inline wrap_duplicates(target_systems::Tuple, source_systems::Tuple) = Tuple(target_system in source_systems ? SortWrapper(target_system) : target_system for target_system in target_systems)

@inline wrap_duplicates(target_system, source_system) = target_system == source_system ? SortWrapper(target_system) : target_system

@inline wrap_duplicates(target_system, source_systems::Tuple) = target_system in source_systems ? SortWrapper(target_system) : target_system

@inline wrap_duplicates(target_systems::Tuple, source_system) = Tuple(target_system == source_system ? SortWrapper(target_system) : target_system for target_system in target_systems)


# access functions for SortWrapper
get_n_bodies(sys::SortWrapper) = get_n_bodies(sys.system)

Base.setindex!(sys::SortWrapper,val,i) = setindex!(sys.system,val,sys.index[i])
Base.setindex!(sys::SortWrapper,val,i,parameter::ScalarPotential) = setindex!(sys.system,val,sys.index[i],parameter)
#Base.setindex!(sys::SortWrapper,val,i,parameter::VectorPotential) = setindex!(sys.system,val,sys.index[i],parameter)
Base.setindex!(sys::SortWrapper,val,i,parameter::Velocity) = setindex!(sys.system,val,sys.index[i],parameter)
Base.setindex!(sys::SortWrapper,val,i,parameter::VelocityGradient) = setindex!(sys.system,val,sys.index[i],parameter)

Base.getindex(sys::SortWrapper, i) = getindex(sys.system, sys.index[i])
Base.getindex(sys::SortWrapper, i, parameter::Position) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Radius) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Strength) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::ScalarPotential) = getindex(sys.system, sys.index[i], parameter)
#Base.getindex(sys::SortWrapper, i, parameter::VectorPotential) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Velocity) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::VelocityGradient) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Vertex, i_vertex) = Base.getindex(sys.system, sys.index[i], parameter, i_vertex)

Base.eltype(sys::SortWrapper) = Base.eltype(sys.system)

body_to_multipole!(system::SortWrapper, branch, bodies_index, harmonics, expansion_order) =
    body_to_multipole!(system.system, branch, system.index[bodies_index], harmonics, expansion_order)
