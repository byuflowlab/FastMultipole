# overload for SortWrapper
Base.length(sys::SortWrapper) = length(sys.system)

Base.setindex!(sys::SortWrapper,val,i) = setindex!(sys.system,val,sys.index[i])
Base.setindex!(sys::SortWrapper,val,i,parameter::ScalarPotential) = setindex!(sys.system,val,sys.index[i],parameter)
Base.setindex!(sys::SortWrapper,val,i,parameter::VectorPotential) = setindex!(sys.system,val,sys.index[i],parameter)
Base.setindex!(sys::SortWrapper,val,i,parameter::Velocity) = setindex!(sys.system,val,sys.index[i],parameter)
Base.setindex!(sys::SortWrapper,val,i,parameter::VelocityGradient) = setindex!(sys.system,val,sys.index[i],parameter)

Base.getindex(sys::SortWrapper, i) = getindex(sys.system, sys.index[i])
Base.getindex(sys::SortWrapper, i, parameter::Position) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Radius) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::ScalarStrength) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::VectorStrength) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::ScalarPotential) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::VectorPotential) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Velocity) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::VelocityGradient) = getindex(sys.system, sys.index[i], parameter)
Base.getindex(sys::SortWrapper, i, parameter::Vertex, i_vertex) = Base.getindex(sys.system, sys.index[i], parameter, i_vertex)

Base.eltype(sys::SortWrapper) = Base.eltype(sys.system)

B2M!(system::SortWrapper, branch, bodies_index, harmonics, expansion_order) =
    B2M!(system.system, branch, system.index[bodies_index], harmonics, expansion_order)

direct!(target_system, target_index, source_system::SortWrapper, source_index) =
    direct!(target_system, target_index, source_system.system, source_system.index[source_index])