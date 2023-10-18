# overload for SortWrapper
Base.length(sys::SortWrapper) = length(sys.system)

Base.setindex!(sys::SortWrapper,val,i) = setindex!(sys.system,val,sys.index[i])

Base.setindex!(sys::SortWrapper,val,i,parameter::ScalarPotential) = setindex!(sys.system,val,sys.index[i],parameter::ScalarPotential)

Base.setindex!(sys::SortWrapper,val,i,parameter::VectorPotential) = setindex!(sys.system,val,sys.index[i],parameter::VectorPotential)

Base.setindex!(sys::SortWrapper,val,i,parameter::Velocity) = setindex!(sys.system,val,sys.index[i],parameter::Velocity)

Base.setindex!(sys::SortWrapper,val,i,parameter::VelocityGradient) = setindex!(sys.system,val,sys.index[i],parameter::VelocityGradient)

Base.getindex(sys::SortWrapper, i) = getindex(sys.system, sys.index[i])

Base.getindex(sys::SortWrapper, i, parameter) = getindex(sys.system, sys.index[i], parameter)

B2M!(system::SortWrapper, branch, bodies_index, harmonics, expansion_order) =
    B2M!(system.system, branch, system.index[bodies_index], harmonics, expansion_order)

direct!(target_system, target_index, source_system::SortWrapper, source_index) =
    direct!(target_system, target_index, source_system.system, source_system.index[source_index])