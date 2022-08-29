#####
##### KERNELS
#####
abstract type Kernel{TF,dims} end

(K::Kernel{TF,dims})(n, x_source, q_source, x_target) where {TF,dims} =
    K.potential_derivatives[n...](x_target - x_source) * q_source

(K::Kernel{TF,dims})(x_source, q_source, x_target) where {TF,dims} =
    K.potential_derivatives[1](x_target - x_source) * q_source

function kernel end

function kernel! end
