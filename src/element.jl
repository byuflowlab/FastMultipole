#####
##### abstract Element definition
#####
abstract type Element{TF} end

abstract type VectorPotential{TF} <: Element{TF} end

abstract type ScalarPotential{TF} <: Element{TF} end
