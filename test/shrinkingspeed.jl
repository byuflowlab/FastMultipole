using Test
import Statistics
S = Statistics
import Symbolics
sym = Symbolics
using StaticArrays

import FLOWFMM
fmm = FLOWFMM

using LegendrePolynomials
import PyPlot
plt = PyPlot
using LaTeXStrings
using ProgressMeter

test_dir = @__DIR__

include(joinpath(test_dir, "gravitational.jl"))

panels_runs = [3,3]
panels = [1,10,100,1000,10000,100000,1000000]
metrics = zeros(panels_runs[1],5)

expansion_order = 4
theta = 4
options1 = fmm.Options(expansion_order, 1, 4.0,false,false)
options2 = fmm.Options(expansion_order, 1, 4.0,true,false)
options3 = fmm.Options(expansion_order, 1, 4.0,true,true)

runs = [1,10,100,1000,10000,100000]
control = zeros(runs[panels_runs[2]])
shrinking = zeros(runs[panels_runs[2]])
tree_shrinking = zeros(runs[panels_runs[2]])
tree_shrinking_second = zeros(runs[panels_runs[2]])
second_pass = zeros(runs[panels_runs[2]])

first_body = rand(5, 10)
first_element = Gravitational(first_body) 
first_tree_control= fmm.Tree((first_element,), options1)
first_tree = fmm.Tree((first_element,), options2)
first_tree_second = fmm.Tree((first_element,), options3)
first_shrinking = fmm.update_radius((first_element,),first_tree.branches,1)
first_secondpass = fmm.update_radius((first_element,),first_tree.branches,1;second_pass=false)


println(string(panels_runs[1]) * " Iterations to Run")

for j = 1:panels_runs[1]
    k = panels_runs[1]+1-j
    local bodies = rand(5, panels[k])
    local elements = Gravitational(bodies) 
    local tree = fmm.Tree((elements,), options1)
    @showprogress .001 for i = 1:runs[panels_runs[2]]
        control[i] = @elapsed fmm.Tree((elements,), options1)
        shrinking[i] = @elapsed fmm.update_radius((elements,),tree.branches,1)
        second_pass[i] = @elapsed fmm.update_radius((elements,),tree.branches,1;second_pass=false)
        tree_shrinking[i] = @elapsed fmm.Tree((elements,), options2)
        tree_shrinking_second[i] = @elapsed fmm.Tree((elements,), options3)
   end

    metrics[k,1] = S.mean(control) 
    metrics[k,2] = S.mean(shrinking) 
    metrics[k,3] = S.mean(second_pass) 
    metrics[k,4] = S.mean(tree_shrinking)
    metrics[k,5] = S.mean(tree_shrinking_second)
end

plt.close("all")
plt.loglog(panels[1:panels_runs[1]],metrics[:,:])
plt.title("Relative Speed of Calculations")
plt.legend(["Tree Creation","Shrinking Only","Second Pass","TC with Shrinking","TC with Shrinking/2Pass"])
plt.grid("on")
plt.xlabel("Number of Bodies")
plt.ylabel("Average (1000) Time for Single Calculation")