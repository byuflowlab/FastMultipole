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

test_dir = @__DIR__

include(joinpath(test_dir, "gravitational.jl"))


@testset "Update Radius" begin
    xs = [
        1.2 1.1 0.8;
        0.8 0.9 0.2;
        0.1 0.2 0.9;
        0.1 0.3 0.2;
        0.2 0.25 0.4
    ]
    
    radii = [
        0.8,
        1.1,
        2.2,
        0.5,
        1.9
    ]

    ms = [
        0.8,
        1.1,
        2.2,
        0.5,
        1.9
    ]


    correct_bodies = [
        .1 .2 .8 2.31533193;
        .2 .25 .4 1.9;
        .8 .9 .2 1.1;
        .1 .2 .9 2.2;
        1.2 1.1 .8 .8;
        .1 .3 .2 .5;
        .2 .25 .4 1.9;
    ]
    
    bodies = zeros(8,length(radii))
    
    for i in 1:length(radii)
        bodies[1:3,i] .= xs[i,1:3]
        bodies[4,i] = radii[i]
        bodies[5,i] = ms[i]
    end
    elements = Gravitational(bodies) # change to volumetric 
    
    expansion_order = 24
    theta = 4
    options = fmm.Options(expansion_order, 1, 4.0,true,true)
    tree = fmm.Tree((elements,), options)

    # Test step_through_bodies for correct values
    for j = 1:length(correct_bodies[:,1])
        for i = 1:3
            # println(string(round(tree.branches[j].center[i],digits=5)) * "    " * string(round(correct_bodies[j,i],digits=5)))
            @test round(tree.branches[j].center[i],digits=5) == round(correct_bodies[j,i],digits=5)
            # @show round(tree.branches[j].center[i],digits=5) round(correct_bodies[j,i],digits=5)
        end
        # println(string(round(tree.branches[j].radius,digits=5)) * "     " * string(round(correct_bodies[j,4],digits=5)))
        # println()

        @test round(tree.branches[j].radius[1],digits=5) == round(correct_bodies[j,4],digits=5)
    end
end


