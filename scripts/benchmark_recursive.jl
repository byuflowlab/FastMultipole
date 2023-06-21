function forloop(a1, n)
    for i in 1:n-1
        a1 = cos(a1) * exp(0.01)
    end
    return a1
end

function recurse(a1,n)
    recurse(a1,1,n)
end

function recurse(a1,i,n)
    if i == n
        return a1
    else
        a1 = cos(a1) * exp(0.01)
        recurse(a1,i+1,n)
    end
end

using BenchmarkTools

vf = @benchmark forloop(1.0,100000)
vr = @benchmark recurse(1.0,100000)
vf = forloop(1.0,100000)
vr = recurse(1.0,100000)

@show vf vr