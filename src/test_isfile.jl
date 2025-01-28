if isfile("error.jl")
    println("Found it!")
else
    println("Nope.")
end

@show readdir(".")
