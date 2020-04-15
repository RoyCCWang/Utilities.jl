
import Pkg

#println("pre instantiate")
Pkg.instantiate()
#println("post instantiate")

using Test
#using LinearAlgebra
using Utilities

tests = ["toy_check"]

for t in tests
    @testset "$t" begin
        include("$(t).jl")
    end
end
