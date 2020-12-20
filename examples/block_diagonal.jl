
import Random
using LinearAlgebra
import Printf


include("../src/matrix.jl")


Random.seed!(25)

##### block inverse.
D_x = 3
D_y = D_x

N = 4

y_set = collect( randn(D_y, D_x) for n = 1:N )

Y = makeblockdiagonalmatrix(y_set)

inv_y_set = similar(y_set)
invertblockmatrix!(inv_y_set, y_set)
inv_y =  makeblockdiagonalmatrix(inv_y_set)

Y_inv = inv(Y)

println("Inverse:")
println("discrepancy between inv_y and Y_inv is ", norm(Y_inv - inv_y))
println()
