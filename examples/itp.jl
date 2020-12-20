
import PyPlot
import Random
using LinearAlgebra

import StaticArrays
import Interpolations

import Calculus

include("../src/function_interpolation.jl")
include("../src/patterns.jl")
include("../src/visualize_2D.jl")

PyPlot.close("all")
fig_num = 1

Random.seed!(25)

D = 2

N_array = [50; 50]

x_a = [-3.5; 2.5]
x_b = [3.5; 4.5]


x_ranges = collect( LinRange(x_a[d], x_b[d], N_array[d]) for d = 1:D )
X_nD = ranges2collection(x_ranges, Val(2))

f = xx->sinc.(3*norm(xx))
A = f.(X_nD)

f_itp, d_f_itp, d2_f_itp = setupcubicitp(A, x_ranges, 1.0)
#f_itp, d_f_itp!, d2_f_itp! = setupcubicitpmutator(A, x_ranges, 1.0)

fig_num = visualizemeshgridpcolor(x_ranges, f_itp.(X_nD),
                vec(X_nD), ".",
               fig_num, "f_itp", "x1", "x2")

df_ND = xx->Calculus.gradient(f, xx)

println("Lower bounds: ", x_a)
println("Upper bounds: ", x_b)
println()

x0 = randn(D) + (x_a+x_b) ./2
df_ND(x0)
println("In bounds:")
println("x0 = ", x0)
#println("d_f_itp!(x0) = ", d_f_itp!(x0))
println("d_f_itp(x0) = ", d_f_itp(x0))
println("df_ND(x0)   = ", df_ND(x0))
println()

x0 = randn(D)
df_ND(x0)
println("Out of bounds:")
println("x0 = ", x0)
println("d_f_itp(x0) = ", d_f_itp(x0))
println("df_ND(x0)   = ", df_ND(x0))
println()
