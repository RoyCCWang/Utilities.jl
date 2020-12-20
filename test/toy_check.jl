
x = randn()

LHS = sin(x)
RHS = cos(x-pi/2)

@test isnumericallyclose(LHS, RHS)

println("LHS = ", LHS)
println("RHS = ", RHS)
