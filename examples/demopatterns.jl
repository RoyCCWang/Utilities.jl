
function demoscaletoChebyinterval()

    # ground truth: we know u for x[i] is i.
    x = LinRange(-1,1,5)
    println(x)

    c = collect(1:5)
    u = coordinatetolattice1D(collect(x),x[1],x[end], length(x))

    x_rec = lattice1Dtocoordinate(u,x[1],x[end], length(x))
    println("u is ", u)
    println("x_rec is ", xx)
end
