using SymPh.Maths

x = FormVariable{0, Primal}("x")
y = FormVariable{0, Primal}("y")

z = x + y
println(z)
println(z isa Form{0, Primal})

dx = ExteriorDerivative(x)
println(dx)
println(dx isa Form{1, Primal})

dz = ExteriorDerivative(z)
println(dz)
println(dz isa Form{1, Primal})

v = VectorVariable{Primal}("v")
ivdx = InteriorProduct(v,dx)
println(ivdx)
println(ivdx.name)
println(ivdx isa Form{0, Primal})
