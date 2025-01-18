using SymbolicPhysics.Maths

x = FormVariable{0,Primal}("x")

xexp = explicit(x)
println(xexp)
