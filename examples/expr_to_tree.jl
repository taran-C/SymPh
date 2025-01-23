using SymbolicPhysics.Arrays

x = ArrayVariable("x")

a = Arrays.Multiplication("a", Arrays.RealValue(2), x)
b = Arrays.Multiplication("b", x,x)
d = Arrays.Substraction("d", Arrays.RealValue(3), b)
e = Arrays.Negative("e", d)

c = Arrays.Addition("c", a, b)

y = Arrays.Addition("y", c, e)

println("Developed Expression :")
println(string(y)*"\n")

tree = to_deptree!(y, Set{String}(["c", "d", "b", "a"]))
#println(tree)

println("Corresponding Sequence :")
seq = to_sequence!(tree)
println(string(seq)*"\n")

println("Generated code :")
func!, funcstr = Arrays.to_kernel(seq)
println(funcstr)
