using SymbolicPhysics.Arrays

#Defining an expression. Could be done in a prettier way without naming everything wich allows for infix operator usage but makes saving each variable harder
x = ArrayVariable("x")

a = Arrays.Multiplication("a", Arrays.RealValue(2), x)
b = Arrays.Multiplication("b", x[1,0],x[-1,0])
d = Arrays.Substraction("d", Arrays.RealValue(3), b[0,1])
e = Arrays.Negative("e", d)
c = Arrays.Addition("c", a, b)
y = Arrays.Addition("y", c, e)
z = Arrays.Multiplication("z", a,c)

#Showing the developed expression, only directly depending on arrays
println("Developed Expression :")
println(string(y)*"\n")

#Converting the expression into a Dependency tree, with at its root an empty node allowing us to fuse multiple expressions at once
tree = to_deptree!(Set{String}(["c"]), y, z)
#println(tree)

#Converting again into a sequence of blocks showing which calculations are to be done in which order
println("Corresponding Sequence :")
seq = to_sequence!(tree)
println(string(seq)*"\n")

#And finally generating a function from this sequence
println("Generated code :")
func!, funcstr = to_kernel(seq)
println(funcstr)

