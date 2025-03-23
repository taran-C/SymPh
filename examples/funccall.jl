using SymbolicPhysics.Arrays

f() = print("test")

a = ArrayVariable("a")
y = 2*a
y.name = "y"
c = 1+y
c.name = "ci"
z = 3+c[1,0]
z.name = "z"
x = FuncCall("x", "f", [z, y], 0, 0) #Find a way to do this in a prettier way than a fname string
b = 2*x[2,0]
b.name = "b"

t = to_deptree!(Set{String}(), b)
println(Arrays.to_graphviz(t))

s = to_sequence!(t)
println(string(s))

func!, str = to_kernel(s)
println(str)
