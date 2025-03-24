using SymbolicPhysics.Arrays
import SymbolicPhysics: @Let, State

f(out, a, b) = print("test")

@Let a = ArrayVariable()
@Let y = 2*a
@Let c = 1+y
@Let z = 3+c[1,0]
@Let x = FuncCall("f", [y, c]) #Find a way to do this in a prettier way than a fname string
@Let b = 2*x[2,0]

t = to_deptree!(Set{String}(), b)
println(Arrays.to_graphviz(t))

s = to_sequence!(t)
println(string(s))

func!, str = to_kernel(s)
println(str)

mesh = Arrays.Mesh(15, 15, 3, ones(15,15), 1, 1)
state = State(mesh)

func!(mesh; c = state.c, x = state.x, b = state.b, a = state.a, y = state.y)
