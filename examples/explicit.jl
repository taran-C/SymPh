using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

for i in 1:1
	#Defining our equation
	a = FormVariable{0,Primal}("a")
	da = ExteriorDerivative(a)
	dda = ExteriorDerivative(da)

	#-----This part will be encapsulated into a function that automatically optimizes ----------
	#Transforming the Forms expression into an Expression on arrays
	exprs = explicit(dda)
	println("Developped expression :")
	println(string(exprs))

	#Transforming our Expression into a dependency tree
	tree = Arrays.to_deptree!(Set{String}(["da_x", "da_y"]), exprs)
	println("Tree view")
	println(string(tree))
	println("Graphviz view of tree")
	println(Arrays.to_graphviz(tree))
	
	#Transforming our dependency tree into a sequence of expressions to compute
	seq = Arrays.to_sequence!(tree)
	println("Corresponding Sequence :")
	println(string(seq)*"\n")

	#Generating the final function
	func!, funcstr = Arrays.to_kernel(seq)
	println("Generated code :")
	println(funcstr)
	#-------------------------------------------------------------------------------------------

	#Testing the function
	nx = 10
	ny = 10
	nh = 2

	A = zeros(nx, ny)
	for i in 1:nx, j in 1:ny
		A[i,j] = i+j
	end

	DAx = zeros(nx, ny)
	DAy = zeros(nx, ny)

	DDA = zeros(nx, ny)

	func!(;nx=nx, ny=ny, nh=nh, a=A, da_x=DAx, da_y = DAy, dda = DDA)

	display(DAx)
	display(DAy)
	display(DDA)
end
