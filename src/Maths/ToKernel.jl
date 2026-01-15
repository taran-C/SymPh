export to_kernel

"""
	to_kernel(exprs...; explparams = ExplicitParam(), verbose=false, bcs = [])

Takes one or multiple expressions and returns a kernel that computes them.
The kernel has signature `compute!(mesh::Mesh, state::State)`.

# Keyword Arguments
- `explparams::ExplicitParam` : An object holding the parameters like which interpolation to use, etc... (see [`ExplicitParam`](@ref))
- `verbose::Bool` : Wether or not to print the functions being generated, the dependency tree, etc...
- `bcs` : A list of objects representing our boundary conditions (WIP)
"""
function to_kernel(exprs...; explparams = ExplicitParam(), verbose=0, bcs = [])
	#Transforming the Forms expression into an Expression on arrays (TODO probably a more elegant way to do this)
	math_exprs = []
	for expr in exprs
		expl_expr = explicit(expr; param = explparams)
		if expl_expr isa AbstractArray
			math_exprs = vcat(math_exprs, expl_expr)
		else
			push!(math_exprs, expl_expr)
		end
	end
	
	if verbose >= 1
		println("Developped expression :")
		println(string(exprs))
	end

	if verbose >= 2
		println("Explicit expression :")
		println(string(math_exprs))
	end

	#Transforming our Expression into a dependency tree
	tree = Arrays.to_deptree!(math_exprs)
	
	if verbose >= 1
		println("Tree view :")
		println(string(tree))

		println("Graphviz view of tree:")
		println(Arrays.to_graphviz(tree))
	end


	#Transforming our dependency tree into a sequence of expressions to compute
	seq = Arrays.to_sequence!(tree)

	if verbose >= 1
		println("Corresponding Sequence :")
		println(string(seq)*"\n")
	end

	#Handling the BCs
	fill = []
	for bc in bcs
		if bc isa Vect
			push!(fill, bc.name * "_X")
			push!(fill, bc.name * "_Y")
		elseif degree(bc) in (0, 2)
			push!(fill, bc.name)
		elseif degree(bc) == 1
			push!(fill, bc.name * "_i")
			push!(fill, bc.name * "_j")
		end
	end
	
	#Generating the final function
	if verbose >= 2
		println("Generated functions :\n")
	end
	func!,  vars = Arrays.to_kernel(seq, fill; verbose = verbose)

	return func!
end
