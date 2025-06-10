export to_kernel

function to_kernel(exprs...; save = [], explparams = ExplicitParam(), verbose=false, bcs = [])
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
	
	if verbose
		println("Developped expression :")
		println(string(exprs))
	end

	if verbose
		println("Explicit expression :")
		println(string(math_exprs))
	end

	#Transforming our Expression into a dependency tree
	tree = Arrays.to_deptree!(Set{String}(save), math_exprs)
	
	if verbose
		println("Tree view :")
		println(string(tree))

		println("Graphviz view of tree:")
		println(Arrays.to_graphviz(tree))
	end


	#Transforming our dependency tree into a sequence of expressions to compute
	seq = Arrays.to_sequence!(tree)

	if verbose
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
	if verbose
		println("Generated functions :\n")
	end
	func!,  vars = Arrays.to_kernel(seq, fill; verbose = verbose)

	return func!
end
