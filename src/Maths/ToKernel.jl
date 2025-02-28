export to_kernel

function to_kernel(exprs...; verbose=false)
	#Transforming the Forms expression into an Expression on arrays (TODO probably a more elegant way to do this)
	math_exprs = []
	for expr in exprs
		expl_expr = explicit(expr)
		if expl_expr isa AbstractArray
			math_exprs = vcat(math_exprs, expl_expr)
		else
			push!(math_exprs, expl_expr)
		end
	end

	#Transforming our Expression into a dependency tree
	tree = Arrays.to_deptree!(Set{String}(["zeta"]), math_exprs)
	
	#Transforming our dependency tree into a sequence of expressions to compute
	seq = Arrays.to_sequence!(tree)
	
	#Generating the final function
	func!, funcstr = Arrays.to_kernel(seq)
	
	if verbose
		println("Developped expression :")
		println(string(exprs))
		
		println("Tree view")
		println(string(tree))
		println("Graphviz view of tree")
		println(Arrays.to_graphviz(tree))

		println("Corresponding Sequence :")
		println(string(seq)*"\n")

		println("Generated code :")
		println(funcstr)
	end
	
	return func!
end
