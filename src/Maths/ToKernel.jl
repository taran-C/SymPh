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
	fill = Dict()
	for bc in bcs
		if bc isa Vect
			if primality(bc) == Dual
				fill[bc.name * "_X"] = (+1, 0, 0, 0)
				fill[bc.name * "_Y"] = (0, 0, +1, 0)
			elseif primality(bc) == Primal
				fill[bc.name * "_X"] = (0, 0, -1, 0)
				fill[bc.name * "_Y"] = (-1, 0, 0, 0)
			end
		elseif degree(bc) == 0
			if primality(bc) == Dual
				fill[bc.name] = (0, 0, 0, 0) # second part : (ldec, rdec, bdec, tdec), offset wrt main msk, eh 
			elseif primality(bc) == Primal
				fill[bc.name] = (-1, 0, -1, 0)
			end
		elseif degree(bc) == 1
			if primality(bc) == Dual
				fill[bc.name * "_i"] = (+1, 0, 0, 0)
				fill[bc.name * "_j"] = (0, 0, +1, 0)
			elseif primality(bc) == Primal
				fill[bc.name * "_i"] = (0, 0, -1, 0)
				fill[bc.name * "_j"] = (-1, 0, 0, 0)
			end
		elseif degree(bc) == 2
			if primality(bc) == Dual
				fill[bc.name] = (+1, 0, +1, 0)
			elseif primality(bc) == Primal
				fill[bc.name] = (0, 0, 0, 0)
			end
		end
	end
	
	#Generating the final function
	func!,  vars = Arrays.to_kernel(seq, fill)

	return func!
end
