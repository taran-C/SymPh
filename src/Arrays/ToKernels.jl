export to_kernel

#TODO define each loop as a function used in the returned one to allow for passing arguments (strongly typed) + multi-threading with e-g ManagedLoops (or idk which other one) 
#TODO THIS IS TERRIBLE ! (?) ACTUALLY CONSTRUCT SOMETHING INSTEAD OF RELYING ON STRINGS
function to_kernel(seq::Sequence)
	calls = []

	for b in seq.blocks
		if b isa CallBlock
			push!(calls, b.expr.func)
		elseif b isa LoopBlock
			push!(calls, eval(Meta.parse(generate_loop_call(seq, b))))
		end
	end

	vars = get_terms(seq)

	function kernel!(mesh, state; var_repls=[])	
		#Generate kwargs from state and terms
		kwargs = []

		for var in vars
			if var in keys(var_repls)
				push!(kwargs, Pair(Symbol(var), getproperty(state, Symbol(var_repls[var]))))
			else
				push!(kwargs, Pair(Symbol(var), getproperty(state, Symbol(var))))
			end
		end

		for call in calls
			call(mesh; kwargs...)
		end
	end

	return kernel!, vars
end

"""
generate_loop_call(seq, block)

	Generate a function performing one of our blocks
"""
function generate_loop_call(seq, block)
	str = "(mesh::Mesh; "

	for term in get_terms(seq)
		str = str*term*"::Matrix{Float64}, "
	end

	str = chop(str; tail = 2)
	str = str * ") -> begin\n"

	str = str * "\tfor i in 1+mesh.nh:mesh.nx-mesh.nh, j in 1+mesh.nh:mesh.ny-mesh.nh\n"
			
	for key in keys(block.exprs)
		str = str * "\t\t" *key * "[i,j] = $(string(block.exprs[key]))\n"
	end
	str = str * "\tend\nend\n"
end
