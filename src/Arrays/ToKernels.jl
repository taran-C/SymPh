using ManagedLoops: @loops, @vec
using SIMD

export to_kernel

#TODO define each loop as a function used in the returned one to allow for passing arguments (strongly typed) + multi-threading with e-g ManagedLoops (or idk which other one) 
#TODO THIS IS TERRIBLE ! (?) ACTUALLY CONSTRUCT SOMETHING INSTEAD OF RELYING ON STRINGS
function to_kernel(seq::Sequence)
	calls = []

	vars = get_terms(seq)
	
	for b in seq.blocks
		if b isa CallBlock
			push!(calls, (b.expr.func,:call))
		elseif b isa LoopBlock
			push!(calls, (eval(Meta.parse(generate_loop_call(seq, vars, b))), :loop))
		end
	end

	function kernel!(mesh, state; var_repls=[])	
		#Generate kwargs from state and terms
		args = []
		kwargs = []

		for var in vars
			if var in keys(var_repls)
				push!(kwargs, Pair(Symbol(var), getproperty(state, Symbol(var_repls[var]))))
				push!(args, getproperty(state, Symbol(var_repls[var])))
			else
				push!(kwargs, Pair(Symbol(var), getproperty(state, Symbol(var))))
				push!(args, getproperty(state, Symbol(var)))
			end
		end

		for call in calls
			if call[2] == :loop
				call[1](mesh.mgr, mesh, args...)
			elseif call[2] == :call
				call[1](mesh; kwargs...)
			end
		end
	end

	return kernel!, vars
end

"""
generate_loop_call(seq, block)

	Generate a function performing one of our blocks
"""
function generate_loop_call(seq, vars, block)
	str = "@loops function compute_$(join(keys(block.exprs), "_"))(_, mesh::Mesh, "
	
	for var in vars
		str = str*var*"::Matrix{Float64}, "
	end

	str = chop(str; tail = 2)
	str = str * ")\n"

	str = str * "\tlet (irange, jrange) = (1+mesh.nh:mesh.nx-mesh.nh, 1+mesh.nh:mesh.ny-mesh.nh)\n\t\t@vec for i in irange, j in jrange\n"
			
	for key in keys(block.exprs)
		str = str * "\t\t\t" * key * "[i, j] = $(string(block.exprs[key]))\n"
	end
	str = str * "\t\tend\n\tend\nend\n"
	
	return str
end
