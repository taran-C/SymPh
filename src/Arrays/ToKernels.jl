using ManagedLoops: @loops, @vec
using SIMD

export to_kernel

"""
	Converts a sequence into a computing kernel
"""
function to_kernel(seq::Sequence, fill)
	calls = []

	vars = get_terms(seq)
	
	for b in seq.blocks
		if b isa CallBlock
			push!(calls, (b.expr.func,:call, [b.name]))
		elseif b isa LoopBlock
			str, keys = generate_loop_call(seq, vars, b)
			push!(calls, (eval(Meta.parse(str)), :loop, keys))
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
			
			#Halo filling
			for key in call[3]
				if haskey(fill, key)
					nx = mesh.nx
					ny = mesh.ny
					nh = mesh.nh
					
					if key in keys(var_repls)
						q = getproperty(state, Symbol(var_repls[key]))
					else
						q = getproperty(state, Symbol(key))
					end
					
					if mesh.xperio
						copy_x!(q, nx, ny, nh, fill[key])
					end
					if mesh.yperio
						copy_y!(q, nx, ny, nh, fill[key])
					end
				end
			end
		end
	end

	return kernel!, vars
end

#dec : (ldec, rdec, bdec, tdec) TODO only works for dual grid rn (need to be able to copy less than full halo for bigger grid)
function copy_x!(q, nx, ny, nh, dec)
        #horizontal edges
        for i = 1:nh, j = 1:ny
                q[i+dec[1],j] = q[i+nx-2*nh+dec[2], j]
                q[nx-nh+i-dec[2],j] = q[i+nh+dec[1], j]
        end
end
function copy_y!(q, nx, ny, nh, dec)
        #vertical edges
        for i = 1:nx, j = 1:nh
                q[i,j+dec[3]] = q[i, j+ny-2*nh+dec[4]]
                q[i,ny-nh+j-dec[4]] = q[i, j+nh+dec[3]]
        end
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
	
	return str, keys(block.exprs)
end
