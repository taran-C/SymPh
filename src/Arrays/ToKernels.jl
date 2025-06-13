using ManagedLoops: @loops, @vec
using SIMD

export to_kernel

"""
	to_kernel(seq::Sequence, fill; verbose = false)

Converts `seq` into a computing kernel, `fill` is the list of names of value that need their halo filled.
"""
function to_kernel(seq::Sequence, fill; verbose = false)
	calls = []

	vars = get_terms(seq)
	
	for b in seq.blocks
		if b isa CallBlock
			push!(calls, (b.expr.func,:call, [b.name]))
			if verbose
				println("Function call to " * b.name)
			end
		elseif b isa LoopBlock
			str, keys = generate_loop_call(seq, vars, b)
			if verbose
				println(str)
			end
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
				if key in fill
					ni = mesh.ni
					nj = mesh.nj
					nh = mesh.nh
					
					if key in keys(var_repls)
						q = getproperty(state, Symbol(var_repls[key]))
					else
						q = getproperty(state, Symbol(key))
					end
					
					#=	
					if key == "dtb"
						neumann_center(q, ni, nj, nh)
					end
					
					if key == "dhb_j"
						dirichlet_jedge_dual(q, ni, nj, nh; t = 0.03125, b = 0.03125, l = 0.03125, r = 0.03125)
					end
					if key == "dhb_i"
						dirichlet_iedge_dual(q, 0, 0, ni, nj, nh)
					end

					if key == "u_i"
						neumann_iedge_dual(q, ni, nj, nh)
					end
					if key == "u_j"
						neumann_jedge_dual(q, ni, nj, nh)
					end
					if key == "U_X"
						neumann_iedge_dual(q, ni, nj, nh)
					end
					if key == "U_Y"
						neumann_jedge_dual(q, ni, nj, nh)
					end
					=#

					if mesh.iperio
						copy_i!(q, ni, nj, nh)
					end
					if mesh.jperio
						copy_j!(q, ni, nj, nh)
					end
				end
			end
		end
	end

	return kernel!, vars
end

#Boundary Conditions-------------------------------------------------------
#TODO dispatch on form types, separate left from right and top from bottom
function dirichlet_i_center(q, ni, nj, nh; l=0, r=0, t=0, b=0)
	#i
	for j in nh+1:nj-nh
		q[nh, j] = 2*a - q[nh+1, j]
		q[ni-nh+1,j] = 2*b - q[ni-nh, j]
	end
	#j
	for i in nh+1:ni-nh
		q[i,nh] = 2*a - q[i, nh+1]
		q[i,nj-nh+1] = 2*b - q[i,nj-nh]
	end
end

function dirichlet_iedge_dual(q, a, b, ni, nj, nh)
	for i in nh+1:ni-nh
		q[i,1:nh] .= 2*a - q[i, nh+1]
		q[i,nj-nh+1:end] .= 2*b - q[i,nj-nh]
	end
	for j in nh+1:nj-nh
		q[1:nh+1, j] .= 2*a - q[nh+2, j]
		q[ni-nh+1:end,j] .= 2*b - q[ni-nh, j]
	end
end
function dirichlet_jedge_dual(q, ni, nj, nh; l = 0, r = 0, t = 0, b=0)
	for i in nh+1:ni-nh
		q[i,1:nh+1] .= 2*b - q[i, nh+2]
		q[i,nj-nh+1:end] .= 2*t - q[i,nj-nh]
	end
	for j in nh+1:nj-nh
		q[1:nh, j] .= 2*l - q[nh+1, j]
		q[ni-nh+1:end,j] .= 2*r - q[ni-nh, j]
	end
end

function neumann_center(q, ni, nj, nh; l = 0, r = 0, t = 0, b = 0)
	#i
	for j in nh+1:nj-nh
		q[1:nh, j] .= q[nh+1, j] - l
		q[ni-nh+1:end,j] .= q[ni-nh, j] + r
	end

	#j
	for i in nh+1:ni-nh
		q[i,1:nh] .= q[i, nh+1] - b
		q[i,nj-nh+1:end] .= q[i,nj-nh] + t
	end
end

function neumann_iedge_dual(q, ni, nj, nh; a=0, b=0)
	for i in nh+1:ni-nh
		q[i,1:nh] .= q[i, nh+1] - a
		q[i,nj-nh+1:end] .= q[i,nj-nh] + b
	end
	for j in nh+1:nj-nh
		q[1:nh+1, j] .= q[nh+2, j] - a
		q[ni-nh+1:end,j] .= q[ni-nh, j] + b
	end
end
function neumann_jedge_dual(q, ni, nj, nh; l = 0, r = 0, t = 0, b=0)
	for i in nh+1:ni-nh
		q[i,1:nh+1] .= q[i, nh+2] - b
		q[i,nj-nh+1:end] .= q[i,nj-nh] + t
	end
	for j in nh+1:nj-nh
		q[1:nh, j] .= q[nh+1, j] - l
		q[ni-nh+1:end,j] .= q[ni-nh, j] + r
	end
end

function copy_i!(q, ni, nj, nh)
        #horizontal edges
        for i = 1:nh, j = 1:nj
                q[i,j] = q[i+ni-2*nh, j]
                q[ni-nh+i,j] = q[i+nh, j]
        end
end
function copy_j!(q, ni, nj, nh)
        #vertical edges
        for i = 1:ni, j = 1:nh
                q[i,j] = q[i, j+nj-2*nh]
                q[i,nj-nh+j] = q[i, j+nh]
        end
end
#--------------------------------------------------------------------


"""
	generate_loop_call(seq, vars, block)

Generate a function performing one of our blocks
"""
function generate_loop_call(seq, vars, block)
	str = "@loops function compute_$(join(keys(block.exprs), "_"))(_, mesh::Mesh, "
	
	for var in vars
		str = str*var*"::Matrix{Float64}, "
	end

	str = chop(str; tail = 2)
	str = str * ")\n"
	str = str * "\tlet (irange, jrange) = (mesh.nh:mesh.ni-mesh.nh+1, mesh.nh:mesh.nj-mesh.nh+1)\n\t\t@vec for i in irange, j in jrange\n"
			
	for key in keys(block.exprs)
		str = str * "\t\t\t" * key * "[i, j] = $(string(block.exprs[key]))\n"
	end
	str = str * "\t\tend\n\tend\nend\n"
	
	return str, keys(block.exprs)
end
