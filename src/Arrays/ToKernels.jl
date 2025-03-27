export to_kernel

#TODO define each loop as a function used in the returned one to allow for passing arguments (strongly typed) + multi-threading with e-g ManagedLoops (or idk which other one) 
#TODO THIS IS TERRIBLE ! (?) ACTUALLY CONSTRUCT SOMETHING INSTEAD OF RELYING ON STRINGS
function to_kernel(seq::Sequence)
	#Signature
	str = "(mesh::Mesh; "

	for term in get_terms(seq)
		str = str*term*"::AbstractArray{Float64}, "
	end

	str = chop(str; tail = 2)
	str = str * ") -> begin\n"

	#Mesh fields forwarding
	#for fieldname in fieldnames(Mesh)
	#	str = str * "\t$(fieldname)=mesh.$(fieldname)\n"
	#end
	#str = str*"\n"

	#Computation loops
	for b in seq.blocks
		if b isa CallBlock
			str = str * "\t"*b.expr.func*"("*b.expr.name*", "
			for arg in b.expr.args
				str = str * arg.name * ", "
			end
			str = chop(str, head = 0, tail = 2)
			str = str * ")\n\n"
		elseif b isa LoopBlock
			str = str * "\tfor i in 1+mesh.nh:mesh.nx-mesh.nh, j in 1+mesh.nh:mesh.ny-mesh.nh\n"
			
			for key in keys(b.exprs)
				str = str * "\t\t" *key * "[i,j] = $(string(b.exprs[key]))\n"
			end
			str = str * "\tend\n\n"
		end
	end

	str = str * "end\n"
		
	return eval(Meta.parse(str)), str, get_terms(seq)
end
