export to_kernel

#TODO THIS IS TERRIBLE ! (?) ACTUALLY CONSTRUCT SOMETHING INSTEAD OF RELYING ON STRINGS
function to_kernel(seq::Sequence)
	str = "(;nh,nx,ny, "

	for term in get_terms(seq)
		str = str*term*", "
	end
	str = chop(str; tail = 2)
	str = str * ") -> begin\n"
	for b in seq.blocks
		str = str * "\tfor i in 1+nh:nx-nh, j in 1+nh:ny-nh\n"
		
		for key in keys(b.exprs)
			str = str * "\t\t" *key * "[i,j] = $(string(b.exprs[key]))\n"
		end
		str = str * "\tend\n\n"
	end

	str = str * "end\n"
	
	return eval(Meta.parse(str)), str
end
