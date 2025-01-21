export Block, Sequence

struct Block
	exprs::Dict{String, Expression}
end

function get_terms(b::Block)
	return Set(collect(keys(b.exprs)))
end

function string(b::Block)
	s = "["
	for key in keys(b.exprs)
		s = s*key*" = "*string(b.exprs[key])*",\n "
	end
	if length(b.exprs)>0
		s = chop(s; tail=3)
	end
	s = s*"]"
	return s
end

struct Sequence
	vars::Vector{String}
	blocks::Vector{Block}
end

function get_terms(seq::Sequence)
	terms = Set(seq.vars)
	for b in seq.blocks
		terms = union(terms, get_terms(b))
	end
	return terms
end

function string(seq::Sequence)
	s = "{vars : $(string(seq.vars))\n"
	for block in seq.blocks
		s = s*string(block)*";\n\n"
	end
	if length(seq.blocks)>0
		s = chop(s; tail=3)
	end
	s = s*"}"
	return s
end
