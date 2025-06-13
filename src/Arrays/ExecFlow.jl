export Block, Sequence

"""
	Block
	
Any operation or block of operations that can be done at once
"""
abstract type Block end

"""
	LoopBlock
	
A block of operations to be computed in a loop on our domain
"""
struct LoopBlock
	exprs::Dict{String, Expression}
end

"""
	get_terms(b::LoopBlock)

Get a set of string representing the arrays being assigned in the loop
"""
function get_terms(b::LoopBlock)
	return unique(collect(keys(b.exprs)))
end

function string(b::LoopBlock)
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

"""
	CallBlock
	
A block representing a single function call
"""
struct CallBlock
	name::String
	expr::FuncCall
end

function get_terms(b::CallBlock)
	return Set((b.name,))
end

string(b::CallBlock) = "[" * string(b.expr) * "]"

"""
	Sequence
	
A sequence of blocks of operations
"""
struct Sequence
	vars::Vector{String}
	blocks#::Vector{Block}
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
