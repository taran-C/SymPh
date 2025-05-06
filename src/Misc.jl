export @Let
"""
@Let name foo

Evaluates the expression foo and sets its name to name, then assigns a local variable named appropriately
"""
macro Let(expr)
	@assert expr.head == :(=) "the let macro should only be used on assignations"
	
	id = 2
	#Checking for kwargs
	rhs = expr.args[2]
	if length(rhs.args) > 1
		if typeof(rhs.args[2]) == Expr
			if typeof(rhs.args[2].args[1]) == Expr
				if rhs.args[2].args[1].head == :kw
					id = 3
				end
			end
		end
	end
	
	#Inserting as first argument (right after kwargs)
	insert!(expr.args[2].args, id, string(expr.args[1]))
	#dump(expr)
	return :($(esc(expr.args[1])) = $(esc(expr)))
end
