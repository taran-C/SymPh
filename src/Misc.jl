export @Let
"""
@Let name = foo

Evaluates the expression `foo` and sets its name to `name`, then assigns a local variable named appropriately
"""
macro Let(expr)
	@assert expr.head == :(=) "the let macro should only be used on assignations"
	
	rhs = expr.args[2]
	
	#Checking for already existing kwargs
	kw = false
	if length(rhs.args) > 1 
		if typeof(rhs.args[2]) == Expr 
			if typeof(rhs.args[2].args[1]) == Expr 
				if rhs.args[2].args[1].head == :kw
					append!(rhs.args[2].args[1].args, [string(expr.args[1], true)])
					kw = true
				end
			end
		end			
	end
	
	#Inserting as first argument (right after kwargs)
	#insert!(expr.args[2].args, id, string(expr.args[1]))
	#insert!(expr.args[2].args, id+1, true)
	if !kw
		insert!(rhs.args, 2, Expr(:parameters, Expr(:kw, :name, string(expr.args[1])), Expr(:kw, :save, true)))
	end

	return :($(esc(expr.args[1])) = $(esc(expr)))
end
