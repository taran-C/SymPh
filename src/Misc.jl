export @Let
"""
@Let name foo

Evaluates the expression foo and sets its name to name, then assigns a local variable named appropriately
"""
macro Let(expr)
	@assert expr.head == :(=) "the let macro should only be used on assignations"
	insert!(expr.args[2].args, 2, string(expr.args[1]))
	return :($(esc(expr.args[1])) = $(esc(expr)))
end
