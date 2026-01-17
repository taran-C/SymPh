export DepNode
export addchild!, shave!, to_sequence!, to_deptree!

#TODO check if AbstractTree.jl can be useful
"""
	DepNode(name::String, expr::Union{Nothing, Expression})

A node of a dependency tree of expressions made of expressions
"""
mutable struct DepNode
	name::String
	expr::Union{Nothing, Expression}
	childs::Set{DepNode}
	parent::Union{Nothing, DepNode}
	function DepNode(name::String, expr::Union{Nothing, Expression})
		return new(name, expr, Set{DepNode}(), nothing)
	end
end

function addchild!(parent, child)
	#Exclude repeated dependencies TODO check if necessary and most importantly if not dangerous (depends on unicity of names)
	for c in parent.childs
		if c.name == child.name
			return
		end
	end

	push!(parent.childs, child)
	child.parent = parent
end

function removechild!(parent, child)
	delete!(parent.childs, child)
	child.parent = nothing
end

function shave!(tree::DepNode; exprs = Dict{String, Expression}())
	if length(tree.childs) == 0
		removechild!(tree.parent, tree)
		exprs[tree.name] = tree.expr
	else
		for child in tree.childs
			shave!(child; exprs=exprs)
		end
	end

	return exprs
end

#TODO this is destructive for the dependency tree, should it be ?
function to_sequence!(tree::DepNode)
	#Elements of mesh are not variables but are forwarded
	vars = setdiff(collect(keys(shave!(tree))), "mesh." .* string.(fieldnames(Mesh)))
	blocks = []
	
	while length(tree.childs) > 0
		exprs = shave!(tree)

		loopexprs = Dict{String, Expression}()

		for expr in exprs
			if expr.second isa FuncCall
				push!(blocks, CallBlock(expr.first,expr.second))
			else
				loopexprs[expr.first] = expr.second
			end
		end
		if length(loopexprs) > 0
			push!(blocks, LoopBlock(loopexprs))
		end
	end

	return Sequence(vars, blocks)
end

#TODO fix bug where if the name of the result is present in node_names it is recalculated from its own value
to_deptree!(node_names::Set{String}, expr::Expression) = to_deptree!(node_names, (expr,)) #unravel to allow flexibility

function to_deptree!(exprs)
	root = DepNode("root", nothing)
	
	for expr in exprs
		expr.save = false
		expr_to_node!(expr, root)
	end

	return root
end	

function expr_to_node!(expr::Expression, parent)
	if expr isa ArrayVariable
		#This is ugly and this whole way of creating the tree is kinda terrible
		n = DepNode(expr.name, expr[-expr.depi, -expr.depj])
		addchild!(parent, n)
		
		return n
	else 
		n = DepNode(expr.name, nothing)
		n.expr = go_deeper(expr, n)
		n.expr = n.expr[-n.expr.depi, -n.expr.depj]
		addchild!(parent, n)

		return n
	end
end

#TODO set depi/depj to 0 when creating a new node (so node isn't offset WRT itself)
function go_deeper(expr::Expression, parent)
	if expr isa FuncCall
		n = DepNode(expr.name, nothing)
		n.expr = deepcopy(expr)
		n.expr = n.expr[-n.expr.depi, -n.expr.depj]
		addchild!(parent, n)
		for (i, arg) in enumerate(n.expr.args)
			expr_to_node!(arg, n)
			n.expr.args[i] = ArrayVariable(arg.name, false, arg.depi, arg.depj)
		end
		return ArrayVariable(expr.name, false, expr.depi, expr.depj)
	elseif expr.save
		expr = deepcopy(expr)
		expr.save = false
		expr_to_node!(expr, parent)
		return ArrayVariable(expr.name, false, expr.depi, expr.depj)
	elseif expr isa BinaryOperator
		return typeof(expr)(expr.name, expr.save, go_deeper(expr.left, parent), go_deeper(expr.right, parent), expr.depi, expr.depj)
	elseif expr isa UnaryOperator
		return typeof(expr)(expr.name, expr.save, go_deeper(expr.expr, parent), expr.depi, expr.depj)
	elseif expr isa TernaryOperator
		return TernaryOperator(expr.name, expr.save, go_deeper(expr.a, parent), go_deeper(expr.b, parent), go_deeper(expr.c, parent), expr.depi, expr.depj)
	elseif expr isa ArrayVariable
		expr_to_node!(expr, parent)
		return ArrayVariable(expr.name, expr.save, expr.depi, expr.depj)
	else
		return expr
	end
end

function to_graphviz(tree::DepNode)
	s=""
	for c in tree.childs
		if length(c.name) < 4 || SubString(c.name, 1, 4) != "mesh"
			s = s*"$(tree.name) -> $(c.name)\n"
			s = s*to_graphviz(c)
		end
	end
	return s
end

function string(tree::DepNode, lev = 0; show_expr = false)
	#keeping track of indentation level and using it before each expression
		s = "$("\t"^lev)$(tree.name)\n"
	if show_expr && tree.expr != nothing
		s = "$("\t"^(lev+1))expr : $(string(tree.expr)) @ [$(tree.expr.depi),$(tree.expr.depj)]\n"
	end

	if length(tree.childs) > 0
		s *= "$("\t"^(lev+1))childs:\n"
	end
	for c in tree.childs
		if length(c.name) < 4 || SubString(c.name, 1, 4) != "mesh"
			s = s*string(c, lev+2)
		end
	end

	return s
end
