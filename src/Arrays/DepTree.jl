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

function to_deptree!(node_names::Set{String}, exprs)
	root = DepNode("root", nothing)
	
	for expr in exprs
		expr_to_node!(expr, node_names, root)
	end

	return root
end	

function expr_to_node!(expr::Expression, node_names::Set{String}, parent)
	if expr isa ArrayVariable
		#This is ugly and this whole way of creating the tree is kinda terrible
		n = DepNode(expr.name, expr[-expr.depi, -expr.depj])
		addchild!(parent, n)
		
		return n
	else 
		n = DepNode(expr.name, nothing)
		n.expr = go_deeper(expr, node_names, n)
		n.expr = n.expr[-n.expr.depi, -n.expr.depj]
		addchild!(parent, n)

		return n
	end
end

#TODO set depi/depj to 0 when creating a new node (so node isn't offset WRT itself)
#Also WHY is there nondeterministic behavior ? -> Due to the way Sets are handled i guess, iteration on it seems to be chaotic
function go_deeper(expr::Expression, node_names::Set{String}, parent)
	if in(expr.name, node_names)
		nnames = copy(node_names)
		delete!(nnames, expr.name)	
		expr_to_node!(expr, nnames, parent)
		return ArrayVariable(expr.name, expr.depi, expr.depj)
	elseif expr isa FuncCall
		n = DepNode(expr.name, nothing)
		n.expr = deepcopy(expr)
		n.expr = n.expr[-n.expr.depi, -n.expr.depj]
		addchild!(parent, n)
		for (i, arg) in enumerate(n.expr.args)
			expr_to_node!(arg, node_names, n)
			n.expr.args[i] = ArrayVariable(arg.name, arg.depi, arg.depj)
		end
		return ArrayVariable(expr.name, expr.depi, expr.depj)
	elseif expr isa BinaryOperator
		return typeof(expr)(expr.name, go_deeper(expr.left, node_names, parent), go_deeper(expr.right, node_names, parent), expr.depi, expr.depj)
	elseif expr isa UnaryOperator
		return typeof(expr)(expr.name, go_deeper(expr.expr, node_names, parent), expr.depi, expr.depj)
	elseif expr isa TernaryOperator
		return TernaryOperator(expr.name, go_deeper(expr.a, node_names, parent), go_deeper(expr.b, node_names, parent), go_deeper(expr.c, node_names, parent), expr.depi, expr.depj)
	elseif expr isa ArrayVariable
		expr_to_node!(expr, node_names, parent)
		return ArrayVariable(expr.name, expr.depi, expr.depj)
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
