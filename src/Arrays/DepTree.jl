export DepNode
export addchild!, shave!, to_sequence!, to_deptree!

#TODO check if AbstractTree.jl can be useful
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
		n = DepNode(expr.name, expr[-expr.depx, -expr.depy])
		addchild!(parent, n)
		
		return n
	else 
		n = DepNode(expr.name, nothing)
		n.expr = go_deeper(expr, node_names, n)
		n.expr = n.expr[-n.expr.depx, -n.expr.depy]
		addchild!(parent, n)

		return n
	end
end

#TODO set depx/depy to 0 when creating a new node (so node isn't offset WRT itself)
#Also WHY is there nondeterministic behavior ? -> Due to the way Sets are handled i guess, iteration on it seems to be chaotic
function go_deeper(expr::Expression, node_names::Set{String}, parent)
	if in(expr.name, node_names)
		nnames = copy(node_names)
		delete!(nnames, expr.name)	
		expr_to_node!(expr, nnames, parent)
		return ArrayVariable(expr.name, expr.depx, expr.depy)
	elseif expr isa FuncCall
		nnames = copy(node_names)
		n = DepNode(expr.name, nothing)
		n.expr = expr
		n.expr = n.expr[-n.expr.depx, -n.expr.depy]
		addchild!(parent, n)
		for (i, arg) in enumerate(expr.args)
			expr_to_node!(arg, nnames, n)
			expr.args[i] = ArrayVariable(arg.name, arg.depx, arg.depy)
		end
		return ArrayVariable(expr.name, expr.depx, expr.depy)
	elseif expr isa BinaryOperator
		return typeof(expr)(expr.name, go_deeper(expr.left, node_names, parent), go_deeper(expr.right, node_names, parent), expr.depx, expr.depy)
	elseif expr isa UnaryOperator
		return typeof(expr)(expr.name, go_deeper(expr.expr, node_names, parent), expr.depx, expr.depy)
	elseif expr isa TernaryOperator
		return TernaryOperator(expr.name, go_deeper(expr.a, node_names, parent), go_deeper(expr.b, node_names, parent), go_deeper(expr.c, node_names, parent), expr.depx, expr.depy)
	elseif expr isa ArrayVariable
		expr_to_node!(expr, node_names, parent)
		return ArrayVariable(expr.name, expr.depx, expr.depy)
	else
		return expr
	end
end

function to_graphviz(tree::DepNode)
	s=""
	for c in tree.childs
		s = s*"$(tree.name) -> $(c.name)\n"
		s = s*to_graphviz(c)
	end
	return s
end

function string(tree::DepNode, lev = 0)
	#keeping track of indentation level and using it before each expression
	
	if tree.expr != nothing
		s = "$("\t"^lev)$(tree.name) :\n$("\t"^(lev+1))expr : $(string(tree.expr)) @ [$(tree.expr.depx),$(tree.expr.depy)]\n$("\t"^(lev+1))childs:\n"
	else 
		s = "$("\t"^lev)$(tree.name) :\n$("\t"^(lev+1))childs:\n"
	end

	for c in tree.childs
		s = s*string(c, lev+2)
	end

	return s
end
