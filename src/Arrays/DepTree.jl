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
	root = DepNode("root", nothing)
	addchild!(root, tree)

	vars = collect(keys(shave!(tree)))
	blocks = []
	
	while length(root.childs) > 0
		exprs = shave!(tree)
		push!(blocks, Block(exprs))
	end

	return Sequence(vars, blocks)
end

function to_deptree!(expr::Expression, node_names::Set{String}; parent=nothing)
	if expr isa ArrayVariable
		n = DepNode(expr.name, expr)
		if parent == nothing
			return n
		else
			addchild!(parent, n)
			
		end
	else
		if parent == nothing #Root of our tree
			n = DepNode(expr.name, nothing)
			n.expr = go_deeper(expr, node_names, n)
			return n	
		else
			n = DepNode(expr.name, nothing)
			n.expr = go_deeper(expr, node_names, n)
			addchild!(parent, n)
		end
	end
end

function go_deeper(expr::Expression, node_names::Set{String}, parent)
	if in(expr.name, node_names)
		nnames = copy(node_names)
		delete!(nnames, expr.name)	
		to_deptree!(expr, nnames; parent)
		return ArrayVariable(expr.name)
	elseif expr isa BinaryOperator
		return typeof(expr)(expr.name, go_deeper(expr.left, node_names, parent), go_deeper(expr.right, node_names, parent))
	elseif expr isa UnaryOperator
		return typeof(expr)(expr.name, go_deeper(expr.expr, node_names, parent))
	elseif expr isa ArrayVariable
		to_deptree!(expr, node_names; parent)
		return ArrayVariable(expr.name)
	else
		return expr
	end
end

