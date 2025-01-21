export DepNode
export addchild!, shave!, to_sequence!

#TODO check if AbstractTree.jl can be useful
mutable struct DepNode
	name::String
	expr::Union{Nothing, Arrays.Expression}
	childs::Set{DepNode}
	parent::Union{Nothing, DepNode}
	function DepNode(name::String, expr::Union{Nothing, Arrays.Expression})
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

function shave!(tree::DepNode; exprs = Dict{String, Arrays.Expression}())
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
		push!(blocks, Arrays.Block(exprs))
	end

	return Arrays.Sequence(vars, blocks)
end
