"""
	BoundaryCondition

Objects that encodes what happens to variables at the top/bottom/left/right of the domain
"""
struct BoundaryCondition
	name::Symbol
	top
	bottom
	left
	right
end

function apply_bc(state, bc::BoundaryCondition)

end
