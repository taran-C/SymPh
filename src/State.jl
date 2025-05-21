#TODO check if this is compatible with multi-threading (should be) (since the array is gathered before being used in the subfunction parameters)
export State

mutable struct State
	mesh
	fields::Dict{Symbol, Array{Float64,2}}
	function State(mesh)
		return new(mesh, Dict{Symbol, Array{Float64,2}}())
	end
end
function Base.getproperty(obj::State, sym::Symbol)
	fields = getfield(obj, :fields)
	mesh = getfield(obj, :mesh)
	
	if !(sym in fieldnames(State))
		if sym in keys(fields)
			return fields[sym]
		else
			fields[sym] = zeros(Float64, mesh.ni, mesh.nj)
			return fields[sym]
		end
	else
		return getfield(obj, sym)
	end
end

export reset_state

function reset_state(state::State)
	state.fields = Dict{Symbol, Array{Float64,2}}()
end
