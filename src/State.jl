#TODO check if this is compatible with multi-threading (should be) (since the array is gathered before being used in the subfunction parameters)
struct State
	size #To be replaced by mesh object
	name
	fields::Dict{Symbol, Array{Float64,2}}
	function State(nx, ny, name)
		return new((nx,ny), name, Dict{Symbol, Array{Float64,2}}())
	end
end
function Base.getproperty(obj::State, sym::Symbol)
	fields = getfield(obj, :fields)
	size = getfield(obj, :size)

	if !(sym in fieldnames(State))
		if sym in keys(fields)
			return fields[sym]
		else
			fields[sym] = zeros(size...) #TODO replace by an alloc function in the mesh
			return fields[sym]
		end
	else
		return getfield(obj, sym)
	end
end
	
