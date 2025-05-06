export rk3step!, rk4step!, euler_forwardstep!

"""
rk3step!(dt, mesh, state, progs)
	
	dt : time increment
	mesh : Mesh object
	state : State object
	progs : Name of prognostic variables
"""
function rk3step!(rhs!, dt, mesh, state, progs)	
	
	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "1"))
	for p in progs
		getproperty(state, Symbol(p)) .+= dt * getproperty(state, Symbol("dt" * p * "1"))
	end

	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "2"))
	for p in progs
		getproperty(state, Symbol(p)) .-= dt * 3/4 * getproperty(state, Symbol("dt" * p * "1"))
		getproperty(state, Symbol(p)) .+= dt * 1/4 * getproperty(state, Symbol("dt" * p * "2"))
	end

	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "3"))
	for p in progs
		getproperty(state, Symbol(p)) .-= dt * 1/12 * getproperty(state, Symbol("dt" * p * "1"))
		getproperty(state, Symbol(p)) .-= dt * 1/12 * getproperty(state, Symbol("dt" * p * "2"))
		getproperty(state, Symbol(p)) .+= dt * 2/3 * getproperty(state, Symbol("dt" * p * "3"))
	end
end

"""
rk4step!(dt, mesh, state, progs)
	
	dt : time increment
	mesh : Mesh object
	state : State object
	progs : Name of prognostic variables
"""
function rk4step!(rhs!, dt, mesh, state, progs)	
	
	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "1"))
	for p in progs
		getproperty(state, Symbol(p)) .+= dt * 1/2 * getproperty(state, Symbol("dt" * p * "1"))
	end

	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "2"))
	for p in progs
		getproperty(state, Symbol(p)) .-= dt * 1/2 * getproperty(state, Symbol("dt" * p * "1"))
		getproperty(state, Symbol(p)) .+= dt * 1/2 * getproperty(state, Symbol("dt" * p * "2"))
	end

	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "3"))
	for p in progs
		getproperty(state, Symbol(p)) .-= dt * 1/2 * getproperty(state, Symbol("dt" * p * "2"))
		getproperty(state, Symbol(p)) .+= dt * getproperty(state, Symbol("dt" * p * "3"))
	end
	
	rhs!(mesh, state; var_repls = Dict{String, String}("dt" .* progs .=> "dt" .* progs .* "4"))
	for p in progs
		getproperty(state, Symbol(p)) .+= dt * 1/6 * getproperty(state, Symbol("dt" * p * "1"))
		getproperty(state, Symbol(p)) .+= dt * 1/3 * getproperty(state, Symbol("dt" * p * "2"))
		getproperty(state, Symbol(p)) .-= dt * 2/3 * getproperty(state, Symbol("dt" * p * "3"))
		getproperty(state, Symbol(p)) .+= dt * 1/6 * getproperty(state, Symbol("dt" * p * "4"))
	end
end

"""
euler_forwardstep!(dt, mesh, state, progs)

	dt : time increment
	mesh : Mesh object
	state : State object
	progs : Name of prognostic variables
"""
function euler_forwardstep!(rhs!, dt, mesh, state, progs)
	rhs!(mesh, state)
	for p in progs
		getproperty(state, Symbol(p)) .+= dt * getproperty(state, Symbol("dt" * p))
	end
end
