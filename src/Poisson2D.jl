using LinearAlgebra, SparseArrays

function laplacian(mesh, bc)
	@assert bc in ["dirichlet", "neumann"]
	
	G = zeros(Int32, mesh.nx, mesh.ny)
	G[mesh.msk .> 0] .= collect(1:length(G[mesh.msk .> 0]))
	G[mesh.msk .== 0] .= -1
	
	#TODO adapt to variable dx/dy
	dx2 = 1/mesh.dx ^ 2
	dy2 = 1/mesh.dy ^ 2
	
	N = sum(G .> -1)

	vals = zeros(0)
	rows = zeros(0)
	cols = zeros(0)
	counter = 1

	function add_entry(val, r, c, counter)
		push!(vals, val)
		push!(rows, r)
		push!(cols, c)
		return counter +1
	end

	for j in 1:mesh.ny
		for i in 1:mesh.nx
			I = G[i,j]

			if I>-1
				s = 0
				
				west = i>1 ? G[i-1, j] : -1 #(G[j, -1] if mesh.xperiodic else -1)
				east = i<mesh.nx-1 ? G[i+1, j] : -1 #(G[j, 0] if mesh.xperiodic else -1)
				south = j>1 ? G[i, j-1] : -1 #(G[-1, i] if mesh.yperiodic else -1)
				north = j<mesh.ny-1 ? G[i, j+1] : -1 #(G[0, i] if mesh.yperiodic else -1)
				
				if west > -1
					counter = add_entry(dx2, I, west, counter)
					s += dx2
				end

				if east > -1
					counter = add_entry(dx2, I, east, counter)
					s += dx2
				end

				if north > -1
					counter = add_entry(dy2, I, north, counter)
					s += dy2
				end

				if south > -1
					counter = add_entry(dy2, I, south, counter)
					s += dy2
				end

				if bc == "dirichlet"
					counter = add_entry(-2*(dx2+dy2), I, I, counter)
				elseif bc == "neumann"
					counter = add_entry(-s, I, I, counter)
				end
			end
		end
	end

	return factorize(sparse(rows, cols, vals, N, N))
end

"""
	solve_poisson(out, b)

Solves the poisson problem Ax = b and sets out to be x
"""
function get_poisson_solver(mesh, bc)
	@time A = laplacian(mesh, bc)
	
	#Create a function that solves Ax = b
	function func!(out, b)
		out[mesh.msk .> 0] = A\b[mesh.msk .> 0]
	end

	return func!
end
