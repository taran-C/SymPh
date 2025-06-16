using LinearAlgebra, SparseArrays

export Poisson2D, solve_poisson

mutable struct Poisson2D
	bc
	form_type
	poisson_solver
	order
end
Poisson2D(bc, form_type; order = 2) = Poisson2D(bc, form_type, nothing, order)
function solve_poisson(poisson::Poisson2D, mesh, out, b)
	if poisson.poisson_solver == nothing
		poisson.poisson_solver = get_poisson_solver(mesh, b, poisson.bc, poisson.form_type; o = poisson.order)
	end
	poisson.poisson_solver(out, b)
end
function ap_op!(A, q, mesh; out = zero(q), idx = nothing)
    #Only apply to inner points (desingularize)
    if idx == nothing
        @views out .= A * q
    else
        @views out[idx] .= A * q[idx]
    end
    return out
end
using LinearSolve
using Statistics
function ap_inv!(A, q, mesh; out = zero(q), idx = nothing)
    	#Only apply to inner points (desingularize)
    	if idx == nothing
        	@views out .= A \ q
    	else
		#TODO ALLOCATES, BAD, use a workspace
		out[idx] .= Krylov.solution(workspace)
	end
	#display(mean(out[idx]))
	#out[idx] .-= mean(out[idx])
	return out
end
function desingularize_op(A)
    idx = findall(!=(0),diag(A))
    N = size(A)[1]
    N1 = length(idx)
    
    row = idx
    col = collect(1:N1)
    data = ones(N1)
    
    J = sparse(row, col, data, N, N1)

    return J' * A * J, idx
end

function get_op_from_coeffs(mesh, coeffs; factor=nothing, idec = 0, jdec = 0)
	ni,nj,nh = mesh.ni, mesh.nj, mesh.nh
	N = ni * nj
	Is = zeros(0)
	Js = zeros(0)
	Vs = zeros(0)

	#ONLY FOR VERTICES, REARRANGE TO ACTUALLY PASS MASK TO LAPLACIAN
	idec = 1
	jdec = 1

	function add_entry(I, J, V)
		push!(Is, I)
		push!(Js, J)
		push!(Vs, V)
	end
	ijtoI(i,j) = (j-1) * ni + i
    
	for i in 1:ni, j in 1:nj
		I = ijtoI(i,j)

		#Check if all points in domain
		all_in = true
		for point in coeffs
			if ((i + point[1][1] < nh+1 + idec) & !mesh.iperio) || 
				((j + point[1][2] < nh+1 + jdec) & !mesh.jperio) || 
				((i + point[1][1] > ni-nh) & !mesh.iperio) || 
				((j + point[1][2] > nj-nh) & !mesh.jperio)
				all_in = false
			end
		end

		#Add the points (no need to check periodicity here, it is taken into account higher)
		if all_in
			for point in coeffs
				if (i + point[1][1] > nh) & (i + point[1][1] < ni-nh+1)
					ieff = i + point[1][1]
				elseif i + point[1][1]<= nh
					ieff = i + point[1][1] + (ni-2*nh)
				elseif i + point[1][1] >= ni-nh+1
					ieff = i + point[1][1] - (ni-2*nh)
				end
				if (j + point[1][2] > nh) & (j + point[1][2] < nj-nh+1)
					jeff = j + point[1][2]
				elseif j + point[1][2] <= nh
					jeff = j + point[1][2] + (nj-2*nh)
				elseif j + point[1][2] >= nj-nh+1
					jeff = j + point[1][2] - (nj-2*nh)
				end
                
				JJ = (jeff-1) * ni + ieff
                
				if factor == nothing
					fac = 1
				else
					fac = factor[i,j]
				end
				add_entry(I, JJ, fac * point[2])
			end
		end
	end
	return sparse(Is, Js, Vs, N, N)
end

#FVtoFD TODO allow choice
Ii(mesh) = get_op_from_coeffs(mesh, Dict((-1,0)=>0, (0,0)=>1, (1,0)=>0))
Iim(mesh) = get_op_from_coeffs(mesh, Dict((-1,0)=>-1/24, (0,0)=>25/24))
Iip(mesh) = get_op_from_coeffs(mesh, Dict((1,0)=>-1/24, (0,0)=>25/24))

Ij(mesh) = get_op_from_coeffs(mesh, Dict((0,-1)=>0, (0,0)=>1, (0,1)=>0))
Ijm(mesh) = get_op_from_coeffs(mesh, Dict((0,-1)=>-1/24, (0,0)=>25/24))
Ijp(mesh) = get_op_from_coeffs(mesh, Dict((0,1)=>-1/24, (0,0)=>25/24))

get_fvtofd_i(mesh) = Iip(mesh) + Iim(mesh) - Ii(mesh)
get_fvtofd_j(mesh) = Ijp(mesh) + Ijm(mesh) - Ij(mesh)

get_fvtofd_ij(mesh) = get_fvtofd_j(mesh) * get_fvtofd_i(mesh)

#Differential TODO figure out non periodic
get_diff_i(mesh; factor = nothing) = get_op_from_coeffs(mesh, Dict((-1,0) => -1)) + get_op_from_coeffs(mesh, Dict((0,0) => 1))
get_diff_j(mesh; factor = nothing) = get_op_from_coeffs(mesh, Dict((0,-1) => -1)) + get_op_from_coeffs(mesh, Dict((0,0) => 1))
#get_diff_j(mesh; factor = nothing) = get_op_from_coeffs(mesh, Dict((0,-1) => -1, (0,0) => 1))

#Codifferential TODO idem
get_codiff_i(mesh; factor = ones(mesh.ni, mesh.nj)) = get_op_from_coeffs(mesh, Dict((0,0)=>-1, (1,0)=>1); factor)
get_codiff_j(mesh; factor = ones(mesh.ni, mesh.nj)) = get_op_from_coeffs(mesh, Dict((0,0)=>-1, (0,1)=>1); factor)

function get_laplacian(mesh; o=2)
	if mesh.iperio & mesh.jperio
		return (get_codiff_i(mesh; factor=mesh.dy ./mesh.dx) * get_fvtofd_i(mesh)) * (get_diff_i(mesh) * get_fvtofd_i(mesh)) .* (mesh.dy/mesh.dx) + (get_codiff_j(mesh; factor=mesh.dx ./mesh.dy) * get_fvtofd_j(mesh)) * (get_diff_j(mesh) * get_fvtofd_j(mesh)) .* (mesh.dx/mesh.dy) #WORKS FOR PERIODIC
	elseif o==4 
		return (get_diff_i(mesh; factor=mesh.dy ./mesh.dx)' * get_fvtofd_i(mesh)) * (get_diff_i(mesh) * get_fvtofd_i(mesh)) + (get_diff_j(mesh; factor = mesh.dx./mesh.dy)' * get_fvtofd_j(mesh)) * (get_diff_j(mesh) * get_fvtofd_j(mesh)) 
		#WORKS FOR NON PERIODIC but is of order 0
	else
		return get_diff_i(mesh; factor=mesh.dy ./ mesh.dx)' * get_diff_i(mesh) + get_diff_j(mesh; factor = mesh.dx ./mesh.dy)' * get_diff_j(mesh) 
		#WORKS FOR NON PERIODIC (order 2) 
	end
end

laplacian(mesh, msk, bc, location; o) = get_laplacian(mesh; o)

"""
	solve_poisson(out, b)

Solves the poisson problem Ax = b and sets out to be x
"""
function get_poisson_solver(mesh, b, bc, form_type; o = 2)
	@assert bc in ["dirichlet", "neumann"]
	@assert form_type in ["0p", "0d", "2p", "2d"]

	if form_type == "0p" #TODO check
		location = "vertex"
		msk = mesh.msk0p
	elseif form_type == "0d"
		location = "center"
		msk = mesh.msk0d
	elseif form_type == "2p"
		location = "center"
		msk = mesh.msk2p
	elseif form_type == "2d"
		location = "vertex"
		msk = mesh.msk2d
	end

	A,idx = desingularize_op(laplacian(mesh, msk, bc, location; o))
	display(Matrix(A))
	A = -A#factorize(A)#lu(-A)
	prob = LinearProblem(A, b[idx])
	linsolve = init(prob, KrylovJL_CRAIGMR())
	#=
	if bc=="dirichlet"
		if location == "vertex"
			#Only for dirichlet BC at vertices
			A[A .== -2] .= -4
			A[A .== -3] .= -4
		elseif location == "center"
			#Dirichlet centers
			A[A .== -2] .= -6
			A[A .== -3] .= -5
		end
	end
	=#
	A = lu(A)
	#Create a function that solves Ax = b
	function func!(out, b)
		#ap_inv!(linsolve, b, mesh; out, idx)
		#=linsolve.b = b[idx]
		sol = solve!(linsolve)
		out[idx] .= sol.u
		=#
		@views out[idx] .= A \ b[idx]
		
		#println(mean(out[idx]))
		#out[idx] .-= mean(out[idx])
	end

	return func!
end
