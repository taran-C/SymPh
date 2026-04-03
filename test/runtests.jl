using SymPh
using SymPh.Maths
import SymPh.Arrays
using LoopManagers: PlainCPU, VectorizedCPU, MultiThread

using Test

#Addition
function test_add_0p(mesh)
	#Equation
	@Let a = FormVariable{0, Primal}()
	@Let b = FormVariable{0, Primal}()

	@Let sum0p = a+b
	
	comp! = to_kernel(sum0p)

	state = State(mesh)

	state.a .= rand()*mesh.msk0p
	state.b .= rand()*mesh.msk0p

	@invokelatest comp!(mesh, state)

	@info "Addition of two Primal 0-forms"
	return state.sum0p ≈ (state.a .+ state.b)
end

function test_add_1p(mesh)
	#Equation
	@Let a = FormVariable{1, Primal}()
	@Let b = FormVariable{1, Primal}()

	@Let sum1p = a+b
	
	comp! = to_kernel(sum1p)

	state = State(mesh)

	state.a_i .= rand()*mesh.msk1pi
	state.a_j .= rand()*mesh.msk1pj
	state.b_i .= rand()*mesh.msk1pi
	state.b_j .= rand()*mesh.msk1pj

	@invokelatest comp!(mesh, state)

	@info "Addition of two Primal 1-forms"
	return (state.sum1p_i ≈ (state.a_i .+ state.b_i)) & (state.sum1p_j ≈ (state.a_j .+ state.b_j))
end

function test_add_2p(mesh)
	#Equation
	@Let a = FormVariable{2, Primal}()
	@Let b = FormVariable{2, Primal}()

	@Let sum2p = a+b
	
	comp! = to_kernel(sum2p)

	state = State(mesh)

	state.a .= rand()*mesh.msk2p
	state.b .= rand()*mesh.msk2p

	@invokelatest comp!(mesh, state)

	@info "Addition of two Primal 2-forms"
	return state.sum2p ≈ (state.a .+ state.b)
end

function test_add_0d(mesh)
	#Equation
	@Let a = FormVariable{0, Dual}()
	@Let b = FormVariable{0, Dual}()

	@Let sum0d = a+b
	
	comp! = to_kernel(sum0d)

	state = State(mesh)

	state.a .= rand()*mesh.msk0d
	state.b .= rand()*mesh.msk0d

	@invokelatest comp!(mesh, state)

	@info "Addition of two Dual 0-forms"
	return state.sum0d ≈ (state.a .+ state.b)
end

function test_add_1d(mesh)
	#Equation
	@Let a = FormVariable{1, Dual}()
	@Let b = FormVariable{1, Dual}()

	@Let sum1d = a+b
	
	comp! = to_kernel(sum1d)

	state = State(mesh)

	state.a_i .= rand()*mesh.msk1di
	state.a_j .= rand()*mesh.msk1dj
	state.b_i .= rand()*mesh.msk1di
	state.b_j .= rand()*mesh.msk1dj

	@invokelatest comp!(mesh, state)

	@info "Addition of two Dual 1-forms"
	return (state.sum1d_i ≈ (state.a_i .+ state.b_i)) & (state.sum1d_j ≈ (state.a_j .+ state.b_j))
end

function test_add_2d(mesh)
	#Equation
	@Let a = FormVariable{2, Dual}()
	@Let b = FormVariable{2, Dual}()

	@Let sum2d = a+b
	
	comp! = to_kernel(sum2d)

	state = State(mesh)

	state.a .= rand()*mesh.msk2d
	state.b .= rand()*mesh.msk2d

	@invokelatest comp!(mesh, state)

	@info "Addition of two Dual 2-forms"
	return state.sum2d ≈ (state.a .+ state.b)
end

@testset "SymPh" begin
	#mesh
	nh=3
	ni=20
	nj=15
	simd = VectorizedCPU(16)
	
	mesh=Arrays.CartesianMesh(ni,nj,nh,simd)

	#tests
	@test test_add_0p(mesh)
	@test test_add_1p(mesh)
	@test test_add_2p(mesh)
	@test test_add_0d(mesh)
	@test test_add_1d(mesh)
	@test test_add_2d(mesh)
end
