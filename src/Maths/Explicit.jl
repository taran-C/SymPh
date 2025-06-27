import ..Arrays
export explicit, ExplicitParam

"""
	ExplicitParam(;interp = Arrays.upwind, fvtofd = Arrays.fvtofd2, fdtofv = Arrays.fdtofv2)
	
The parameters used when necessary to choose numerical methods
"""
struct ExplicitParam
	interp
	fvtofd
	fdtofv
	laporder

	function ExplicitParam(;interp = Arrays.upwind, fvtofd = Arrays.fvtofd2, fdtofv = Arrays.fdtofv2, laporder = 2)
		return new(interp, fvtofd, fdtofv, laporder)
	end
end


#------------------------Vectors----------------------------------------------------------------

"""
	explicit(expr::Expression, param = ExplicitParam())

Returns a `SymPh.Arrays` `Expression` representing the operation `expr` on an array, pulling numerical methods from `param`
"""
function explicit(vect::VectorVariable{P}; param = ExplicitParam()) where {P}
	return [Arrays.ArrayVariable(vect.name*"_X"), Arrays.ArrayVariable(vect.name*"_Y")]
end

#----------------------FormVariables---------------------------------------------------------------

function explicit(form::FormVariable{0, P}; param = ExplicitParam()) where {P}
	return Arrays.ArrayVariable(form.name)
end

function explicit(form::FormVariable{1, P}; param = ExplicitParam()) where {P}
	return [Arrays.ArrayVariable(form.name*"_i"), Arrays.ArrayVariable(form.name*"_j")]
end

function explicit(form::FormVariable{2,P}; param = ExplicitParam()) where {P}
	return Arrays.ArrayVariable(form.name)
end

#----------------------FuncCall-----------------------------------------------------------------------

function explicit(form::FuncCall{0,P}; param = ExplicitParam()) where {P}
	#Expliciting arguments
	argexprs = []
	for arg in form.args
		push!(argexprs, explicit(arg; param))
	end

	call = Arrays.FuncCall(form.name, form.func, argexprs, 0, 0)
	display(call)
	return call
end
#TODO not great, find a way to only call one function
function explicit(form::FuncCall{1,P}; param = ExplicitParam()) where {P}
	#Expliciting arguments
	argexprs = []
	for arg in form.args
		push!(argexprs, explicit(arg; param))
	end

	call_i = Arrays.FuncCall(form.name*"_i", form.func[1], argexprs, 0, 0)
	call_j = Arrays.FuncCall(form.name*"_j", form.func[2], argexprs, 0, 0)
	return (call_i, call_j)
end
function explicit(form::FuncCall{2,P}; param = ExplicitParam()) where {P}
	#Expliciting arguments
	argexprs = []
	for arg in form.args
		push!(argexprs, explicit(arg; param))
	end

	call = Arrays.FuncCall(form.name, form.func, argexprs, 0, 0)
	display(call)
	return call
end

#-----------------------------Addition-----------------------------------------------------------------------------
function explicit(form::Addition{0,P}; param = ExplicitParam()) where {P}
	return Arrays.Addition(form.name, explicit(form.left; param = param), explicit(form.right; param = param))
end

function explicit(form::Addition{1,P}; param = ExplicitParam()) where {P}
	ls = explicit(form.left; param = param)
	rs = explicit(form.right; param = param)

	return [Arrays.Addition(form.name*"_i", ls[1], rs[1]), Arrays.Addition(form.name*"_j", ls[2], rs[2])]
end

function explicit(form::Addition{2,P}; param = ExplicitParam()) where {P}
	return Arrays.Addition(form.name, explicit(form.left; param = param), explicit(form.right; param = param))
end

#-----------------------------Substraction---------------------------------------------------------------------------
function explicit(form::Substraction{0,P}; param = ExplicitParam()) where {P}
	return Arrays.Addition(form.name, explicit(form.left; param = param), explicit(form.right; param = param))
end

function explicit(form::Substraction{1,P}; param = ExplicitParam()) where {P}
	ls = explicit(form.left; param = param)
	rs = explicit(form.right; param = param)

	return [Arrays.Substraction(form.name*"_i", ls[1], rs[1]), Arrays.Substraction(form.name*"_j", ls[2], rs[2])]
end

function explicit(form::Substraction{2,P}; param = ExplicitParam()) where {P}
	return Arrays.Substraction(form.name, explicit(form.left; param = param), explicit(form.right; param = param))
end

#----------------------Negative---------------------------------------------------------------------------------------
function explicit(form::Negative{0,P}; param = ExplicitParam()) where {P}
	return Arrays.Negative(form.name, explicit(form.form; param = param))
end

function explicit(form::Negative{1,P}; param = ExplicitParam()) where {P}
	fexpr = explicit(form.form; param = param)
	return [Arrays.Negative(form.name*"_i", fexpr[1]), Arrays.Negative(form.name*"_j", fexpr[2])]
end

function explicit(form::Negative{2,P}; param = ExplicitParam()) where {P}
	return Arrays.Negative(form.name, explicit(form.form; param = param))
end

#------------------------Division------------------------------------------------------------------------
function explicit(form::Division{2,Dual}; param = ExplicitParam())
	left = explicit(form.left; param = param)
	right = explicit(form.right; param = param)

	res = left / Arrays.avg4pt(right, -1, -1) * Arrays.msk2d
	res.name = form.name

	return res
end

#----------------------RealProducts---------------------------------------------------------------------
function explicit(form::RealProdForm{0,D}; param = ExplicitParam()) where {D}
	res = form.real * explicit(form.form; param = param)
	res.name = form.name

	return res
end

function explicit(form::RealProdForm{1,D}; param = ExplicitParam()) where {D}
	exprs = explicit(form.form; param = param)

	l = form.real * exprs[1]
	r = form.real * exprs[2]

	l.name = form.name*"_i"
	r.name = form.name*"_j"

	return [l, r]
end

function explicit(form::RealProdForm{2,D}; param = ExplicitParam()) where {D}
	res = form.real * explicit(form.form; param = param)
	res.name = form.name

	return res
end

#--------------------------------------Wedge--------------------------------------------------------------------
function explicit(form::Wedge{1, 0, 1, Dual}; param = ExplicitParam())
	left = explicit(form.left; param = param)
	righti, rightj = explicit(form.right; param = param)

	if true #form.interp == nothing
		interp = param.interp
	else
		interp = form.interp
	end

	#TODO use interp here !
	li = 0.5 * (left[0,0] + left[-1,0])
	lj = 0.5 * (left[0,0] + left[0,-1])
	#li = interp(righti, left, Arrays.o2di, "left", "i")
	#lj = interp(rightj, left, Arrays.o2dj, "left", "j")
	
	resi = li * righti * Arrays.msk1di
	resj = lj * rightj * Arrays.msk1dj
	
	resi.name = form.name*"_i"
	resj.name = form.name*"_j"
	
	return [resi, resj]
end

#------------------------------ExteriorDerivative--------------------------------------------------------------
function explicit(form::ExteriorDerivative{1, Primal}; param = ExplicitParam())
	expr = explicit(form.form; param = param)

	#Actual differentiation
	d_i = (expr[1,0] - expr[0,0]) * Arrays.msk1pi
	d_j = (expr[0,1] - expr[0,0]) * Arrays.msk1pj
	
	#Renaming
	d_i.name = form.name * "_i"
	d_j.name = form.name * "_j"

	return [d_i, d_j]
end

function explicit(form::ExteriorDerivative{1, Dual}; param = ExplicitParam())
	expr = explicit(form.form; param = param)
	
	#Actual differentiation
	d_i = (expr[0,0] - expr[-1,0]) * Arrays.msk1di
	d_j = (expr[0,0] - expr[0,-1]) * Arrays.msk1dj
	
	#Renaming
	d_i.name = form.name * "_i"
	d_j.name = form.name * "_j"

	return [d_i, d_j]
end

function explicit(form::ExteriorDerivative{2, Primal}; param = ExplicitParam())
	exprs = explicit(form.form; param = param)

	#Actual differentiation
	dq = ((exprs[2][1,0]-exprs[2][0,0])-(exprs[1][0,1]-exprs[1][0,0])) * Arrays.msk2p
	
	#Renaming
	dq.name = form.name

	return dq
end

function explicit(form::ExteriorDerivative{2, Dual}; param = ExplicitParam())
	exprs = explicit(form.form; param = param)
	
	#Actual differentiation
	dq = ((exprs[2][0,0]-exprs[2][-1,0])-(exprs[1][0,0]-exprs[1][0,-1])) * Arrays.msk2d
	
	#Renaming
	dq.name = form.name

	return dq
end

#---------------------------------Codifferential-----------------------------------------------------
function explicit(form::Codifferential{0, Dual}; param = ExplicitParam())
	exprs = explicit(form.form; param = param)

	dq = ((exprs[1][1,0]-exprs[1][0,0]) + (exprs[2][0,1]-exprs[2][0,0])) * Arrays.msk0d
	dq.name = form.name

	return dq
end

function explicit(form::Codifferential{1, Dual}; param = ExplicitParam())
	expr = explicit(form.form; param = param)

	du = -(expr[0,1] - expr[0,0]) * Arrays.msk1di
	dv = (expr[1,0] - expr[0,0]) * Arrays.msk1dj

	du.name = form.name * "_i"
	dv.name = form.name * "_j"

	return[du, dv]
end

#----------------------------------InteriorProduct---------------------------------------------------------

function explicit(form::InteriorProduct{0, Dual, Dual}; param = ExplicitParam())
	fx, fy = explicit(form.form; param = param)
	U, V = explicit(form.vect; param = param)

	fu = fx * U
	fv = fy * V
	
	if form.interp == nothing
		interp = param.interp
	else
		interp = form.interp
	end
	
	Uint = interp(U[0,0], fu, Arrays.o1di, "right", "i")
	Vint = interp(V[0,0], fv, Arrays.o1dj, "right", "j")
	
	qout = (Uint + Vint) * Arrays.msk0d
	
	qout.name = form.name

	return qout
end

function explicit(form::InteriorProduct{1, Dual, Primal}; param = ExplicitParam()) 
	fexpr = explicit(form.form; param = param)
	uexpr, vexpr = explicit(form.vect; param = param)

	if form.interp == nothing
		interp = param.interp
	else
		interp = form.interp
	end

	uout = -vexpr * interp(vexpr, fexpr, Arrays.o2pj, "left", "j") * Arrays.msk1pi
	vout = uexpr * interp(uexpr, fexpr, Arrays.o2pi, "left", "i") * Arrays.msk1pj

	uout.name = form.name*"_i"
	vout.name = form.name*"_j"

	return [uout, vout]
end

#TODO better interp here to have 4th order rsw
function explicit(form::InteriorProduct{1, Dual, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)
	uexpr, vexpr = explicit(form.vect; param = param)

	if form.interp == nothing
		interp = param.interp
	else
		interp = form.interp
	end

	#udec = Arrays.avg4pt(uexpr, 1, -1)
	#vdec = Arrays.avg4pt(vexpr, -1, 1)

	udec = Arrays.centered4(vexpr, Arrays.centered4(uexpr, uexpr, Arrays.o1di, "right", "i"), Arrays.o1dj, "left", "j")
	vdec = Arrays.centered4(vexpr, Arrays.centered4(uexpr, vexpr, Arrays.o1di, "left", "i"), Arrays.o1dj, "right", "j")

	#TODO transp velocity dec or not
	iout = -vdec * interp(vdec, fexpr, Arrays.o2dj, "right", "j") * Arrays.msk1di
	jout = udec * interp(udec, fexpr, Arrays.o2di, "right", "i") * Arrays.msk1dj
	
	iout.name = form.name*"_i"
	jout.name = form.name*"_j"

	return [iout, jout]
end

#---------------------------------Sharp---------------------------------------------------------------------

function explicit(vec::Sharp{D}; param = ExplicitParam()) where D #TODO separate Primal and dual areas (could be very different, especially for non square grids)
	iexpr, jexpr = explicit(vec.form; param = param)
	if vec.fvtofd == nothing
		fvtofd = param.fvtofd
	else
		fvtofd = vec.fvtofd
	end
	if vec.fdtofv == nothing
		fdtofv = param.fdtofv
	else
		fdtofv = vec.fdtofv
	end

	#TODO per object configurable fvtofd function
	xout = fdtofv(fvtofd(iexpr, Arrays.msk1di, "i"), Arrays.msk1di, "j") / Arrays.dx / Arrays.dx * Arrays.msk1di 
	yout = fdtofv(fvtofd(jexpr, Arrays.msk1dj, "j"), Arrays.msk1dj, "i") / Arrays.dy / Arrays.dy * Arrays.msk1dj

	xout.name = vec.name*"_X"
	yout.name = vec.name*"_Y"
	return [xout, yout]
end

#--------------------------------Hodge-------------------------------------------------------------------

#UNTESTED
function explicit(form::Hodge{0, Primal}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	res = param.fvtofd(param.fvtofd(fexpr, Arrays.msk2d, "i"), Arrays.msk2d, "j") / Arrays.dx / Arrays.dy * Arrays.msk0p
	
	res.name = form.name
	return res
end

function explicit(form::Hodge{0, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	res = param.fvtofd(param.fvtofd(fexpr, Arrays.msk2p, "i"), Arrays.msk2p, "j") / Arrays.dx / Arrays.dy * Arrays.msk0d
	
	res.name = form.name
	return res
end

function explicit(form::Hodge{2, Primal}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	res = param.fdtofv(param.fdtofv(fexpr, Arrays.msk0d, "i"), Arrays.msk0d, "j") * Arrays.dx * Arrays.dy * Arrays.msk2p
	
	res.name = form.name
	return res
end

#UNTESTED
function explicit(form::Hodge{2, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	res = param.fdtofv(param.fdtofv(fexpr, Arrays.msk0p, "i"), Arrays.msk0p, "j") * Arrays.dx * Arrays.dy * Arrays.msk2d
	
	res.name = form.name
	return res
end

#----------------------------------InverseLaplacian---------------------------------------------------------------
#TODO have bcs influence that

function explicit(form::InverseLaplacian{0, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	poisson = Poisson2D("dirichlet", "0d"; order = param.laporder)

	function poiss_dirich_0d(mesh;kwargs...)
		#args = Dict(kwargs)
		solve_poisson(poisson, mesh, kwargs[Symbol(form.name)], kwargs[Symbol(form.form.name)])
	end

	return Arrays.FuncCall(form.name, poiss_dirich_0d, [fexpr], 0, 0)
end

function explicit(form::InverseLaplacian{2, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)
	
	poisson = Poisson2D("dirichlet", "2d"; order = param.laporder)

	function poiss_dirich_2d(mesh;kwargs...)
		#args = Dict(kwargs)
		solve_poisson(poisson, mesh, kwargs[Symbol(form.name)], kwargs[Symbol(form.form.name)])
	end

	return Arrays.FuncCall(form.name, poiss_dirich_2d, [fexpr], 0, 0)
end
