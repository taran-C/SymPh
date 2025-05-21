import ..Arrays
export explicit, ExplicitParam

struct ExplicitParam
	interp
	fvtofd

	function ExplicitParam(;interp = Arrays.upwind, fvtofd = Arrays.fvtofd2)
		return new(interp, fvtofd)
	end
end


#Vectors
function explicit(vect::VectorVariable{P}; param = ExplicitParam()) where {P}
	return [Arrays.ArrayVariable(vect.name*"_X"), Arrays.ArrayVariable(vect.name*"_Y")]
end

#FormVariables
function explicit(form::FormVariable{0, P}; param = ExplicitParam()) where {P}
	return Arrays.ArrayVariable(form.name)
end

function explicit(form::FormVariable{1, P}; param = ExplicitParam()) where {P}
	return [Arrays.ArrayVariable(form.name*"_i"), Arrays.ArrayVariable(form.name*"_j")]
end

function explicit(form::FormVariable{2,P}; param = ExplicitParam()) where {P}
	return Arrays.ArrayVariable(form.name)
end

#Addition
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

#Substraction
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

#Negative
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

#Division
function explicit(form::Division{2,Dual}; param = ExplicitParam())
	left = explicit(form.left; param = param)
	right = explicit(form.right; param = param)

	res = left / Arrays.avg4pt(right, -1, -1) * Arrays.msk2d
	res.name = form.name

	return res
end

#RealProducts
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

#ExteriorDerivative
function explicit(form::ExteriorDerivative{1, Primal}; param = ExplicitParam())
	expr = explicit(form.form; param = param)

	#Finite diff in both direction
	expi = param.fvtofd(expr, Arrays.msk0p, "i")
	expj = param.fvtofd(expr, Arrays.msk0p, "j")

	#Actual differentiation
	d_i = (expi[1,0] - expi[0,0]) * Arrays.msk1pi
	d_j = (expj[0,1] - expj[0,0]) * Arrays.msk1pj
	
	#Renaming
	d_i.name = form.name * "_i"
	d_j.name = form.name * "_j"

	return [d_i, d_j]
end

function explicit(form::ExteriorDerivative{1, Dual}; param = ExplicitParam())
	expr = explicit(form.form; param = param)
	
	#Finite diff in both direction
	expi = param.fvtofd(expr, Arrays.msk0d, "i")
	expj = param.fvtofd(expr, Arrays.msk0d, "j")

	#Actual differentiation
	d_x = (expi[0,0] - expi[-1,0]) * Arrays.msk1di
	d_y = (expj[0,0] - expj[0,-1]) * Arrays.msk1dj
	
	#Renaming
	d_x.name = form.name * "_i"
	d_y.name = form.name * "_j"

	return [d_x, d_y]
end

function explicit(form::ExteriorDerivative{2, Primal}; param = ExplicitParam())
	exprs = explicit(form.form; param = param)

	#Finite diff in both direction
	expi = param.fvtofd(exprs[1], Arrays.msk1di, "j")
	expj = param.fvtofd(exprs[2], Arrays.msk1dj, "i")

	#Actual differentiation
	dq = ((expj[1,0]-expj[0,0])-(expi[0,1]-expi[0,0])) * Arrays.msk2p
	
	#Renaming
	dq.name = form.name

	return dq
end

function explicit(form::ExteriorDerivative{2, Dual}; param = ExplicitParam())
	exprs = explicit(form.form; param = param)
	
	#Finite diff in both direction
	expi = param.fvtofd(exprs[1], Arrays.msk1di, "j")
	expj = param.fvtofd(exprs[2], Arrays.msk1dj, "i")

	#Actual differentiation
	dq = ((expj[0,0]-expj[-1,0])-(expi[0,0]-expj[0,-1])) * Arrays.msk2d
	
	#Renaming
	dq.name = form.name

	return dq
end

#Codifferential #TODO handle FV to FD for codifferential
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

#InteriorProduct
function explicit(form::InteriorProduct{0, Dual, Dual}; param = ExplicitParam())
	fx, fy = explicit(form.form; param = param)
	U, V = explicit(form.vect; param = param)

	fu = fx * U
	fv = fy * V
	
	if form.interp == Nothing
		interp = param.interp
	else
		interp = form.interp
	end
	
	Uint = interp(U[0,0] + U[1,0], fu, Arrays.o1di, "right", "i")
	Vint = interp(V[0,0] + V[0,1], fv, Arrays.o1dj, "right", "j")
	
	qout = (Uint + Vint) * Arrays.msk0d
	
	qout.name = form.name

	return qout
end

function explicit(form::InteriorProduct{1, Dual, Primal}; param = ExplicitParam()) 
	fexpr = explicit(form.form; param = param)
	uexpr, vexpr = explicit(form.vect; param = param)

	if form.interp == Nothing
		interp = param.interp
	else
		interp = form.interp
	end

	finti = interp(uexpr, fexpr, Arrays.o1pi, "left", "i")
	fintj = interp(vexpr, fexpr, Arrays.o1pj, "left", "j")

	uout = -vexpr * fintj * Arrays.msk1pi
	vout = uexpr * finti * Arrays.msk1pj

	uout.name = form.name*"_i"
	vout.name = form.name*"_j"

	return [uout, vout]
end

function explicit(form::InteriorProduct{1, Dual, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)
	uexpr, vexpr = explicit(form.vect; param = param)

	if form.interp == Nothing
		interp = param.interp
	else
		interp = form.interp
	end

	udec = Arrays.avg4pt(uexpr, 1, -1)
	vdec = Arrays.avg4pt(vexpr, -1, 1)

	#udec = interp(vexpr, interp(uexpr, uexpr, Arrays.o1di, "right", "x"), Arrays.o1dj, "left", "y")
	#vdec = interp(vexpr, interp(uexpr, vexpr, Arrays.o1di, "left", "x"), Arrays.o1dj, "right", "y")

	#TODO transp velocity dec or not
	iout = -vdec * interp(vdec, fexpr, Arrays.o2dj, "right", "i") * Arrays.msk1di
	jout = udec * interp(udec, fexpr, Arrays.o2di, "right", "j") * Arrays.msk1dj
	
	iout.name = form.name*"_i"
	jout.name = form.name*"_j"

	return [iout, jout]
end

#Sharp
function explicit(vec::Sharp{D}; param = ExplicitParam()) where D #TODO separate Primal and dual areas (could be very different, especially for non square grids)
	iexpr, jexpr = explicit(vec.form; param = param)

	#TODO per object configurable fvtofd function
	#xout = param.fvtofd(iexpr, Arrays.msk1di, "i") / Arrays.dx / Arrays.dx * Arrays.msk1di 
	#yout = param.fvtofd(jexpr, Arrays.msk1dj, "j") / Arrays.dy / Arrays.dy * Arrays.msk1dj

	xout = iexpr / Arrays.dx / Arrays.dx
	yout = jexpr / Arrays.dy / Arrays.dy

	xout.name = vec.name*"_X"
	yout.name = vec.name*"_Y"
	return [xout, yout]
end

#Hodge
function explicit(form::Hodge{0, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	#res = param.fvtofd(param.fvtofd(fexpr, Arrays.msk2p, "i"), Arrays.msk2p, "j") / Arrays.dx /Arrays.dy * Arrays.msk0d
	res = fexpr / Arrays.dx / Arrays.dy
	res.name = form.name
	return res
end

#InverseLaplacian
function explicit(form::InverseLaplacian{0, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)

	poisson = Poisson2D("dirichlet", "0d")

	function poiss_dirich_0d(mesh;kwargs...)
		args = Dict(kwargs)
		solve_poisson(poisson, mesh, kwargs[Symbol(form.name)], kwargs[Symbol(form.form.name)])
	end

	return Arrays.FuncCall(form.name, poiss_dirich_0d, [fexpr], 0, 0)
end

function explicit(form::InverseLaplacian{2, Dual}; param = ExplicitParam())
	fexpr = explicit(form.form; param = param)
	#TODO check with periodic condition +bc
	poisson = Poisson2D("dirichlet", "2d")

	function poiss_dirich_2d(mesh;kwargs...)
		args = Dict(kwargs)
		solve_poisson(poisson, mesh, kwargs[Symbol(form.name)], kwargs[Symbol(form.form.name)])
	end

	return Arrays.FuncCall(form.name, poiss_dirich_2d, [fexpr], 0, 0)
end
