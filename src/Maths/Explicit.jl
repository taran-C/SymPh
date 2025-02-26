import ..Arrays
export explicit

#TODO change this to a param (or something) object
interp = "upwind"

#Vectors
function explicit(vect::VectorVariable{P}) where {P}
	return [Arrays.ArrayVariable(vect.name*"_X"), Arrays.ArrayVariable(vect.name*"_Y")]
end

#FormVariables
function explicit(form::FormVariable{0, P}) where {P}
	return Arrays.ArrayVariable(form.name)
end

function explicit(form::FormVariable{1, P}) where {P}
	return [Arrays.ArrayVariable(form.name*"_x"), Arrays.ArrayVariable(form.name*"_y")]
end

function explicit(form::FormVariable{2,P}) where {P}
	return Arrays.ArrayVariable(form.name)
end

#Addition
function explicit(form::Addition{0,P}) where {P}
	return Arrays.Addition(form.name, explicit(form.left), explicit(form.right))
end

function explicit(form::Addition{1,P}) where {P}
	ls = explicit(form.left)
	rs = explicit(form.right)

	return [Arrays.Addition(form.name*"_x", ls[1], rs[1]), Arrays.Addition(form.name*"_y", ls[2], rs[2])]
end

function explicit(form::Addition{2,P}) where {P}
	return Arrays.Addition(form.name, explicit(form.left), explicit(form.right))
end

#Substraction
function explicit(form::Substraction{0,P}) where {P}
	return Arrays.Addition(form.name, explicit(form.left), explicit(form.right))
end

function explicit(form::Substraction{1,P}) where {P}
	ls = explicit(form.left)
	rs = explicit(form.right)

	return [Arrays.Substraction(form.name*"_x", ls[1], rs[1]), Arrays.Substraction(form.name*"_y", ls[2], rs[2])]
end

function explicit(form::Substraction{2,P}) where {P}
	return Arrays.Substraction(form.name, explicit(form.left), explicit(form.right))
end

#Negative
function explicit(form::Negative{0,P}) where {P}
	return Arrays.Negative(form.name, explicit(form.form))
end

function explicit(form::Negative{1,P}) where {P}
	fexpr = explicit(form.form)
	return [Arrays.Negative(form.name*"_x", fexpr[1]), Arrays.Negative(form.name*"_y", fexpr[2])]
end

function explicit(form::Negative{2,P}) where {P}
	return Arrays.Negative(form.name, explicit(form.form))
end

#RealProducts
function explicit(form::RealProdForm{0,D}) where {D}
	return form.real * explicit(form.form)
end

function explicit(form::RealProdForm{1,D}) where {D}
	exprs = explicit(form.form)

	return [form.real * exprs[1], form.real * exprs[2]]
end

function explicit(form::RealProdForm{2,D}) where {D}
	return form.real * explicit(form.form)
end

#ExteriorDerivative
function explicit(form::ExteriorDerivative{1, Primal})
	expr = explicit(form.form)
	d_x = (expr[1,0] - expr[0,0]) * Arrays.msk1px
	d_x.name = form.name * "_x"

	d_y = (expr[0,1] - expr[0,0]) * Arrays.msk1py
	d_y.name = form.name * "_y"

	return [d_x, d_y]
end

function explicit(form::ExteriorDerivative{1, Dual})
	expr = explicit(form.form)
	d_x = (expr[0,0] - expr[-1,0]) * Arrays.msk1dx
	d_x.name = form.name * "_x"

	d_y = (expr[0,0] - expr[0,-1]) * Arrays.msk1dy
	d_y.name = form.name * "_y"

	return [d_x, d_y]
end

function explicit(form::ExteriorDerivative{2, Primal})
	exprs = explicit(form.form)
	dq = ((exprs[2][1,0]-exprs[2][0,0])-(exprs[1][0,1]-exprs[1][0,0])) * Arrays.msk2p
	dq.name = form.name

	return dq
end

function explicit(form::ExteriorDerivative{2, Dual})
	exprs = explicit(form.form)
	dq = ((exprs[2][0,0]-exprs[2][-1,0])-(exprs[1][0,0]-exprs[1][0,-1])) * Arrays.msk2d
	dq.name = form.name

	return dq
end

#InteriorProduct
#TODO Not tested !!
function explicit(form::InteriorProduct{0, Dual, Dual})
	fx, fy = explicit(form.form)
	U, V = explicit(form.vect)

	fu = fx * uexpr
	fv = fy * vexpr
	
	if interp == "upwind"
		Uint = Arrays.upwind(U[0,0] + U[1,0], fu, Arrays.o2px, "left", "x")
		Vint = Arrays.upwind(V[0,0] + V[0,1], fv, Arrays.o2py, "left", "y")
	elseif interp == "2ptavg"
		Uint = 0.5 * (fu[0,0] + fu[-1,0])
		Vint = 0.5 * (fv[0,0] + fv[0,-1])
	end
	
	qout = (Uint + Vint) * Arrays.msk0d

	return 
end

function explicit(form::InteriorProduct{1, Dual, Primal}) 
	fexpr = explicit(form.form)
	uexpr, vexpr = explicit(form.vect)


	if interp == "upwind"
		fintx = Arrays.upwind(uexpr, fexpr, Arrays.o2px, "left", "x")
		finty = Arrays.upwind(vexpr, fexpr, Arrays.o2py, "left", "y")
	elseif interp == "2ptavg"
		fintx = 0.5 * (fexpr[0,0] + fexpr[-1,0])
		finty = 0.5 * (fexpr[0,0] + fexpr[0,-1])
	end

	uout = -vexpr * finty * Arrays.msk1px
	vout = uexpr * fintx * Arrays.msk1py

	uout.name = form.name*"_x"
	vout.name = form.name*"_y"

	return [uout, vout]
end

function explicit(form::InteriorProduct{1, Dual, Dual})
	fexpr = explicit(form.form)
	uexpr, vexpr = explicit(form.vect)

	udec = Arrays.avg4pt(uexpr, 1, -1)
	vdec = Arrays.avg4pt(vexpr, -1, 1)
	
	if interp == "upwind"
		xout = -vdec * Arrays.upwind(vexpr, fexpr, Arrays.o2dy, "right", "y") * Arrays.msk1dx
		yout = udec * Arrays.upwind(uexpr, fexpr, Arrays.o2dx, "right", "x") * Arrays.msk1dy
	elseif interp == "2ptavg"
		xout = -vdec * 0.5 * (fexpr[0,0]+fexpr[0,1]) * Arrays.msk1dx
		yout = udec * 0.5 * (fexpr[0,0]+fexpr[1,0]) * Arrays.msk1dy
	else
		@assert false "TODO"
	end

	return [xout, yout]
end

#Sharp
function explicit(vec::Sharp{D}) where D #TODO separate Primal and dual areas (could be very different, especially for non square grids)
	xexpr, yexpr = explicit(vec.form)

	return [xexpr/Arrays.dx, yexpr/Arrays.dy]
end

#Hodge
function explicit(form::Hodge{0, Dual})
	fexpr = explicit(form.form)

	return fexpr / Arrays.A
end

#InnerProduct
function explicit(form::InnerProduct{2, Primal})
	ax, ay = explicit(form.left)
	bx, by = explicit(form.right)
	return 0.5*(ax*bx + ax[1,0]*bx[1,0] + ay*by + ay[0,1]*by[0,1]) * Arrays.msk2p
end
