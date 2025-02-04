import ..Arrays
export explicit

#Vectors
function explicit(vect::VectorVariable{P}) where {P}
	return (Arrays.ArrayVariable(vect.name*"_X"), Arrays.ArrayVariable(vect.name*"_Y"))
end

#0-Forms
function explicit(form::FormVariable{0, P}) where {P}
	return Arrays.ArrayVariable(form.name)
end

function explicit(form::Addition{0,P}) where {P}
	return Arrays.Addition(form.name, explicit(form.left), explicit(form.right))
end

function explicit(form::ExteriorDerivative{1, Primal})
	expr = explicit(form.form)
	d_x = (expr[1,0] - expr[0,0]) * Arrays.mskx
	d_x.name = form.name * "_x"

	d_y = (expr[0,1] - expr[0,0]) * Arrays.msky
	d_y.name = form.name * "_y"

	return (d_x, d_y)
end

#1-Forms
function explicit(form::FormVariable{1, P}) where {P}
	return (Arrays.ArrayVariable(form.name*"_x"), Arrays.ArrayVariable(form.name*"_y"))
end

function explicit(form::Addition{1,P}) where {P}
	ls = explicit(form.left)
	rs = explicit(form.right)

	return (Arrays.Addition(form.name*"_x", ls[1], rs[1]), Arrays.Addition(form.name*"_y", ls[2], rs[2]))
end

function explicit(form::ExteriorDerivative{2, Primal})
	exprs = explicit(form.form)
	dq = ((exprs[2][1,0]-exprs[2][0,0])-(exprs[1][0,1]-exprs[1][0,0])) * Arrays.mskv
	dq.name = form.name

	return dq
end

#2-Forms
function explicit(form::FormVariable{2,P}) where {P}
	return Arrays.ArrayVariable(form.name)
end

function explicit(form::Addition{2,P}) where {P}
	return Arrays.Addition(form.name, explicit(form.left), explicit(form.right))
end

function explicit(form::InteriorProduct{1, Dual, Primal}) #TODO implement interpolations
	fexpr = explicit(form.form)
	uexpr, vexpr = explicit(form.vect)

	uout = -vexpr * (0.5*(fexpr[0,0]+fexpr[0,-1])) * Arrays.mskx
	vout = uexpr * (0.5*(fexpr[0,0]+fexpr[-1,0])) * Arrays.msky

	uout.name = form.name*"_x"
	vout.name = form.name*"_y"

	return (uout, vout)
end
