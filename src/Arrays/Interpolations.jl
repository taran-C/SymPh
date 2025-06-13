F = Float64 #TODO Incorporate into settings

"""
Interpolation Orders
"""
#TODO old code, make more compatible
getrangeleft(q, step) = 1+3*step:length(q)-2*step
getrangeright(q, step) = 1+2*step:length(q)-3*step

#to compute q[i-1/2]
stencil6left(q, i, step) = (q[i-3*step], q[i-2*step], q[i-step], q[i], q[i+step], q[i+2*step])
stencil4left(q, i, step) = (q[i-2*step], q[i-step], q[i], q[i+step])
stencil2left(q, i, step) = (q[i-step], q[i])
#to compute q[i+1/2]
stencil6right(q, i, step) = (q[i-2*step], q[i-step], q[i], q[i+step], q[i+2*step], q[i+3*step])
stencil4right(q, i, step) = (q[i-step], q[i], q[i+step], q[i+2*step])
stencil2right(q, i, step) = (q[i], q[i+step])

function get_order_left(msk, step, order)
    let irange = getrangeleft(msk, step)
        for i in irange
            s6 = sum(stencil6left(msk, i, step))
            s4 = sum(stencil4left(msk, i, step))
            s2 = sum(stencil2left(msk, i, step))
            order[i] = floor(s6/6)*6 + (1-floor(s6/6))*(floor(s4/4)*4 + (1-floor(s4/4))*floor(s2/2)*2)
        end
    end
end

function get_order_right(msk, step, order)
    let irange = getrangeright(msk, step)
        for i in irange
            s6 = sum(stencil6right(msk, i, step))
            s4 = sum(stencil4right(msk, i, step))
            s2 = sum(stencil2right(msk, i, step))
            order[i] = floor(s6/6)*6 + (1-floor(s6/6))*(floor(s4/4)*4 + (1-floor(s4/4))*floor(s2/2)*2)
        end
    end
end

"""
Interpolations
"""

#Upwind interpolation
up3(qmm::Expression, qm::Expression, qp::Expression) = (5*qm+2*qp-qmm)/6
up5(qmmm::Expression, qmm::Expression, qm::Expression, qp::Expression, qpp::Expression) = (2*qmmm - 13*qmm + 47*qm + 27*qp - 3*qpp)/60

function upinterp(U::Expression, a::Expression, lr::String, dir::String, o::Integer)
	@assert lr in ["left", "right"]
	@assert dir in ["i", "j"]
	if lr == "right" 
		if o==1
			if dir=="i"
				return TernaryOperator(U>0, a[0,0], a[1,0])
			elseif dir=="j"
				return TernaryOperator(U>0, a[0,0], a[0,1])
			end
		elseif o==3
			if dir=="i"
				return TernaryOperator(U>0, up3(a[-1,0], a[0,0], a[1,0]), up3(a[2,0], a[1,0], a[0,0]))
			elseif dir=="j"
				return TernaryOperator(U>0, up3(a[0,-1], a[0,0], a[0,1]), up3(a[0,2], a[0,1], a[0,0]))
			end
		elseif o==5
			if dir=="i"
				return TernaryOperator(U>0, up5(a[-2,0], a[-1,0], a[0,0], a[1,0], a[2,0]), up5(a[3,0], a[2,0], a[1,0], a[0,0], a[-1,0]))
			elseif dir=="j"
				return TernaryOperator(U>0, up5(a[0,-2], a[0,-1], a[0,0], a[0,1], a[0,2]), up5(a[0,3], a[0,2], a[0,1], a[0,0], a[0,-1]))
			end
		end
	elseif lr == "left"
		if o==1
                        if dir=="i"
                                return TernaryOperator(U>0, a[-1,0], a[0,0])
                        elseif dir=="j"
                                return TernaryOperator(U>0, a[0,-1], a[0,0])
                        end
                elseif o==3
                        if dir=="i"
                                return TernaryOperator(U>0, up3(a[-2,0], a[-1,0], a[0,0]), up3(a[1,0], a[0,0], a[-1,0]))
                        elseif dir=="j"
                                return TernaryOperator(U>0, up3(a[0,-2], a[0,-1], a[0,0]), up3(a[0,1], a[0,0], a[0,-1]))
                        end
                elseif o==5
                        if dir=="i"
                                return TernaryOperator(U>0, up5(a[-3,0], a[-2,0], a[-1,0], a[0,0], a[1,0]), up5(a[2,0], a[1,0], a[0,0], a[-1,0], a[-2,0]))
                        elseif dir=="j"
				return TernaryOperator(U>0, up5(a[0,-3], a[0,-2], a[0,-1], a[0,0], a[0,1]), up5(a[0,2], a[0,1], a[0,0], a[0,-1], a[0,-2]))
                        end
                end
	end
end


#TODO Specific interpolations for forms (choosing the right orders, directions...)
"""
	upwind(U::Expression, a::Expression, o::Expression, lr::String, dir::String)

Five point upwind interpolation

# Arguments
- U : The transportant velocity
- a : The object to upwind
- o : The order of interpolation at that point
- lr : If `lr==\"left\"`, compute `a[i-0.5]`, else `a[i+0.5]`
- dir : along `i` or `j`
"""
upwind(U::Expression, a::Expression, o::Expression, lr::String, dir::String) = TernaryOperator(o > 4, upinterp(U, a, lr, dir, 5), 
						       TernaryOperator(o > 2, upinterp(U, a, lr, dir, 3),
						       TernaryOperator(o > 0, upinterp(U, a, lr, dir, 1), RealValue(0.0))))

"""
	avg2pt(U::Expression, a::Expression, o::Expression, lr::String, dir::String)

Two point average interpolation

# Arguments
- U : The transportant velocity, unused
- a : The object to interpolate
- o : The order of interpolation at that point, unused
- lr : If `lr==\"left\"`, compute `a[i-0.5]`, else `a[i+0.5]`
- dir : along `i` or `j`
"""
function avg2pt(U::Expression, a::Expression, o::Expression, lr::String, dir::String)
	@assert lr in ["left", "right"]
	@assert dir in ["i", "j"]
	
	if lr == "right"
		if dir == "i"
			return 0.5 * (a[1,0]+a[0,0])
		elseif dir == "j"
			return 0.5 * (a[0,1]+a[0,0])
		end
	elseif lr == "left"
		if dir == "i"
			return 0.5 * (a[-1,0]+a[0,0])
		elseif dir == "j"
			return 0.5 * (a[0,-1]+a[0,0])
		end
	end
end

#Weno TODO BROKEN
"""
	weno(U::Expression, a::Expression, o::Expression, lr::String, dir::String)

Weno interpolation of fifth order

# Arguments
- U : The transportant velocity
- a : The object to upwind
- o : The order of interpolation at that point
- lr : If `lr==\"left\"`, compute `a[i-0.5]`, else `a[i+0.5]`
- dir : along `i` or `j`
"""
function weno(U::Expression, a::Expression, o::Expression, lr::String, dir::String)
	@assert lr in ["left", "right"]
	@assert dir in ["i", "j"]

	if lr == "left" #TODO USE STENCIL HERE
		if dir == "i"
			qmmm = a[-3,0]
			qmm = a[-2,0]
			qm = a[-1,0]
			qp = a
			qpp = a[1,0]
			qppp = a[2,0]
		else
			qmmm = a[0,-3]
			qmm = a[0,-2]
			qm = a[0,-1]
			qp = a
			qpp = a[0,1]
			qppp = a[0,2]
		end
	else
		if dir == "j"
			qmmm = a[-2,0]
			qmm = a[-1,0]
			qm = a
			qp = a[1,0]
			qpp = a[2,0]
			qppp = a[3,0]
		else
			qmmm = a[0,-2]
			qmm = a[0,-1]
			qm = a
			qp = a[0,1]
			qpp = a[0,2]
			qppp = a[0,3]
		end
	end

	return TernaryOperator(o>5,
			TernaryOperator(U>0, wen5(qmmm, qmm, qm, qp, qpp), wen5(qppp, qpp, qp, qm, qmm)),
			TernaryOperator(o>3,
				TernaryOperator(U>0, wen3(qmm, qm, qp), wen3(qpp, qp, qm)),
				TernaryOperator(o>1, 
					TernaryOperator(U>0, qm, qp),
					RealValue(0.0)
				)
			)
		)

end

function wen5(qmm, qm, q0, qp, qpp)
	    """
	    5-points non-linear left-biased stencil reconstruction

	    qmm----qm-----q0--x--qp----qpp

	    An improved weighted essentially non-oscillatory scheme for hyperbolic
	    conservation laws, Borges et al, Journal of Computational Physics 227 (2008)
	    """
	# factor 6 missing here
	qi1 = 2 * qmm - 7 * qm + 11 * q0
	qi2 = -qm + 5 * q0 + 2 * qp
	qi3 = 2 * q0 + 5 * qp - qpp

	k1, k2 = 13/12, 0.25
	beta1 = k1 * (qmm - 2 * qm + q0)^2 + k2 * (qmm - 4 * qm + 3 * q0)^2
	beta2 = k1 * (qm - 2 * q0 + qp)^2 + k2 * (qm - qp)^2
	beta3 = k1 * (q0 - 2 * qp + qpp)^2 + k2 * (3 * q0 - 4 * qp + qpp)^2

	tau5 = abs(beta1 - beta3)

	g1, g2, g3 = 0.1, 0.6, 0.3
	w1 = g1 * (1 + tau5 / (beta1 + eps(F)))
	w2 = g2 * (1 + tau5 / (beta2 + eps(F)))
	w3 = g3 * (1 + tau5 / (beta3 + eps(F)))

	# factor 6 is hidden below
	return (w1 * qi1 + w2 * qi2 + w3 * qi3) / (6 * (w1 + w2 + w3))
end

function wen3(qm, q0, qp)
	qi1 = -qm / 2 + 3 * q0 / 2
	qi2 = (q0 + qp) / 2

	beta1 = (q0 - qm)^2
	beta2 = (qp - q0)^2
	tau = abs(beta2 - beta1)

	g1, g2 = 1/3, 2/3
	w1 = g1 * (1 + tau / (beta1 + eps(F)))
	w2 = g2 * (1 + tau / (beta2 + eps(F)))

	return (w1 * qi1 + w2 * qi2) / (w1 + w2)	
end


#Averages :
export avg4pt
avg4pt(U::Expression, di, dj) = 0.25 * (U[0,0] + U[di,0] + U[0,dj] + U[di,dj])

interpi(a::Expression) = 0.5*(a[1,0]+a[0,0])
interpj(a::Expression) = 0.5*(a[0,1]+a[0,0])
interpdiag(a::Expression) = interpx(interpy(a))
