"""
Interpolation Orders
"""
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
Interpolations TODO add weno and configurable interpolations
"""
interpx(a::Expression) = 0.5*(a[1,0]+a[0,0])
interpy(a::Expression) = 0.5*(a[0,1]+a[0,0])
interpdiag(a::Expression) = interpx(interpy(a))

up3(qmm::Expression, qm::Expression, qp::Expression) = (5*qm+2*qp-qmm)/6
up5(qmmm::Expression, qmm::Expression, qm::Expression, qp::Expression, qpp::Expression) = (2*qmmm - 13*qmm + 47*qm + 27*qp - 3*qpp)/60

function flx(U::Expression, a::Expression, lr::String, dir::String, o::Integer)
	@assert lr in ["left", "right"]
	@assert dir in ["x", "y"]
	if lr == "right" 
		if o==1
			if dir=="x"
				return TernaryOperator(U>0, a[0,0], a[1,0])
			elseif dir=="y"
				return TernaryOperator(U>0, a[0,0], a[0,1])
			end
		elseif o==3
			if dir=="x"
				return TernaryOperator(U>0, up3(a[-1,0], a[0,0], a[1,0]), up3(a[2,0], a[1,0], a[0,0]))
			elseif dir=="y"
				return TernaryOperator(U>0, up3(a[0,-1], a[0,0], a[0,1]), up3(a[0,2], a[0,1], a[0,0]))
			end
		elseif o==5
			if dir=="x"
				return TernaryOperator(U>0, up5(a[-2,0], a[-1,0], a[0,0], a[1,0], a[2,0]), up5(a[3,0], a[2,0], a[1,0], a[0,0], a[-1,0]))
			elseif dir=="y"
				return TernaryOperator(U>0, up5(a[0,-2], a[0,-1], a[0,0], a[0,1], a[0,2]), up5(a[0,3], a[0,2], a[0,1], a[0,0], a[0,-1]))
			end
		end
	elseif lr == "left"
		if o==1
                        if dir=="x"
                                return TernaryOperator(U>0, a[-1,0], a[0,0])
                        elseif dir=="y"
                                return TernaryOperator(U>0, a[0,-1], a[0,0])
                        end
                elseif o==3
                        if dir=="x"
                                return TernaryOperator(U>0, up3(a[-2,0], a[-1,0], a[0,0]), up3(a[1,0], a[0,0], a[-1,0]))
                        elseif dir=="y"
                                return TernaryOperator(U>0, up3(a[0,-2], a[0,-1], a[0,0]), up3(a[0,1], a[0,0], a[0,-1]))
                        end
                elseif o==5
                        if dir=="x"
                                return TernaryOperator(U>0, up5(a[-3,0], a[-2,0], a[-1,0], a[0,0], a[1,0]), up5(a[2,0], a[1,0], a[0,0], a[-1,0], a[-2,0]))
                        elseif dir=="y"
				return TernaryOperator(U>0, up5(a[0,-3], a[0,-2], a[0,-1], a[0,0], a[0,1]), up5(a[0,2], a[0,1], a[0,0], a[0,-1], a[0,-2]))
                        end
                end
	end
end


#TODO Specific interpolations for forms (choosing the right orders, directions...)
upwind(U::Expression, a::Expression, o::Expression, lr::String, dir::String) = TernaryOperator(o > 4, flx(U, a, lr, dir, 5), 
						       TernaryOperator(o > 2, flx(U, a, lr, dir, 3),
						       TernaryOperator(o > 0, flx(U, a, lr, dir, 1), RealValue(0))))

#Averages :
export avg4pt
avg4pt(U::Expression, dx, dy) = 0.25 * (U[0,0] + U[dx,0] + U[0,dy] + U[dx,dy])
