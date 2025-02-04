"""
Interpolations TODO add weno and configurable interpolations
"""
interpx(a::Expression) = 0.5*(a[1,0]+a[0,0])
interpy(a::Expression) = 0.5*(a[0,1]+a[0,0])
interpdiag(a::Expression) = interpx(interpy(a))

up3(qmm::Expression, qm::Expression, qp::Expression) = (5*qm+2*qp-qmm)/6
up5(qmmm::Expression, qmm::Expression, qm::Expression, qp::Expression, qpp::Expression) = (2*qmmm - 13*qmm + 47*qm + 27*qp - 3*qpp)/60

#TODO better handle different direction
flx1x(U::Expression, a::Expression) = TernaryOperator(U>0, a[0,0], a[1,0])
flxup3x(U::Expression, a::Expression) = TernaryOperator(U>0, up3(a[-1,0], a[0,0], a[1,0]), up3(a[2,0], a[1,0], a[0,0]))
flxup5x(U::Expression, a::Expression) = TernaryOperator(U>0, up5(a[-2,0], a[-1,0], a[0,0], a[1,0], a[2,0]), up5(a[3,0], a[2,0], a[1,0], a[0,0], a[-1,0]))

flx1y(U::Expression, a::Expression) = TernaryOperator(U>0, a[0,0], a[0,1])
flxup3y(U::Expression, a::Expression) = TernaryOperator(U>0, up3(a[0,-1], a[0,0], a[0,1]), up3(a[0,2], a[0,1], a[0,0]))
flxup5y(U::Expression, a::Expression) = TernaryOperator(U>0, up5(a[0,-2], a[0,-1], a[0,0], a[0,1], a[0,2]), up5(a[0,3], a[0,2], a[0,1], a[0,0], a[0,-1]))


#TODO Specific interpolations for forms (choosing the right orders, directions...)
upwindx(U::Expression, a::Expression, o::Expression) = TernaryOperator(o > 4, flxup5x(U, a), 
						       TernaryOperator(o > 2, flxup3x(U, a),
						       TernaryOperator(o > 0, flx1x(U, a), RealValue(0))))

upwindy(U::Expression, a::Expression, o::Expression) = TernaryOperator(o > 4, flxup5y(U, a), 
						       TernaryOperator(o > 2, flxup3y(U, a),
						       TernaryOperator(o > 0, flx1y(U, a), RealValue(0))))

#Averages :
export avg4pt
avg4pt(U::Expression, dx, dy) = 0.25 * (U[0,0] + U[dx,0] + U[0,dy] + U[dx,dy])
