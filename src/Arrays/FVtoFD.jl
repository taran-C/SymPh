#TODO handle bc (more order arrays to check if on border ?)
function fvtofd2(q::Expression, h::Expression, dir::String)
	#Simple scaling by edge length
	return q# /h
end

#TODO separate division by h outside so the right one is used when composing two directions (and not h[i+1,j] or whatnot)
#OR NOT ? IDK
function fvtofd4(q::Expression, msk::Expression, dir::String)
	@assert dir in ["i", "j"]
	#Interior
	if dir == "i"
		#=
		return TernaryOperator(msk[-1,0] > 0,
				TernaryOperator(msk[+1,0] > 0,
					(13/12) * q[0,0] - (1/24) * (q[-1,0] + q[+1,0]),
				       (1/24) * (23*q[0,0] + q[-1,0])),
			       (1/24) * (23*q[0,0] + q[+1,0]))
		=#
		return ((13/12) * q[0,0] - (1/24) * (q[-1,0] + q[+1,0]))
		
	elseif dir == "j"
		#=
		return TernaryOperator(msk[0,-1] > 0,
				TernaryOperator(msk[0,+1] > 0,
					(13/12) * q[0,0] - (1/24) * (q[0,-1] + q[0,+1]),
				       (1/24) * (23*q[0,0] + q[0,-1])),
			       (1/24) * (23*q[0,0] + q[0,+1]))
		=#
		return ((13/12) * q[0,0] - (1/24) * (q[0,-1] + q[0,+1]))
	end

        #TODO BC
        #Q[j,1] = (1/24)*(23*q[j,1]+q[j,2])
        #Q[j,-2] = (1/24)*(23*q[j,-2]+q[j,-3])
        #Q[j,0] = Q[j,1]
        #Q[j,-1] = Q[j,-2]
end
