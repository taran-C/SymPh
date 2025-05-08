#TODO handle bc (more order arrays to check if on border ?)
function fvtofd2(q::Expression, h::Expression, dir::String)
	#Simple scaling by edge length
	return q / h
end

function fvtofd4(q::Expression, h::Expression, dir::String)
	@assert dir in ["x", "y"]
	#Interior
	if dir == "x"
		return ((13/12) * q[0,0] - (1/24) * (q[0,-1] + q[0,+1])) / h
	elseif dir == "y"
		return ((13/12) * q[0,0] - (1/24) * (q[-1,0] + q[+1,0])) / h
	end

        #BC
        #Q[j,1] = (1/24)*(23*q[j,1]+q[j,2])
        #Q[j,-2] = (1/24)*(23*q[j,-2]+q[j,-3])
        #Q[j,0] = Q[j,1]
        #Q[j,-1] = Q[j,-2]
end
