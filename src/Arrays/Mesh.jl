struct Mesh
	#Grid
	ni::Int64
	nj::Int64
	
	nh::Int64

	mgr #LoopManager

	#Locations
	xc::AbstractArray{Float64}
	yc::AbstractArray{Float64}

	xy#::AbstractArray{Float64}
	yy#::AbstractArray{Float64}

	xx#::AbstractArray{Float64}
	yx#::AbstractArray{Float64}

	xv::AbstractArray{Float64}
	yv::AbstractArray{Float64}

	#Metric TODO different edge length for primal/dual grids (and primal/dual grids at all anjway)
	dx::AbstractArray{Float64}
	dy::AbstractArray{Float64}
	
	A::AbstractArray{Float64}

	#Masks
	msk0p::AbstractArray{Float64}
	msk0d::AbstractArray{Float64}
	
	msk1pi::AbstractArray{Float64}
	msk1pj::AbstractArray{Float64}
	msk1di::AbstractArray{Float64}
	msk1dj::AbstractArray{Float64}
	
	msk2p::AbstractArray{Float64}
	msk2d::AbstractArray{Float64}

	#Orders
	o1pi::AbstractArray{Float64} #order two, primal, x along x
	o1pj::AbstractArray{Float64}
	o1di::AbstractArray{Float64}
	o1dj::AbstractArray{Float64}
	
	o2pi::AbstractArray{Float64} #order two, primal, along x
	o2pj::AbstractArray{Float64}
	o2di::AbstractArray{Float64}
	o2dj::AbstractArray{Float64}
	
	iperio::Bool
	jperio::Bool
	
	function Mesh(ni, nj, nh, mgr, msk, xc, yc, xv, yv; iperio=false, jperio=false)
		#Metric
		dx, dy, A = compute_metric(ni,nj, xc, yc)

		#Masks
		msks = compute_msks(msk)
	
		#Orders
		ords = compute_orders(msks, iperio, jperio)

		#Creating the mesh
		return new(ni, nj, nh,
			   mgr,
			   xc, yc, 1, 1, 1, 1, xv, yv,
			   dx, dy, A, 
			   msks...,
			   ords...,
			   iperio, jperio)
	end
end

#TODO only compute xc/yc etc if called ?
function get_x(msh::Mesh; location = "c")
	if location == "c"
		return msh.xc
	end
end
function get_y(msh::Mesh; location = "c")
	if location == "c"
		return msh.yc
	end
end


#TODO primal and dual metric
function compute_metric(ni, nj, xc, yc)
	dx, dy = zeros(ni, nj), zeros(ni, nj)
	
	#TODO fix dx and dy (actually compute dx along dirs with norm)
	for i in 2:ni, j in 2:nj
		dx[i,j] = sqrt((xc[i, j]-xc[i-1,j])^2 + (yc[i,j]-yc[i-1,j])^2)
		dy[i,j] = sqrt((xc[i, j]-xc[i,j-1])^2 + (yc[i,j]-yc[i,j-1])^2)
	end
	
	dx[1, :] .= dx[end, :]
	dy[:, 1] .= dy[:, end]
	
	A = dx .* dy

	return dx, dy, A
end

function compute_orders(msks, iperio, jperio)
	_, _, msk1pi, msk1pj, msk1di, msk1dj, msk2p, msk2d = msks
	ni,nj = size(msk2p)

	o1pi = zeros(ni,nj)
	o1pj = zeros(ni,nj)
	o1di = zeros(ni,nj)
	o1dj = zeros(ni,nj)

	o2pi = zeros(ni,nj)
	o2pj = zeros(ni,nj)
	o2di = zeros(ni,nj)
	o2dj = zeros(ni,nj)

	#Primal
	get_order_left(msk1pi, 1, o1pi)
	get_order_left(msk1pj, ni, o1pj)
	get_order_left(msk2p, 1, o2pi)
	get_order_left(msk2p, ni, o2pj)

	#Dual
	get_order_right(msk1di, 1, o1di)
	get_order_right(msk1dj, ni, o1dj)
	get_order_right(msk2d, 1, o2di)
	get_order_right(msk2d, ni, o2dj)

	return (o1pi, o1pj, o1di, o1dj,
		o2pi, o2pj, o2di, o2dj)
end

function compute_msks(msk)
	ni, nj = size(msk)

	msk0p = zeros(ni, nj)
	msk0d = msk
	
	msk1pi = zeros(ni, nj)
	msk1pj = zeros(ni, nj)
	msk1di = zeros(ni, nj)
	msk1dj = zeros(ni, nj)

	msk2p = msk
	msk2d = zeros(ni, nj)

	for i in 2:ni-1, j in 2:nj-1
		#msk0p[i,j] = max(msk[i,j], msk[i-1,j], msk[i,j-1], msk[i-1,j-1])
		msk0p[i,j] = max(msk[i,j], msk[i+1,j], msk[i,j+1], msk[i+1,j+1])
		
		msk1pi[i,j] = max(msk[i,j], msk[i,j-1])
		msk1pj[i,j] = max(msk[i,j], msk[i-1,j])
		msk1di[i,j] = msk[i,j] * msk[i-1,j]
		msk1dj[i,j] = msk[i,j] * msk[i,j-1]
		
		msk1pi[i,j] = msk1dj[i,j]
		msk1pj[i,j] = msk1di[i,j]

		msk2d[i,j] = msk[i,j] * msk[i-1,j] * msk[i,j-1] * msk[i-1,j-1] #TODO only for free-slip (null vorticity) (apparently, a whole lot to check there)
		msk0p[i,j] = msk2d[i,j]
	end
	
	return msk0p, msk0d, msk1pi, msk1pj, msk1di, msk1dj, msk2p, msk2d
end

#The mesh. is terrible in this
dx = ArrayVariable("mesh.dx")
dy = ArrayVariable("mesh.dy")
A = ArrayVariable("mesh.A")

msk0p = ArrayVariable("mesh.msk0p")
msk0d = ArrayVariable("mesh.msk0d")
msk1pi = ArrayVariable("mesh.msk1pi")
msk1pj = ArrayVariable("mesh.msk1pj")
msk1di = ArrayVariable("mesh.msk1di")
msk1dj = ArrayVariable("mesh.msk1dj")
msk2p = ArrayVariable("mesh.msk2p")
msk2d = ArrayVariable("mesh.msk2d")

o1pi = ArrayVariable("mesh.o1pi")
o1pj = ArrayVariable("mesh.o1pj")
o1di = ArrayVariable("mesh.o1di")
o1dj = ArrayVariable("mesh.o1dj")

o2pi = ArrayVariable("mesh.o2pi")
o2pj = ArrayVariable("mesh.o2pj")
o2di = ArrayVariable("mesh.o2di")
o2dj = ArrayVariable("mesh.o2dj")


"""
	CartesianMesh(ni, nj, nh, mgr, Lx = 1, Ly = 1; xperio=false, yperio=false)

A rectangular mesh of extent `Lx`, `Ly`
"""
function CartesianMesh(ni, nj, nh, mgr, Lx = 1, Ly = 1; xperio=false, yperio=false)
		#Locations
		xc, yc, xv, yv = compute_locations_cartesian(ni,nj,nh,Lx,Ly)
		msk = zeros(ni, nj)

		if xperio & yperio
			msk .= 1
		elseif xperio
			msk[1:ni, nh+1:nj-nh] .= 1
		elseif yperio
			msk[nh+1:ni-nh, 1:nj] .= 1
		else
			msk[nh+1:ni-nh, nh+1:nj-nh] .= 1
		end
		
		return Mesh(ni, nj, nh, mgr, msk, xc, yc, xv, yv; iperio = xperio, jperio=yperio)
end

function compute_locations_cartesian(ni, nj, nh, Lx, Ly)
	xc = zeros(ni,nj)
	yc = zeros(ni,nj)
	xv = zeros(ni,nj)
	yv = zeros(ni,nj)
	for i in 1:ni, j in 1:nj
		xc[i,j] = (i-0.5-nh) * Lx/(ni-2*nh)
		yc[i,j] = (j-0.5-nh) * Ly/(nj-2*nh)
		xv[i,j] = (i-1-nh) * Lx/(ni-2*nh)
		yv[i,j] = (j-1-nh) * Ly/(nj-2*nh)
	end
	
	return xc, yc, xv, yv
end

"""
	PolarMesh(ni, nj, nh, mgr, msk, rin = 0.5, rout = 1.5)

An annulus mesh with inner radius `rin` and outer radius `rout`

The i/x direction is the radial, j/y orthoradial
"""
function PolarMesh(ni, nj, nh, mgr, rin = 0.5, rout = 1.5)
		#Locations
		xc, yc, xv, yv = compute_locations_polar(ni,nj,nh,rin, rout)
		msk = zeros(ni, nj)
		msk[nh+1:ni-nh, :] .= 1
		
		return Mesh(ni, nj, nh, mgr, msk, xc, yc, xv, yv; iperio = false, jperio=true)
end

polar2cart(r, theta) = (r * cos(theta), r * sin(theta))

function compute_locations_polar(ni, nj, nh, rin, rout)
	xc = zeros(ni,nj)
	yc = zeros(ni,nj)
	xv = zeros(ni,nj)
	yv = zeros(ni,nj)

	for i in 1:ni, j in 1:nj
		r = rin + (i-0.5-nh) * (rout-rin) / (ni-2*nh)
		theta = (j-0.5-nh) * 2*pi/(nj-2*nh)
		coords = polar2cart(r, theta)
		xc[i,j] = coords[1] 
		yc[i,j] = coords[2]

		r = rin + (i-nh) * (rout-rin) / (ni-2*nh)
		theta = (j-nh) * 2*pi/(nj-2*nh)
		coords = polar2cart(r, theta)
		xv[i,j] = coords[1] 
		yv[i,j] = coords[2]
	end
	
	return xc, yc, xv, yv
end


