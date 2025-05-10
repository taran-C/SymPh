struct Mesh
	#Grid
	nx::Int64
	ny::Int64
	
	nh::Int64

	mgr #LoopManager

	#Locations
	xc::AbstractArray{Float64}
	yc::AbstractArray{Float64}

	xy#::AbstractArray{Float64}
	yy#::AbstractArray{Float64}

	xx#::AbstractArray{Float64}
	yx#::AbstractArray{Float64}

	xv#::AbstractArray{Float64}
	yv#::AbstractArray{Float64}

	#Metric TODO different edge length for primal/dual grids (and primal/dual grids at all anyway)
	dx::AbstractArray{Float64}
	dy::AbstractArray{Float64}
	
	A::AbstractArray{Float64}

	#Masks
	msk0p::AbstractArray{Float64}
	msk0d::AbstractArray{Float64}
	
	msk1px::AbstractArray{Float64}
	msk1py::AbstractArray{Float64}
	msk1dx::AbstractArray{Float64}
	msk1dy::AbstractArray{Float64}
	
	msk2p::AbstractArray{Float64}
	msk2d::AbstractArray{Float64}

	#Orders
	o1px::AbstractArray{Float64} #order two, primal, x along x
	o1py::AbstractArray{Float64}
	o1dx::AbstractArray{Float64}
	o1dy::AbstractArray{Float64}
	
	o2px::AbstractArray{Float64} #order two, primal, along x
	o2py::AbstractArray{Float64}
	o2dx::AbstractArray{Float64}
	o2dy::AbstractArray{Float64}
	
	xperio::Bool
	yperio::Bool
	
	function Mesh(nx, ny, nh, mgr, msk, xc, yc; xperio=false, yperio=false)
		#Metric
		dx, dy, A = compute_metric(nx,ny, xc, yc)

		#Masks
		msk0p, msk0d, msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d = compute_msks(msk)
	
		#Orders
		o1px, o1py, o1dx, o1dy, o2px, o2py, o2dx, o2dy = compute_orders(msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d, xperio, yperio)

		#Creating the mesh
		return new(nx, ny, nh,
			   mgr,
			   xc, yc, 1, 1, 1, 1, 1, 1,
			   dx, dy, A, 
			   msk0p, msk0d, msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d,
			   o1px, o1py, o1dx, o1dy, o2px, o2py, o2dx, o2dy,
			   xperio, yperio)
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
function compute_metric(nx, ny, xc, yc)
	dx, dy = zeros(nx, ny), zeros(nx, ny)
	
	#TODO fix dx and dy (actually compute dx along dirs with norm)
	for i in 2:nx, j in 2:ny
		dx[i,j] = sqrt((xc[i, j]-xc[i-1,j])^2 + (yc[i,j]-yc[i-1,j])^2)
		dy[i,j] = sqrt((xc[i, j]-xc[i,j-1])^2 + (yc[i,j]-yc[i,j-1])^2)
	end
	#dx[2:end, 1:end] .= xc[2:end, 1:end] .- xc[1:end-1, 1:end]
	#dx[1, :] .= dx[end, :]

	#dy[1:end, 2:end] .= yc[1:end, 2:end] .- yc[1:end, 1:end-1]
	#dy[:, 1] .= dy[:, end]
	
	A = dx .* dy

	return dx, dy, A
end

function compute_orders(msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d, xperio, yperio)
	nx,ny = size(msk2p)

	o1px = zeros(nx,ny)
	o1py = zeros(nx,ny)
	o1dx = zeros(nx,ny)
	o1dy = zeros(nx,ny)

	o2px = zeros(nx,ny)
	o2py = zeros(nx,ny)
	o2dx = zeros(nx,ny)
	o2dy = zeros(nx,ny)

	#Primal
	if xperio
		o1px .= 6
	else
		get_order_left(msk1px, 1, o1px)
	end
	if yperio
		o1py .= 6
	else
		get_order_left(msk1py, nx, o1py)
	end

	if xperio
		o2px .= 6
	else
		get_order_left(msk2p, 1, o2px)
	end
	if yperio
		o2py .= 6
	else
		get_order_left(msk2p, nx, o2py)
	end

	#Dual
	if xperio
		o1dx .= 6
	else
		get_order_right(msk1dx, 1, o1dx)
	end
	if yperio
		o1dy .= 6
	else
		get_order_right(msk1dy, nx, o1dy)
	end

	if xperio
		o2dx .= 6
	else
		get_order_right(msk2d, 1, o2dx)
	end
	if yperio
		o2dy .= 6
	else
		get_order_right(msk2d, nx, o2dy)
	end

	return (o1px, o1py, o1dx, o1dy,
		o2px, o2py, o2dx, o2dy)
end

function compute_msks(msk)
	nx, ny = size(msk)

	msk0p = zeros(nx, ny)
	msk0d = msk
	
	msk1px = zeros(nx, ny)
	msk1py = zeros(nx, ny)
	msk1dx = zeros(nx, ny)
	msk1dy = zeros(nx, ny)

	msk2p = msk
	msk2d = zeros(nx, ny)

	for i in 2:nx-1, j in 2:ny-1
		msk0p[i,j] = max(msk[i,j], msk[i-1,j], msk[i,j-1], msk[i-1,j-1])
		
		msk1px[i,j] = max(msk[i,j], msk[i,j-1])
		msk1py[i,j] = max(msk[i,j], msk[i-1,j])
		msk1dx[i,j] = msk[i,j] * msk[i-1,j]
		msk1dy[i,j] = msk[i,j] * msk[i,j-1]

		msk2d[i,j] = msk[i,j] * msk[i-1,j] * msk[i,j-1] * msk[i-1,j-1] #TODO only for free-slip (null vorticity) (apparently, a whole lot to check there)
	end
	
	return msk0p, msk0d, msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d
end

#The mesh. is terrible in this
dx = ArrayVariable("mesh.dx")
dy = ArrayVariable("mesh.dy")
A = ArrayVariable("mesh.A")

msk0p = ArrayVariable("mesh.msk0p")
msk0d = ArrayVariable("mesh.msk0d")
msk1px = ArrayVariable("mesh.msk1px")
msk1py = ArrayVariable("mesh.msk1py")
msk1dx = ArrayVariable("mesh.msk1dx")
msk1dy = ArrayVariable("mesh.msk1dy")
msk2p = ArrayVariable("mesh.msk2p")
msk2d = ArrayVariable("mesh.msk2d")

o1px = ArrayVariable("mesh.o1px")
o1py = ArrayVariable("mesh.o1py")
o1dx = ArrayVariable("mesh.o1dx")
o1dy = ArrayVariable("mesh.o1dy")

o2px = ArrayVariable("mesh.o2px")
o2py = ArrayVariable("mesh.o2py")
o2dx = ArrayVariable("mesh.o2dx")
o2dy = ArrayVariable("mesh.o2dy")


"""
	CartesianMesh(nx, ny, nh, mgr, msk, Lx = 1, Ly = 1; xperio=false, yperio=false)

A rectangular mesh of extent `Lx`, `Ly`
"""
function CartesianMesh(nx, ny, nh, mgr, msk, Lx = 1, Ly = 1; xperio=false, yperio=false)
		#Locations
		xc, yc = compute_locations_cartesian(nx,ny,nh,Lx,Ly)
		
		return Mesh(nx, ny, nh, mgr, msk, xc, yc; xperio = xperio, yperio=yperio)
end

function compute_locations_cartesian(nx, ny, nh, Lx, Ly)
	xc = zeros(nx,ny)
	yc = zeros(nx,ny)
	for i in 1:nx, j in 1:ny
		xc[i,j] = (i-0.5-nh) * Lx/(nx-2*nh)
		yc[i,j] = (j-0.5-nh) * Ly/(ny-2*nh)
	end
	
	return xc, yc
end

"""
	PolarMesh(nx, ny, nh, mgr, msk, rin = 0.5, rout = 1.5)

An annulus mesh with inner radius `rin` and outer radius `rout`

The i/x direction is the radial, j/y orthoradial
"""
function PolarMesh(nx, ny, nh, mgr, msk, rin = 0.5, rout = 1.5)
		#Locations
		xc, yc = compute_locations_polar(nx,ny,nh,rin, rout)
		
		return Mesh(nx, ny, nh, mgr, msk, xc, yc; xperio = false, yperio=true)
end

polar2cart(r, theta) = (r * cos(theta), r * sin(theta))

function compute_locations_polar(nx, ny, nh, rin, rout)
	xc = zeros(nx,ny)
	yc = zeros(nx,ny)

	for i in 1:nx, j in 1:ny
		r = rin + (i-0.5-nh) * (rout-rin) / (nx-2*nh)
		theta = (j-0.5-nh) * 2*pi/(ny-2*nh)
		coords = polar2cart(r, theta)
		xc[i,j] = coords[1] 
		yc[i,j] = coords[2]
	end
	
	return xc, yc
end


