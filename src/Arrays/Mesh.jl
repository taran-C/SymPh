struct Mesh
	#Grid
	nx::Int64
	ny::Int64
	
	nh::Int64

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
	
	function Mesh(nx, ny, nh, msk, Lx = 1, Ly = 1)
		#Locations
		xc, yc = compute_locations(nx,ny,nh,Lx,Ly)
		
		#Metric
		dx, dy, A = compute_metric(nx,ny, nh, Lx, Ly, msk)

		#Masks
		msk0p, msk0d, msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d = compute_msks(msk)
	
		#Orders
		o1px, o1py, o1dx, o1dy, o2px, o2py, o2dx, o2dy = compute_orders(msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d)

		#Creating the mesh
		return new(nx, ny, nh,
			   xc, yc, 1, 1, 1, 1, 1, 1,
			   dx, dy, A, 
			   msk0p, msk0d, msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d,
			   o1px, o1py, o1dx, o1dy, o2px, o2py, o2dx, o2dy)
	end
end

function compute_locations(nx, ny, nh, Lx, Ly)
	xc = zeros(nx,ny)
	yc = zeros(nx,ny)
	for i in 1:nx, j in 1:ny
		xc[i,j] = (i-0.5-nh) * Lx/(nx-2*nh)
		yc[i,j] = (j-0.5-nh) * Ly/(ny-2*nh)
	end
	
	return xc, yc
end

#TODO primal and dual metric
function compute_metric(nx, ny, nh, Lx, Ly, msk)
	dx = ones(nx,ny) .* Lx ./ (nx-2*nh)
	dy = ones(nx,ny) .* Ly ./ (ny-2*nh)
	A = dx .* dy

	return dx, dy, A
end

function compute_orders(msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d)
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
	get_order_left(msk1px, 1, o1px)
	get_order_left(msk1py, nx, o1py)

	get_order_left(msk2p, 1, o2px)
	get_order_left(msk2p, nx, o2py)

	#Dual
	get_order_right(msk1dx, 1, o1dx)
	get_order_right(msk1dy, nx, o1dy)

	get_order_right(msk2d, 1, o2dx)
	get_order_right(msk2d, nx, o2dy)

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
