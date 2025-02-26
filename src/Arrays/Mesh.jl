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
		o2px, o2py, o2dx, o2dy = compute_orders(msk, msk2d)

		#Creating the mesh
		return new(nx, ny, nh,
			   xc, yc, 1, 1, 1, 1, 1, 1,
			   dx, dy, A, 
			   msk0p, msk0d, msk1px, msk1py, msk1dx, msk1dy, msk2p, msk2d,
			   o2px, o2py, o2dx, o2dy)
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

function compute_metric(nx, ny, nh, Lx, Ly, msk)
	dx = msk .* Lx ./ (nx-2*nh)
	dy = msk .* Ly ./ (ny-2*nh)
	A = dx .* dy

	return dx, dy, A
end

function compute_orders(msk, mskv)
	nx,ny = size(msk)

	o2px = zeros(nx,ny)
	o2py = zeros(nx,ny)
	o2dx = zeros(nx,ny)
	o2dy = zeros(nx,ny)

	get_order_right(msk, 1, o2px)
	get_order_right(msk, nx, o2py)

	#TODO different masks for dual grid
	get_order_left(mskv, 1, o2dx)
	get_order_left(mskv, nx, o2dy)

	return o2px, o2py, o2dx, o2dy
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

o2px = ArrayVariable("mesh.o2px")
o2py = ArrayVariable("mesh.o2py")
o2dx = ArrayVariable("mesh.o2dx")
o2dy = ArrayVariable("mesh.o2dy")
