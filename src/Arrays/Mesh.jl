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
	dx#::AbstractArray{Float64}
	dy#::AbstractArray{Float64}
	
	A#::AbstractArray{Float64}

	#Masks
	mskc::AbstractArray{Float64}
	
	mskv::AbstractArray{Float64}
	
	mskx::AbstractArray{Float64}
	msky::AbstractArray{Float64}

	#Orders
	o2px::AbstractArray{Float64} #order two, primal, along x
	o2py::AbstractArray{Float64}
	
	function Mesh(nx, ny, nh, msk, Lx = 1, Ly = 1)
		#Locations
		xc, yc = compute_locations(nx,ny,nh,Lx,Ly)
		
		#Masks
		mskx, msky, mskv = compute_msks(msk)
		
		o2px, o2py = compute_orders(msk)

		#Creating the mesh
		return new(nx, ny, nh,
			   xc, yc, 1, 1, 1, 1, 1, 1,
			   1, 1, 1, 
			   msk, mskv, mskx, msky,
			   o2px, o2py)
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

function compute_orders(msk)
	nx,ny = size(msk)

	o2px = zeros(nx,ny)
	o2py = zeros(nx,ny)

	get_order_left(msk, 1, o2px)
	get_order_left(msk, nx, o2py)

	return o2px, o2py
end

function compute_msks(msk)
	nx, ny = size(msk)

	mskx = zeros(nx, ny)
	msky = zeros(nx, ny)
	mskv = zeros(nx, ny)

	for i in 2:nx-1, j in 2:ny-1
		mskx[i, j] = msk[i+1, j] * msk[i, j]
		msky[i, j] = msk[i, j+1] * msk[i, j]
		mskv[i, j] = msk[i, j] * msk[i+1, j] * msk[i, j+1] * msk[i+1, j+1]
	end
	
	return mskx, msky, mskv
end

#The mesh. is terrible in this
dx = ArrayVariable("mesh.dx")
dy = ArrayVariable("mesh.dy")
A = ArrayVariable("mesh.A")

msk = ArrayVariable("mesh.msk")
mskv = ArrayVariable("mesh.mskv")
mskx = ArrayVariable("mesh.mskx")
msky = ArrayVariable("mesh.msky")

o2px = ArrayVariable("mesh.o2px")
o2py = ArrayVariable("mesh.o2py")
