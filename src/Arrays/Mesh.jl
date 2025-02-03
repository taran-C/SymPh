struct Mesh
	#Grid
	nx
	ny
	
	nh

	#Locations
	xc
	yc

	xy
	yy

	xx
	yx

	xv
	yv

	#Metric TODO different edge length for primal/dual grids (and primal/dual grids at all anyway)
	dx
	dy
	
	A

	#Masks
	mskc
	
	mskv
	
	mskx
	msky
	
	function Mesh(nx, ny, nh, msk, Lx = 1, Ly = 1)
		#Locations
		xc = zeros(nx,ny)
		yc = zeros(nx,ny)
		for i in 1:nx, j in 1:ny
			xc[i,j] = (i-0.5-nh) * Lx/(nx-2*nh)
			yc[i,j] = (j-0.5-nh) * Ly/(ny-2*nh)
		end
		
		#Masks
		mskx = zeros(nx, ny)
		msky = zeros(nx, ny)
		mskv = zeros(nx, ny)

		for i in 2:nx-1, j in 2:ny-1
			mskx[i, j] = msk[i+1, j] * msk[i, j]
			msky[i, j] = msk[i, j+1] * msk[i, j]
			mskv[i, j] = msk[i, j] * msk[i+1, j] * msk[i, j+1] * msk[i+1, j+1]
		end

		return new(nx, ny, nh,
			   xc, yc, 1, 1, 1, 1, 1, 1,
			   1, 1, 1, 
			   msk, mskv, mskx, msky)
	end
end

dx = ArrayVariable("dx")
dy = ArrayVariable("dy")
A = ArrayVariable("A")

msk = ArrayVariable("msk")
mskv = ArrayVariable("mskv")
mskx = ArrayVariable("mskx")
msky = ArrayVariable("msky")
