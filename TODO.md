# To get a functional code
- Replace strings by symbols/objects (would allow to actually pass function)
- Ensure valid and unique names
- Check placing for the 4 needed poisson solvers
- Check every possibility for every operator
- Halo filling + periodic conditions
- Forcing terms

# For the internship
- Higher order differentiation
- Automatic arbitrary order finite-diff + fvtofd
- Create plots for any metric
- Automatically separate diag and rhs in model ?

# Would be cool
- MPI
- Automatic optimization
- Code organization + exports
- Reflatten Loops (generate everything using step ? Nope, would need us to know nx at definition time, not in the philosophy)
- Automatically compute nh (max relative (depx, depy) on all nodes)
- Fuse nodes with same name in tree automatically (lots of duplicates rn)
- eps (machine epsilon) object for weno ?
- ^ (and other) operator, more conditionals for ease of use ?
- Automatic simplification which probably requires to :
- Represent Expressions as a sum of products
- Figure out time
- Replace everything dx/dy by di/dj etc to clarify index coordinates (except in metric where we actually have dx/dy)
