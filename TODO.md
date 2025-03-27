# To get a functional code
- Generate multiple functions called with a loop manager, by manipulating julia Expr objects rather than strings
- Replace strings by symbols/objects (would allow to actually pass function)
- Ensure valid and unique names
- Check placing for the 4 needed poisson solvers
- Check every possibility for every operator

# For the internship
- Higher order differentiation
- GPU/CPU switch

# Would be cool
- MPI
- Automatic optimization
- Code organization + exports
- Reflatten Loops (generate everything using step ? Nope, would need us to know nx at definition time, not in the philosophy)
- Automatically compute nh (max relative (depx, depy) on all nodes)
- Fuse nodes with same name in tree automatically (lots of duplicates rn)
- eps (machine epsilon) object for weno ?
- ^ (and other) operator, more conditionals for ease of use ?
