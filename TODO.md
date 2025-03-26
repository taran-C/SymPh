# To get a functional code
- eps (machine epsilon) object for weno ?
- ^ (and other) operator, more conditionals for ease of use ?
- Actually generate typed kernels with automatic argument passing
- Special ArrayExpression object representing function calls that need sync (poisson problem)

# For the internship
- Higher order differentiation
- GPU/CPU switch

# Would be cool
- MPI
- Automatic optimization
- Code organization + exports
- Ensure valid+unique names + Replace strings by symbols ?
- Reflatten Loops (generate everything using step ? Nope, would need us to know nx at definition time, not in the philosophy)
- Automatically compute nh (max relative (depx, depy) on all nodes)
- Fuse nodes with same name in tree automatically (lots of duplicates rn)
