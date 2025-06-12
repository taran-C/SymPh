# To get a functional code
- Replace strings by symbols/objects (would allow to actually pass function)
- Ensure valid and unique names
- Check placing for the 4 needed poisson solvers + periodicity
- Check every possibility for every operator

# For the internship
- Automatic arbitrary order finite-diff + fvtofd
- Test convergence on polar mesh

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

# Plotting
- Handle NaN values without breaking
- Replacing colors instead of overwriting
- Variable size of arrows depending on intensity
- Ability to define colormaps + Centered around a value to avoid flickering
