using SymbolicPhysics.Arrays

x = ArrayVariable("x")

flow = Sequence([Block(Dict("a" => x[1,0]-x[0,0], "b" => 3*x)), 
		 Block(Dict("c" =>2*x, "d" => x*x))])

println(string(flow))
