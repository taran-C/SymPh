using Documenter, SymPh

pdf = false
local_html = false

if pdf
	format = Documenter.LaTeX(platform="tectonic")
else 
	format = Documenter.HTML(prettyurls = !local_html)
end

makedocs(sitename="SymPh.jl",
	format = format,
	pages = [
		 "Home" => "index.md",
		 "Manual" => Any[
			"Building an equation" => "man/equation.md",
				 ],
		 "Reference" => Any[
			"Differential Geometry Operators" => "man/operators.md",
				    ]
		 ]
	)

if !local_html
	deploydocs(
    		repo = "github.com/taran-C/SymPh.git",
	)
end
