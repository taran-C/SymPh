using Documenter, SymPh

#Parsing command-line arguments
if "pdf" in ARGS
	pdf = true
else
	pdf = false
end
if "local" in ARGS
	local_html = true
else
	local_html = false
end

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
			"Maths" => "ref/maths.md",
			"Arrays" => "ref/arrays.md",
			"Miscellaneous" => "ref/misc.md"
				    ]
		 ]
	)

if !local_html & !pdf
	deploydocs(
    		repo = "github.com/taran-C/SymPh.git",
	)
end
