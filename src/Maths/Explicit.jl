import ..Arrays
export explicit


function explicit(form::Form{0, P}) where {P}
	return Arrays.ArrayVariable(form.name)
end
