import SymbolicPhysics: @Let, State
using SymbolicPhysics.Maths
import SymbolicPhysics.Arrays

#Defining our equation
@Let u = FormVariable{1, Dual}() #Transported velocity
@Let U = Sharp(u) # U = u#

@Let transp = ExteriorDerivative(InteriorProduct(U, u)) + InteriorProduct(U, ExteriorDerivative(u)) #L_U(u)
@Let p = InverseLaplacian(Codifferential(transp)) #

#Time derivative
@Let dtu = -transp - ExteriorDerivative(p) #dtu = -(transport term) - dp

#Defining the parameters needed to explicit
explparams = ExplicitParam(; interp = Arrays.weno)

#Generating the RHS
rhs! = to_kernel(dtu; explparams = explparams)

