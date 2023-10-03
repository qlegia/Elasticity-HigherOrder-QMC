# Adapted from
# https://github.com/jorgensd/dolfinx-tutorial/commit/e3669f2fe760748cbe79ff26d6bec3dd5f7251c6#
# to run the file, 
#  conda activate fenicsx-env
#  python elasticity3D.py
import numpy as np
import ufl

from mpi4py import MPI
from petsc4py.PETSc import ScalarType

from pyQMC.lattice import PolynomialLattice
# QMC points
dim_s = 16   # number of truncated dimensions
m = 3    # N = 2^m points. Choose N=2^m with m = 20 for reference solution
# name of file containing generating vector
qmc_dir = "./pyQMC-master1"
filename = qmc_dir+"/staircase2d_spod_a2_C0.1_SPOD_2dstaircase_t0/SPOD_2dstaircase_t0_m"+str(m)+".json"
# lat is a list containing all the QMC quadrature points
lat = PolynomialLattice(filename, dim_s)

k1 = np.arange(1,dim_s)  #{1,2,3,4,5,...,s}
k2 = np.arange(1,2*dim_s,2)  # {1,3,5,7,9,...,2s+1}

# generating function 
# psi = 1+sum_{j=1}^s y_j (sin(k1*pi*x1)*sin(k2*pi*x2)/ (j^2)
#    y_j is the jth-component of QMC points stored in lat
def psi_str(yp):
    stpsi = '1.0'
    for j in range(dim_s-1): 
        bt_term = 'pow('+str(j+1)+',2)'
        #st_term = '+'+str(yp[j]-0.5)+'*sin('+str(k1[j])+'*M_PI*x[0])*sin('+str(k2[j])+'*M_PI*x[1])/'+bt_term
        st_term = '+'+str(yp[j]-0.5)+'*sin('+str(k1[j])+'*pi*x[0])*sin('+str(k2[j])+'*pi*x[1])/'+bt_term
        stpsi = stpsi + st_term
    print('psi=',stpsi)
    return stpsi

class MyPsi:
    def __init__(self,yp):
        self.yp = yp

    def eval(self,x) :
        psi = 1.0
        for j in range(dim_s-1): 
           term = ( (self.yp[j]-0.5)*np.sin(k1[j]*pi*x[0])*np.sin(k2[j]*pi*x[1]) )/(j+1)**2
           psi = psi + term
        return np.full(x.shape[1], psi)

from dolfinx import mesh, fem, plot, io



# Scaled variable
L = 1
W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma
# We then create the mesh, which will consist of hexahedral elements, along with the function space. We will use the convenience function `VectorFunctionSpace`. However, we also could have used `ufl`s functionality, creating a vector element `element = ufl.VectorElement("CG", mesh.ufl_cell(), 1)
# `, and intitializing the function space as `V = dolfinx.fem.FunctionSpace(mesh, element)`.

nx = ny = 10 
domain = mesh.create_unit_square(MPI.COMM_WORLD,nx,ny) 
V = fem.VectorFunctionSpace(domain, ("P", 2))
# -
# ## Boundary conditions
# As we would like to clamp the boundary at $x=0$, we do this by using a marker function, which locate the facets where $x$ is close to zero by machine prescision.
# +
def clamped_boundary(x):
    return np.isclose(x[0], 0)

fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)
u_D = np.array([0,0], dtype=ScalarType)
bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets), V)

# As we want the traction $T$ over the remaining boundary to be $0$, we create a `dolfinx.Constant`
T = fem.Constant(domain, ScalarType((0, 0)))
# We also want to specify the integration measure $\mathrm{d}s$, which should be the integral over the boundary of our domain. We do this by using `ufl`, and its built in integration measures
ds = ufl.Measure("ds", domain=domain)
# ## Variational formulation
# We are now ready to create our variational formulation in close to mathematical syntax, as for the previous problems.
# +

Ws = fem.FunctionSpace(domain, ("P", 2))
x = ufl.SpatialCoordinate(domain)
sum_u = fem.Constant(domain,ScalarType((0,0)))

from math import sin, pi
for qmc_p in lat:
      # set up lambda = 1 + sum_{j=1}^dim_s 
      psi = MyPsi(qmc_p)
      psi.yp = qmc_p
      lam_x = fem.Function(Ws)
      lam_x.interpolate(psi.eval)
      #
      for qmc_q in lat:
           mu = MyPsi(qmc_q)
           mu.yp = qmc_q
           mu_x =  fem.Function(Ws)
           mu_x.interpolate(mu.eval)

           # Define strain and stress
           def epsilon(u):
               return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
           def sigma(u):
               return lam_x * ufl.nabla_div(u) * ufl.Identity(u.geometric_dimension()) + 2*mu_x*epsilon(u)

           u = ufl.TrialFunction(V)
           v = ufl.TestFunction(V)
           f = fem.Constant(domain, ScalarType((0, -rho*g)))
           a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
           L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

           # -
           # ```{note}
           # Note that we used `nabla_grad` and optionally `nabla_div` for the variational formulation, as 
           # opposed to our previous usage of 
           # `div` and `grad`. This is because for scalar functions $\nabla u$ has a clear meaning
           # $\nabla u = \left(\frac{\partial u}{\partial x}, \frac{\partial u}{\partial y}, \frac{\partial u}{\partial z} \right)$.
           #
           # However, if $u$ is vector valued, the meaning is less clear. Some sources define $\nabla u$ as a matrix with 
           # the elements $\frac{\partial u_j}{\partial x_i}$, while other  sources prefer 
           # $\frac{\partial u_i}{\partial x_j}$. In DOLFINx `grad(u)` is defined as the amtrix with element 
           #$\frac{\partial u_i}{\partial x_j}$. However, as it is common in continuum mechanics to use the other definition, 
           #`ufl` supplies us with `nabla_grad` for this purpose.
           # ```
           #  
           # ## Solve the linear variational problem
           # As in the previous demos, we assemble the matrix and right hand side vector and use PETSc to solve our variational problem
           problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
           uh = problem.solve()

           # ## Stress computation
           # As soon as the displacement is computed, we can compute various stress measures. We will compute the von Mises 
           # stress defined as $\sigma_m=\sqrt{\frac{3}{2}s:s}$ where $s$ is the deviatoric stress tensor 
           # $s(u)=\sigma(u)-\frac{1}{3}\mathrm{tr}(\sigma(u))I$.
           s = sigma(uh) -1./3*ufl.tr(sigma(uh))*ufl.Identity(uh.geometric_dimension())
           von_Mises = ufl.sqrt(3./2*ufl.inner(s, s))
           # The `von_Mises` variable is now an expression that must be projected into an appropriate function space so that 
           # we can visualize it. As `uh` is a linear combination of first order piecewise continuous functions, the von Mises
           # stresses will be a cell-wise constant function.

           V_von_mises = fem.FunctionSpace(domain, ("P", 1))
           stress_expr = fem.Expression(von_Mises, V_von_mises.element.interpolation_points())
           stresses = fem.Function(V_von_mises)
           stresses.interpolate(stress_expr)
