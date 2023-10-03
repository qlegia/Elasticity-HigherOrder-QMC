# Adapted from
# https://github.com/jorgensd/dolfinx-tutorial/commit/e3669f2fe760748cbe79ff26d6bec3dd5f7251c6#
# to run the file, 
#  conda activate fenicsx-env
#  python3 sparse_qmc L k q 
#  python3 sparse_qmc 
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py.PETSc import ScalarType

import sys 
from math import sin, pi

# implement the sparse grid sum (6.6)
# N_1^{(j)} = 2^{L-k},  L-k=9,...,1
# N_2^{(k)} = 2^k,      k=1,...,9
#
# compute
#   est1 = 2^L Q_{N^{L-k,k}} = sum_{j=1}^{2^{L-k}} sum_{i=1}^{2^k} calF(p(j), q(i))
# in partion qq for qq = 1,2,3,4 for L=10,
# let cnt be the counter run over 2 loops of p(j) and q(i) 
#  partition 1 = { 0 < cnt <=256}  
#  partition 2 = { 256 < cnt <= 512}  
#  partition 3 = { 512 < cnt <= 512+256=768}  
#  partition 4 = { 768 < cnt <= 1024}  
#  
# each k will be run in parallel, see sparse_qmc.pbs 
LL = int(sys.argv[1]) # first parameter from the command line
kk = int(sys.argv[2]) # second parameter from the command line
qq = int(sys.argv[3]) # third parameter from the command line

from pyQMC.lattice import PolynomialLattice
# QMC points
dim_s = 256   # number of truncated dimensions
print('dims=',dim_s)
# name of file containing generating vector
qmc_dir = "./pyQMC-master1"
N1 = 2**(LL-kk)  # LL=10 so LL-kk=9,8,7,6,5,4,3,2,1
N2 = 2**kk #                    k=1,2,3,4,5,6,7,8,9
print('N1= ',N1, ' N2= ', N2)
filename = qmc_dir+"/staircase2d_spod_a2_C0.1_SPOD_2dstaircase_t0/SPOD_2dstaircase_t0_m"+str(LL)+".json"
# lat is a list containing all the QMC quadrature points
lat = PolynomialLattice(filename, dim_s)

k1 = np.arange(1,dim_s)  #{1,2,3,4,5,...,s}
k2 = np.arange(1,2*dim_s,2)  # {1,3,5,7,9,...,2s+1}
step = 256

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



# We then create the mesh, which will consist of hexahedral elements, along with the function space. We will use the convenience function `VectorFunctionSpace`. However, we also could have used `ufl`s functionality, creating a vector element `element = ufl.VectorElement("CG", mesh.ufl_cell(), 1)
# `, and intitializing the function space as `V = dolfinx.fem.FunctionSpace(mesh, element)`.

nx = ny = 128 
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


cnt = 0  # counter for both p and q
est1 = 0.0
cmin = 2**LL
cmax = 0
for n1 in range(2**kk):
      #print('n1 =', n1)
      qmc_p = lat.__getitem__(n1*2**(LL-kk))
      # set up lambda = 1 + sum_{j=1}^dim_s
      #print('p=',qmc_p)
      psi = MyPsi(qmc_p)
      psi.yp = qmc_p
      lam_x = fem.Function(Ws)
      lam_x.interpolate(psi.eval)
      #
      for n2 in range(2**(LL-kk)):
         #print('n2 =', n2)
         qmc_q = lat.__getitem__(n2*2**(kk))  
         cnt = cnt + 1
         if ((cnt > (qq-1)*step)&(cnt <= step*qq)):
             mu = MyPsi(qmc_q)
             mu.yp = qmc_q
             mu_x =  fem.Function(Ws)
             mu_x.interpolate(mu.eval)
             #print('q=',qmc_q) 
             # Define strain and stress
             def epsilon(u):
               return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
             def sigma(u):
               return lam_x * ufl.nabla_div(u) * ufl.Identity(u.geometric_dimension()) + 2*mu_x*epsilon(u)
             def rhs(x):
               values = np.zeros((2, x.shape[1]), dtype=ScalarType)
               values[0] = 2*x[0]+10.0
               values[1] = x[1]-3.0
               return values

             u = ufl.TrialFunction(V)
             v = ufl.TestFunction(V)
             #f = fem.Constant(domain, ScalarType((0, -rho*g)))
             f = fem.Function(V)
             f.interpolate(rhs)
             a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx
             L = ufl.dot(f, v) * ufl.dx + ufl.dot(T, v) * ds

             # ## Solve the linear variational problem
             # As in the previous demos, we assemble the matrix and right hand side vector and use PETSc to solve our variational problem
             problem = fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
             uh = problem.solve()
    
             # compute the linear functional applied to uh
             one = fem.Constant(domain, ScalarType((1.0,1.0))) 
             u_one = fem.form(ufl.dot(uh,one)*ufl.dx)
             calL_u = fem.assemble_scalar(u_one)  # integral the whole domain
             est1 = est1 + calL_u
             print('est1 = ', est1)
             if (cnt < cmin):
                  cmin = cnt 
             if (cnt > cmax):
                  cmax = cnt
print('est1= ',est1,' cmin=', cmin, ' cmax=',cmax)
fsav = 'est1_N1_'+str(N1)+'_N2_'+str(N2)+'_q_'+str(qq)
np.savez(fsav, myN1 = N1, myN2 = N2, myk = kk, myL = LL, myq = qq, myS1 = est1)
