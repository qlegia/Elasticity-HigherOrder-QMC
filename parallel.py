from mpi4py import MPI
import numpy

import resource
def using(point=""):
    usage=resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           '''%(point,usage[0],usage[1],
                (usage[2]*resource.getpagesize())/1000000.0 )

def quadrature_MPI(lattice, f, transf = None):
    """
    Evaluation of equal-weight quadrature rule for the given lattice and function f
    in parallel using MPI.
    lattice: Lattice point set.
    f: integrand function
    transf: point transformation to be applied to x in [0,1]
    """
    assert(len(lattice)>0)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    # determine range of indices
    blocksize = int(len(lattice) / size)
    diff = len(lattice)-size*blocksize
    imin = 0
    for i in xrange(rank):
        imin += blocksize
        if i < diff: imin += 1
    imax = imin+blocksize
    if rank < diff: imax+=1
    myrange = xrange(imin,imax)
    #print "rank:",rank,"range:",myrange
    # quadrature loop
    S = 0.
    if transf == None:
        for i in myrange:
            #print using("rank="+str(rank)+", i="+str(i))
            S += f(lattice[i])
    else:
        for i in myrange:
            #print using("rank="+str(rank)+", i="+str(i))
            S += f(transf(lattice[i]))
    #S /= len(lattice)
    #print " -> ", "rank:",rank,"S=",S

    # MPI reduce
    #comm.Barrier()
    sendbuf = numpy.array([S])
    recvbuf = numpy.array([0.])
    comm.Reduce(sendbuf, recvbuf, op=MPI.SUM, root=0)
    if rank == 0: return recvbuf[0]/len(lattice)
    else: return 0

