import ampere,grid,matrices
from numpy import linalg as la
from numpy import dot
import sys
import time

if __name__ == "__main__":
    amp = ampere.ampere('data/AMPERE_north_2015-10-01T16-00-00Z.txt')
    amp_grid = grid.generate_from_data(amp.get_grid())
    amp_jr = amp.get_data()

    g1=grid.uniform(51,25)
   
    t0=time.clock()
    mats = matrices.matrices(amp_grid,g1)
    print time.clock()-t0

    (L1,L2,divJ) = mats.L(8.,amp_jr)
    L2inv = la.pinv(L2)
    jperp=dot(mats.M1-dot(mats.M2,dot(L2inv,L1)),divJ)

    jp = jperp[1::2]
    jt = jperp[::2]

#    figure();pcolormesh(g1.p,g1.t,jt.reshape(g1.nt,g1.np));colorbar()
#    figure();pcolormesh(g1.p,g1.t,jp.reshape(g1.nt,g1.np));colorbar()

