import ampere,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos

if __name__ == "__main__":
    amp = ampere.ampere('data/AMPERE_north_2015-10-01T16-00-00Z.txt')
    amp_grid = grid.generate_from_data(amp.get_grid())
    amp_jr = amp.get_data()

    g1=grid.uniform(10,11)

    j = np.zeros((3,g1.t.shape[0],g1.t.shape[1]))

    t0 = time.clock()
    for iv in np.arange(g1.nt):
        for jv in np.arange(g1.np):
            t = g1.t[iv,jv]
            p = g1.p[iv,jv]

            for iu in np.arange(amp_grid.nt):
                for ju in np.arange(amp_grid.np):

                    ts = amp_grid.t[iu,ju]
                    ps = amp_grid.p[iu,ju]

                    tp = arccos( cos(t)*cos(ts) + sin(t)*sin(ts)*cos(p-ps) )
                    theta_prime_unit_t = 1./sin(tp)*(sin(ts)*cos(t)*cos(ps-p)-cos(ts)*sin(t))
                    phi_prime_unit_t = 1./sin(tp)*sin(ts)*sin(ps-p)

                    # print(theta_prime_unit_t,
                    #       phi_prime_unit_t,
                    #       theta_prime_unit_t**2+phi_prime_unit_t**2)

print time.clock()-t0
