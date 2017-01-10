import ampere_pole,grid
import sys
import time
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array,roll,hstack,newaxis,isscalar,ones,vstack
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate
import time
import solver

if __name__ == "__main__":
    Sigma_const = 8. # scalar for the time being
    interp = True

    # AMPERE stuff
    amp = ampere_pole.ampere_pole('data/AMPERE_north_2015-08-15T20-50-00Z.txt')

    grd = grid.generate_from_data(amp.get_grid())
    jr = amp.get_data()
    if interp:
        ntheta = 50
        nphi = 360

        # Complete full circle with AMPERE data so that interpolated grid doesn't have a huge hole
        # note, however, that we still perform the solve on a cell-centered grid, i.e., in the case of interpolation
        # here we make sure that the resultant grid doesn't wrap around
        p_tmp = hstack((grd.p,2.*pi+grd.p[:,[0]]))
        t_tmp = hstack((grd.t,grd.t[:,[0]]))

        f = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],jr,kx=1,ky=1)  # note, assume uniform ampere grid here
        # this here defines cell vortices (and wraps around)
        new_t = linspace(t_tmp.min(),t_tmp.max(),ntheta)
        new_p = linspace(p_tmp.min(),p_tmp[0,:].max(),nphi+1)
        # now redefine as cell centers
        new_p = 0.5*(new_p[:-1]+new_p[1:])
        grd = grid.grid()
        grd.t,grd.p = meshgrid(new_t,new_p,indexing='ij')
        jr = f(grd.t[:,0],grd.p[0,:])

    jr[abs(jr)<=0.2]=0.
    ############

    Nt,Np = grd.t.shape
    
    s = solver.solver(grd)

    # define conductance as uniform Pedersen for now
    SigmaP = ones_like(jr)*Sigma_const

    # low lat BC values, set to zero for now
    s.setLowLatBC(zeros(Np))

    s.setTerms(SigmaP)

    t0=time.clock()

    # Note, in the second option below we define both RHS and the stencil matrix inside the same
    # loops thus they're both defined in the same function as opposed to
    # the numpy version 1. That has more flexibility as it allows
    # solving the same matrix with a different RHS.

    s.setStencilMatrixNp()  # 1. this is numpy - way faster. 
    s.setRHSNp(jr)
#    s.setStencilMatrixAndRHS(jr)   # 2. This is fortran-ready (explicit loops). 


    print('Time spend constructing matrix (s):',time.clock()-t0)

    # Now, do the solve
    t0 = time.clock()
    pot = s.solve()
    print('Time spend on solve (s):',time.clock()-t0)

    print "Potential min/max (kV): ",pot.min(),pot.max()

    # Plotting

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,pot.reshape(Nt,Np),51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,jr,45)
    colorbar().set_label('Current, microA/m^2')
    contour(grd.p+pi+pi/2,grd.t,pot.reshape(Nt,Np),13,colors='black')
    

