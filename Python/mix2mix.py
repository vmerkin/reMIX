import mix
import grid
import sys
import time
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array,roll,hstack,newaxis,isscalar,ones,vstack,arcsin,arctan2,zeros_like
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate
import time
import solver

filename = '/Users/merkivg1/work/Events/August2010/data/Aug2010_mix_2010-08-04T00-00-00Z.hdf'
interp = False

if __name__ == "__main__":

    ########################
    # read MIX data in
    ########################
    (simtime,
     x,y,
     psi_n,psi_s,
     fac_n,fac_s,
     sigmap_n,sigmap_s,
     sigmah_n,sigmah_s,
     energy_n,energy_s,
     flux_n,flux_s) = mix.get_data(filename)
    
    # get angular coordinates
    theta = arcsin(sqrt(x**2+y**2))
    phi   = (arctan2(y,x)+2*pi)%(2.*pi)

    # remove the periodic overlap; note the pole is already included; note also the transposition to conform to the new solver convention
    theta = theta[0,:]
    phi   = phi[:-1,1]
    jr    = fac_n[:-1,:].T
    SigmaP = sigmap_n[:-1,:].T
    SigmaH = sigmah_n[:-1,:].T
    ########################

    # generate grid 
    grd = grid.grid()
    grd.t,grd.p = meshgrid(theta,phi,indexing='ij')

    if interp:
        ntheta = 100
        nphi = 720

        # Complete full circle with AMPERE data so that interpolated grid doesn't have a huge hole
        # note, however, that we still perform the solve on a cell-centered grid, i.e., in the case of interpolation
        # here we make sure that the resultant grid doesn't wrap around
        p_tmp = hstack((grd.p,2.*pi+grd.p[:,[0]]))
        t_tmp = hstack((grd.t,grd.t[:,[0]]))

        f = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],jr,kx=1,ky=1)  
        fsp = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],SigmaP,kx=1,ky=1)  
        fsh = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],SigmaH,kx=1,ky=1)  
        # this here defines cell vortices (and wraps around)
        new_t = linspace(t_tmp.min(),t_tmp.max(),ntheta)
        new_p = linspace(p_tmp.min(),p_tmp[0,:].max(),nphi+1)
        # now redefine as cell centers
        new_p = 0.5*(new_p[:-1]+new_p[1:])
        grd = grid.grid()
        grd.t,grd.p = meshgrid(new_t,new_p,indexing='ij')
        jr = f(grd.t[:,0],grd.p[0,:])
        SigmaP = fsp(grd.t[:,0],grd.p[0,:])
        SigmaH = fsh(grd.t[:,0],grd.p[0,:])
    ############

    Nt,Np = grd.t.shape
    
    s = solver.solver(grd)

    # low lat BC values, set to zero for now
    s.setLowLatBC(zeros(Np))

    s.setTerms(SigmaP,SigmaH)

    t0=time.clock()

    # Note, in the second option below we define both RHS and the stencil matrix inside the same
    # loops thus they're both defined in the same function as opposed to
    # the numpy version 1. That has more flexibility as it allows
    # solving the same matrix with a different RHS.

#    s.setStencilMatrixNp()  # 1. this is numpy - way faster. 
#    s.setRHSNp(jr)
    s.setStencilMatrixAndRHS(jr)   # 2. This is fortran-ready (explicit loops). 


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
    contourf(grd.p+pi/2,grd.t,pot.reshape(Nt,Np),51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,50.*pi/180.)
    contourf(grd.p+pi/2,grd.t,jr,45)
    colorbar().set_label('Current, microA/m^2')
    contour(grd.p+pi/2,grd.t,pot.reshape(Nt,Np),13,colors='black')
    

