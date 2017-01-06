import ampere,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate

if __name__ == "__main__":
    Sigma = 8. # scalar for the time being
    interp = True

#    amp = ampere.ampere('data/AMPERE_north_2015-10-01T16-00-00Z.txt')
    amp = ampere.ampere('data/AMPERE_north_2015-08-15T12-00-00Z.txt')

    amp_grid = grid.generate_from_data(amp.get_grid())
    amp_jr = amp.get_data()

    # make finer grid
    if interp:
        ntheta = 44
        nphi = 360
        f = interpolate.RectBivariateSpline(amp_grid.t[:,0],amp_grid.p[0,:],amp_jr,kx=1,ky=1)
        new_t = linspace(amp_grid.t[:,0].min(),amp_grid.t[:,0].max(),ntheta)
        new_p = linspace(amp_grid.p[0,:].min(),amp_grid.p[0,:].max(),nphi)
        amp_grid = grid.grid()
        amp_grid.t,amp_grid.p = meshgrid(new_t,new_p,indexing='ij')
        amp_jr = f(amp_grid.t[:,0],amp_grid.p[0,:])

    dt = amp_grid.t[1,0]-amp_grid.t[0,0]
    dp = amp_grid.p[0,1]-amp_grid.p[0,0]

    # could do without it, but being gracious to the area
    # defined for Ampere grid inside "grid", in which case its already defined by now
    if interp:   
        amp_grid.A = sin(amp_grid.t)*dt*dp


    g1t = amp_grid.t+dt/2.
    g1p = amp_grid.p+dp/2.
    g1=  grid.grid(); # grid.uniform(88,30)
    g1.t=g1t
    g1.p=g1p

    j = np.zeros((3,g1.t.shape[0],g1.t.shape[1]))


    t0 = time.clock()

    t = mat(g1.t.ravel()).T
    p = mat(g1.p.ravel()).T

    ts = mat(amp_grid.t.ravel())
    ps = mat(amp_grid.p.ravel())
    dA = amp_grid.A.ravel()
    jr   = amp_jr.ravel()

    tp = arccos( dot(cos(t),cos(ts)) + multiply(dot(sin(t),sin(ts)),cos(repeat(p,ps.size,1)-ps)) )
    theta_prime_unit_t = multiply(1./sin(tp),( multiply(dot(cos(t),sin(ts)),cos(repeat(p,ps.size,1)-ps))-dot(sin(t),cos(ts))))
    phi_prime_unit_t = multiply(1./sin(tp),multiply(sin(repeat(ts,t.size,0)),sin(repeat(p,ps.size,1)-ps)))
    # if jr is in microA/m^2, then this should be in A/m
    jp = dot(multiply(1./tan(tp/2.),phi_prime_unit_t),multiply(jr,dA).T)*6.5/4./pi
    jt = dot(multiply(1./tan(tp/2.),theta_prime_unit_t),multiply(jr,dA).T)*6.5/4./pi

    print 'Time elapsed (s): ',time.clock()-t0

    ep = jp.A.reshape(g1.t.shape)/Sigma*1.e3 # should be in mV/m
    et = jt.A.reshape(g1.t.shape)/Sigma*1.e3
    e = sqrt(ep**2+et**2)
    phi = -et.cumsum(axis=0)*dt*6.5

    # Plotting

    # fixup periodic boundary before plotting
    e_plot = append(e,e[:,[0]],axis=1)
    phi_plot = append(phi,phi[:,[0]],axis=1)
    p_plot = append(g1.p,g1.p[:,[0]]+2*pi,axis=1)
    t_plot = append(g1.t,g1.t[:,[0]],axis=1)
    t1_plot = append(amp_grid.t,amp_grid.t[:,[0]],axis=1)


    figure()
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)


#    pcolormesh(g1.p+pi+pi/2,g1.t,e)
    contourf(p_plot+pi+pi/2,t_plot,e_plot,51)
    colorbar().set_label('|E|, mV/m')
    contour(p_plot+pi+pi/2,t1_plot,phi_plot,17,colors='white')

    figure()
    subplot(111,polar=True)
    pcolormesh(p_plot+pi+pi/2,t1_plot,phi_plot)
    colorbar().set_label('Phi, kV')

    show()
    print 'Phi min/max (kV):', phi.min(),phi.max()
    print '|E| min/max (mV/m):', e.min(),e.max()
