import ampere_pole,grid,matrices
from numpy import linalg as la
import numpy as np
import sys
import time
import point,sphtrig
from numpy import sin,cos,arccos,mat,dot,repeat,multiply,tan,pi,sqrt,append,linspace,meshgrid,zeros,arange,ones_like,array
from matplotlib.pyplot import figure,subplot,pcolormesh,colorbar,contour,show,contourf,rgrids,ylim
from scipy import interpolate
import time

if __name__ == "__main__":
    Sigma_const = 8. # scalar for the time being
    interp = False

    # AMPERE stuff
    amp = ampere_pole.ampere_pole('data/AMPERE_north_2015-08-15T20-50-00Z.txt')

    grd = grid.generate_from_data(amp.get_grid())
    jr = amp.get_data()
    if interp:
        ntheta = 100
        nphi = 720
        f = interpolate.RectBivariateSpline(grd.t[:,0],grd.p[0,:],jr,kx=1,ky=1)
        new_t = linspace(grd.t[:,0].min(),grd.t[:,0].max(),ntheta)
        new_p = linspace(grd.p[0,:].min(),grd.p[0,:].max(),nphi)
        grd = grid.grid()
        grd.t,grd.p = meshgrid(new_t,new_p,indexing='ij')
        jr = f(grd.t[:,0],grd.p[0,:])

    jr[abs(jr)<=0.2]=0.

    dt = grd.t[1:,:]-grd.t[:-1,:]
    dp = grd.p[:,1:]-grd.p[:,:-1]

    (Nt,Np) = grd.t.shape
    ############

    K = lambda i,j: j*Nt+i
    A=zeros((Nt*Np,Nt*Np))
    J = zeros(Nt*Np)

    # inner block
    for i in arange(1,Nt-1):
        for j in arange(0,Np):
            jAm1 = (j-1)%Np
            jAp1 = (j+1)%Np
            jA = j

            jm1 = (j-1)%(Np-1)
            jp1  = (j+1)%(Np-1)
            j      = j%(Np-1)

            ft = 2./(dt[i,j]+dt[i-1,j])
            fp = 2./(dp[i,j]+dp[i,jm1])/sin(grd.t[i,j])**2

            A[K(i,jA),K(i,jA)] = -ft*(1/dt[i,j]+1/dt[i-1,j]) - fp*(1/dp[i,j]+1/dp[i,jm1])  - 0.5*ft*(dt[i,j]/dt[i-1,j]-dt[i-1,j]/dt[i,j])/tan(grd.t[i,j])
            A[K(i,jA),K(i+1,jA)] = ft*(1./dt[i,j]+0.5*dt[i-1,j]/dt[i,j]/tan(grd.t[i,j]))
            A[K(i,jA),K(i-1,jA)] = ft*(1./dt[i-1,j]-0.5*dt[i,j]/dt[i-1,j]/tan(grd.t[i,j]))
            A[K(i,jA),K(i,jAp1)] = fp/dp[i,j]
            A[K(i,jA),K(i,jAm1)] = fp/dp[i,jm1]

            J[K(i,jA)] = jr[i,j]

    # low lat boundary
    for j in arange(0,Np):
        A[K(Nt-1,j),K(Nt-1,j)] = 1.
        J[K(Nt-1,j)] = 0. # FIXME: make arbitrary BC

    # pole boundary
    for j in arange(0,Np):
        A[K(0,j),K(0,j)] = 1.
        for jj in arange(0,Np): 
            A[K(0,j),K(1,jj)] = -dp[1,jj%(Np-1)]/(2.*pi)
        # note, not setting J because it's initializaed to zero anyway

    # Now, do the solve
    import scipy
    from scipy.sparse import linalg
    pot,status=scipy.sparse.linalg.lgmres(A,J)
    pot*=6.5**2*1.e3/Sigma_const

    print "Convergence status: ", status
    print "Potential min/max (kV): ",pot.min(),pot.max()

    # Plotting

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,pot.reshape(Np,Nt).T,51)
    colorbar().set_label('Potential, kV')

    figure() 
    subplot(111,polar=True)
    circles = linspace(10,40,4)*pi/180.
    lbls = [str(elem)+u'\xb0' for elem in [10,20,30,40]] 
    lines, labels = rgrids(circles,lbls)
    ylim(0,44.*pi/180.)
    contourf(grd.p+pi+pi/2,grd.t,jr,51)
    colorbar().set_label('Current, microA/m^2')
    contour(grd.p+pi+pi/2,grd.t,pot.reshape(Np,Nt).T,13,colors='black')
    

